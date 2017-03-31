#!/usr/bin/python

import os
import sys
import math
import csv
import json
import re
import errno
import shutil
import subprocess
import time
import mirtk
from datetime import datetime


##############################################################################
# Spatio-temporal atlas

class SpatioTemporalAtlas(object):
    """Spatio-temporal image and deformation atlas

    Configuration entries such as file and directory paths, atlas time
    points (means) and corresponding temporal kernel width/standard deviation
    (sigma) must be given as dictionary, usually read from JSON file.

    The `construct` function performs the atlas construction steps.
    Template images for each atlas time point specified in the configuration
    are written to the configured `outdir`. At intermediate steps, an average
    image is created for each time point associated with an image in the
    dataset from which the atlas is created. These can be found in the `img`
    subdirectory of each iteration underneath the configured `tmpdir`.

    After atlas construction, average images and longitudinal deformations
    for every continuous time (within the atlas precision/max. resolution)
    can be computed as long as the required deformations are stored.

    """

    def _path(self, paths, name, default):
        """Get absolute path from config dictionary."""
        path = paths.get(name, default)
        return os.path.normpath(os.path.join(self.topdir, path)) if path else None

    def __init__(self, config, step=-1, verbose=1, threads=0):
        """Load spatio-temporal atlas configuration."""
        self.step = step
        self.config = os.path.abspath(config)
        with open(self.config, "rt") as f:
            cfg = json.load(f)
        # Paths of image file and output of prior global normalization step
        paths = cfg.get("paths", {})
        cfgdir = os.path.dirname(self.config)
        self.topdir = os.path.normpath(os.path.join(cfgdir, paths.get("topdir", ".")))
        self.agecsv = self._path(paths, "agecsv", os.path.join(cfgdir, "ages.csv"))
        self.imgcsv = self._path(paths, "imgcsv", os.path.join(cfgdir, "subjects.csv"))
        self.imgage = read_ages(self.agecsv)
        self.imgids = read_imgids(self.imgcsv)
        self.imgdir = self._path(paths, "imgdir", "images")
        self.imgpre = paths.get("imgpre", "")
        self.imgsuf = paths.get("imgsuf", ".nii.gz")
        self.affdir = self._path(paths, "affdir", None)
        self.affpre = paths.get("affpre", "")
        self.affsuf = paths.get("affsuf", ".dof")
        self.tmpdir = self._path(paths, "tmpdir", "cache")
        self.outdir = self._path(paths, "outdir", os.path.join(self.topdir, "templates"))
        self.refimg = self._path(paths, "refimg", os.path.join(self.outdir, "linavg.nii.gz"))
        # Parameters of temporal kernel regressions
        regression = cfg.get("regression", {})
        self.means = [float(t) for t in regression['means']]
        self.sigma = regression.get("sigma", 1.)
        if isinstance(self.sigma, (tuple, list)):
            self.sigma = [float(t) for t in self.sigma]
            if len(self.sigma) != len(self.means):
                raise ValueError("Number of sigma values must equal number of mean values!")
        else:
            self.sigma = [float(self.sigma)] * len(self.means)
        self.epsilon = float(regression.get("epsilon", .001))
        self.precision = int(regression.get('precision', 2))
        # Registration parameters
        registration = cfg.get("registration")
        self.num_bch_terms = 3  # composition should be symmetric, e.g., 2, 3, or 5 terms
        self.residual_svffds = registration.get("residual_svffds", True)
        self.age_specific_svffds = registration.get("age_specific_svffds", True)
        self.age_specific_imgdof = registration.get("age_specific_imgdof", True)
        self.reg_config = self._path(registration, "config", os.path.join(cfgdir, "register.cfg"))
        # Other workflow execution settings
        self.verbose = verbose
        self.threads = threads
        self.condor = cfg.get("condor", {})

    def _run(self, command, args=[], opts={}):
        """Execute MIRTK command."""
        if command in ("edit-dofs"):
            # These commands do not support the -threads option (yet)
            threads = 0
        else:
            threads = self.threads
            if "verbose" not in opts:
                opts["verbose"] = self.verbose - 2
        mirtk.run(command, args=args, opts=opts, verbose=self.verbose - 1, threads=threads)

    def _submit(self, name, script, tasks=-1, opts={}, step=-1, queue=None):
        """Submit batch script."""
        if step < 0:
            step = self.step
        if not queue:
            queue = self.queue
        if tasks == 0:
            return (queue, 0)
        source = 'import sys\n'
        source += 'sys.path.insert(0, "{0}")\n'.format(os.path.dirname(os.path.dirname(mirtk.__file__)))
        source += 'sys.path.insert(0, "{0}")\n'.format(os.path.dirname(__file__))
        source += 'from {0} import SpatioTemporalAtlas\n'.format(os.path.splitext(os.path.basename(__file__))[0])
        source += 'atlas = SpatioTemporalAtlas("{config}", step={step}, threads={threads}, verbose=3)\n'
        if tasks > 0:
            source += 'taskid = int(sys.argv[1])\n'
            source += 'if taskid < 0: sys.stderr.write("Invalid task ID\\n")\n'
        source += script
        source += 'else: sys.stderr.write("Invalid task ID\\n")\n'
        opts.update({
            "config": self.config,
            "step": step,
            "threads": self.threads
        })
        jobname = "i{0:02d}_{1}".format(step, name) if step >= 0 else name
        if queue.lower() in ("condor", "htcondor"):
            if tasks > 0:
                log = os.path.join(self.subdir(step), "log", name + "_$(Cluster).$(Process).log")
            else:
                log = os.path.join(self.subdir(step), "log", name + "_$(Cluster).log")
            jobid = cbatch(name=jobname, script=source, tasks=tasks, opts=opts, log=log,
                           requirements=self.condor.get("requirements", []), verbose=0)
        else:
            if tasks > 0:
                log = os.path.join(self.subdir(step), "log", name + "_%A.%a.log")
            else:
                log = os.path.join(self.subdir(step), "log", name + "_%j.log")
            jobid = sbatch(name=jobname, script=source, tasks=tasks, opts=opts, log=log,
                           threads=self.threads, queue=queue, verbose=0)
        if tasks > 0:
            self.info("Submitted job '{}' (id={}, #jobs={})".format(
                name, jobid[0] if isinstance(jobid, tuple) else jobid, tasks)
            )
        else:
            self.info("Submitted job '{}' (id={})".format(name, jobid))
        return (queue, jobid)

    def wait(self, jobs, interval=60, verbose=5):
        """Wait for batch jobs to complete."""
        if not isinstance(jobs, list):
            jobs = [jobs]
        condor_jobs = []
        slurm_jobs = []
        for queue, jobid in jobs:
            if queue.lower() in ("condor", "htcondor"):
                condor_jobs.append(jobid)
            else:
                slurm_jobs.append(jobid)
        if not cwait(condor_jobs, interval=interval, verbose=verbose):
            raise Exception("Not all HTCondor jobs finished successfully!")
        if not swait(slurm_jobs, interval=interval, verbose=verbose):
            raise Exception("Not all SLURM jobs finished successfully!")

    def info(self, msg, step=-1):
        """Print status message."""
        if self.verbose > 0:
            sys.stdout.write("{:%Y-%b-%d %H:%M:%S} INFO ".format(datetime.now()))
            sys.stdout.write(msg)
            if step >= 0:
                sys.stdout.write(" (step={})".format(step))
            sys.stdout.write("\n")

    def normtime(self, t):
        """Clamp time to be within atlas domain."""
        return round(min(max(t, self.means[0]), self.means[-1]), self.precision)

    def timeindex(self, t):
        """Get discrete time point index."""
        if t <= self.means[0]:
            return 0
        i = len(self.means) - 1
        if t >= self.means[-1]:
            return i
        for j in range(1, len(self.means)):
            if self.means[j] > t:
                i = j - 1
                break
        return i

    def timename(self, t1, t2=None):
        """Get time point string used in file names."""
        t1 = self.normtime(t1)
        if t2 is None:
            t2 = t1
        else:
            t2 = self.normtime(t2)
        if t1 == t2:
            return "t{0:05.{1}f}".format(self.normtime(t1), self.precision)
        return "{0:s}-{1:s}".format(self.timename(t1), self.timename(t2), self.precision)

    def age(self, imgid):
        """Get time associated with the specified image, i.e., mean of temporal Gaussian kernel."""
        return self.normtime(self.imgage[imgid])

    def ages(self):
        """Get set of all ages associated with the images from which atlas is constructed."""
        ages = set()
        for imgid in self.imgids:
            ages.add(self.age(imgid))
        return ages

    def stdev(self, t):
        """Get standard deviation of temporal Gaussian kernel centered at time t."""
        i = self.timeindex(t)
        if i == len(self.sigma) - 1:
            return self.sigma[i]
        else:
            alpha = t - self.means[i]
            return (1. - alpha) * self.sigma[i] + alpha * self.sigma[i + 1]

    def weight(self, t, mean, sigma=0, normalize=True, epsilon=None):
        """Evaluate temporal kernel weight for specified image."""
        if sigma <= 0:
            sigma = self.stdev(mean)
        w = math.exp(- .5 * math.pow((t - mean) / sigma, 2))
        if normalize:
            w /= sigma * math.sqrt(2. * math.pi)
        if epsilon is None:
            epsilon = self.epsilon
        return 0. if w < epsilon else w

    def weights(self, mean=None, sigma=0, zero=False):
        """Get weights of images within local support of given temporal kernel."""
        if mean is None:
            wmean = {}
            for mean in self.means:
                wimg = {}
                for imgid in self.imgids:
                    w = self.weight(self.age(imgid), mean=mean, sigma=sigma)
                    if zero or w > 0.:
                        wimg[imgid] = w
                wmean[mean] = wimg
            return wmean
        else:
            wimg = {}
            for imgid in self.imgids:
                w = self.weight(self.age(imgid), mean=mean, sigma=sigma)
                if zero or w > 0.:
                    wimg[imgid] = w
            return wimg

    def subdir(self, step=-1):
        """Get absolute path of directory of current iteration."""
        path = self.tmpdir
        if step < 0:
            step = self.step
        if step >= 0:
            path = os.path.join(path, "i{0:02d}".format(step))
        return path

    def image(self, imgid):
        """Get absolute path of individual image from which atlas is constructed."""
        return os.path.join(self.imgdir, self.imgpre + imgid + self.imgsuf)
    

    def affdof(self, imgid):
        """Get absolute path of affine transformation of global normalization."""
        if self.affdir:
            return os.path.join(self.affdir, self.affpre + imgid + self.affsuf)
        else:
            return "identity"

    def svffd(self, imgid, t=None, step=-1, path=None, force=False, create=True, batch=False):
        """Register atlas to specified image and return path of SV FFD file."""
        if step < 0:
            step = self.step
        age = self.age(imgid)
        if t is None or not self.age_specific_svffds:
            t = age
        if not path:
            path = os.path.join(self.subdir(step), "dof", "deformation", self.timename(t), "{0}.dof.gz".format(imgid))
        dof = path
        if t == age:
            if force or not os.path.exists(path):
                if step < 1:
                    dof = "identity"
                elif create:
                    args = ["-image", self.image(imgid)]
                    affdof = self.affdof(imgid)
                    if affdof != "identity":
                        args.extend(["-dof", affdof])
                    args.extend(["-image", self.avgimg(t, step=step - 1, create=not batch)])
                    opts = {
                        "parin": self.reg_config,
                        "mask": self.refimg,
                        "dofin": "identity",
                        "dofout": dof
                    }
                    makedirs(os.path.dirname(dof))
                    self._run("register", args=args, opts=opts)
        else:
            dof1 = self.svffd(imgid, t=age, step=step, force=force, create=create and not batch)
            dof2 = self.growth(age, t, step=step - 1, force=force, create=create and not batch)
            if dof1 == "identity" and dof2 == "identity":
                dof = "identity"
            elif dof1 == "identity":
                dof = dof2
            elif dof2 == "identity":
                dof = dof1
            elif create:
                makedirs(os.path.dirname(dof))
                self._run("compose-dofs", args=[dof1, dof2, dof], opts={"bch": self.num_bch_terms, "global": False})
        return dof

    def svffds(self, step=-1, force=False, create=True, queue=None):
        """Register all images to their age-specific template."""
        if step < 0:
            step = self.step
        self.info("Register images to template of corresponding age", step=step)
        script = ""
        tasks = 0
        dofs = {}
        for imgid in self.imgids:
            dof = self.svffd(imgid, step=step, force=force, create=create and not queue)
            if dof != "identity" and queue and (force or not os.path.exists(dof)):
                remove_or_makedirs(dof)
                script += 'elif taskid == {taskid}: atlas.svffd("{imgid}", batch=True)\n'.format(taskid=tasks, imgid=imgid)
                tasks += 1
            dofs[imgid] = dof
        if create and queue:
            return (dofs, self._submit("register".format(step), script=script, tasks=tasks, step=step, queue=queue))
        return dofs

    def doftable(self, t, step=-1, force=False, create=True, batch=False):
        """Write table with image to atlas deformations of images with non-zero weight."""
        dofs = []
        t = self.normtime(t)
        weights = self.weights(t, zero=False)
        all_dofs_are_identity = True
        for imgid in self.imgids:
            if imgid in weights:
                if not self.age_specific_svffds or t == self.age(imgid):
                    dof = self.svffd(imgid, t=t, step=step, force=force, create=create and not batch)
                else:
                    dof = self.svffd(imgid, t=t, step=step, force=force, create=create, batch=batch)
                dofs.append((dof, weights[imgid]))
                if dof != "identity":
                    all_dofs_are_identity = False
        if all_dofs_are_identity:
            table = "identity"
        elif len(dofs) > 0:
            table = os.path.join(self.subdir(step), "etc", "{}-dofs.tsv".format(self.timename(t)))
            if create and (force or not os.path.exists(table)):
                makedirs(os.path.dirname(table))
                with open(table, "wt") as f:
                    for dof, w in dofs:
                        f.write("{}\t{}\n".format(dof, w))
        else:
            raise ValueError("No image has non-zero weight for time {0}!".format(t))
        return (table, "", "")

    def avgdof(self, t, path=None, step=-1, force=False, create=True, batch=False):
        """Get mean cross-sectional SV FFD transformation at given time."""
        if not path:
            path = os.path.join(self.subdir(step), "dof", "average", "{0}.dof.gz".format(self.timename(t)))
        if create and (force or not os.path.exists(path)):
            table, prefix, suffix = self.doftable(t, step=step, force=force, batch=batch)
            if table == "identity":
                path = "identity"
            else:
                makedirs(os.path.dirname(path))
                self._run("average-dofs", args=[path], opts={"dofnames": table, "prefix": prefix, "suffix": suffix, "invert": None})
        return path

    def avgdofs(self, step=-1, force=False, create=True, queue=None):
        """Compute all average SV FFDs needed for (parallel) atlas construction."""
        if step < 0:
            step = self.step
        self.info("Compute residual average deformations at each age", step=step)
        script = ""
        tasks = 0
        dofs = {}
        for t in self.ages():
            dof = self.avgdof(t, step=step, force=force, create=create and not queue)
            if queue and (force or not os.path.exists(dof)):
                remove_or_makedirs(dof)
                script += 'elif taskid == {taskid}: atlas.avgdof({t}, batch=True)\n'.format(taskid=tasks, t=t)
                tasks += 1
            dofs[t] = dof
        if create and queue:
            return (dofs, self._submit("avgdofs", script=script, tasks=tasks, step=step, queue=queue))
        return dofs

    def growth(self, t1, t2, step=-1, force=False, create=True, batch=False):
        """Make composite SV FFD corresponding to longitudinal change from t1 to t2."""
        if step < 0:
            step = self.step
        t1 = self.normtime(t1)
        t2 = self.normtime(t2)
        if t1 == t2 or not self.residual_svffds:
            return "identity"
        dof = os.path.join(self.subdir(step), "dof", "growth", "{0}.dof.gz".format(self.timename(t1, t2)))
        if step < 1:
            return path if os.path.exists(dof) else "identity"
        if create and (force or not os.path.exists(dof)):
            dofs = [
                self.avgdof(t1, step=step, force=force, create=not batch),
                self.growth(t1, t2, step=step - 1, force=force, create=not batch),
                self.avgdof(t2, step=step, force=force, create=not batch)
            ]
            if dofs[1] == "identity":
                del dofs[1]
            makedirs(os.path.dirname(dof))
            self._run("compose-dofs", args=dofs + [dof], opts={"scale": -1., "bch": self.num_bch_terms, "global": False})
        return dof

    def compose(self, step=-1, ages=[], allpairs=False, force=False, create=True, queue=None):
        """Compose longitudinal deformations with residual average deformations."""
        if step < 0:
            step = self.step
        self.info("Update all pairs of longitudinal deformations", step=step)
        script = ""
        tasks = 0
        dofs = {}
        if not ages:
            ages = self.ages()
        for t1 in ages:
            dofs[t1] = {}
            for t2 in ages:
                if allpairs or (t1 != t2 and self.weight(t1, mean=t2) > 0.):
                    dof = self.growth(t1, t2, step=step, force=force, create=create and not queue)
                    if dof != "identity" and queue and (force or not os.path.exists(dof)):
                        remove_or_makedirs(dof)
                        script += 'elif taskid == {taskid}: atlas.growth({t1}, {t2}, batch=True)\n'.format(taskid=tasks, t1=t1, t2=t2)
                        tasks += 1
                    dofs[t1][t2] = dof
        if create and queue:
            return (dofs, self._submit("compose", script=script, tasks=tasks, step=step, queue=queue))
        return dofs

    def imgdof(self, imgid, t, step=-1, decomposed=False, force=False, create=True, batch=False):
        """Compute composite image to atlas transformation."""
        if step < 0:
            step = self.step
        if step > 0:
            dof = os.path.join(self.subdir(step), "dof", "composite", self.timename(t), "{0}.dof.gz".format(imgid))
            if decomposed or (create and (force or not os.path.exists(dof))):
                dofs = [
                    self.affdof(imgid),
                    self.svffd(imgid, t, step=step, force=force, create=create and not batch)
                ]
                if not self.age_specific_svffds and self.age_specific_imgdof:
                    growth = self.growth(self.age(imgid), t, step=step - 1, force=force, create=create and not batch)
                    if growth != "identity":
                        dofs.append(growth)
                if self.residual_svffds:
                    dofs.append(self.avgdof(t, step=step, force=force, create=create and not batch))
                if decomposed:
                    dof = dofs
                elif create:
                    makedirs(os.path.dirname(dof))
                    self._run("compose-dofs", args=dofs + [dof])
        else:
            dof = self.affdof(imgid)
            if decomposed:
                dof = [dof] + ['identity'] * 2
        return dof

    def imgtable(self, t, step=-1, decomposed=False, force=False, create=True, batch=False):
        """Compute average and composite transformations to deform each image to the specified time point."""
        table = os.path.join(self.subdir(step), "etc", "{}-imgs.tsv".format(self.timename(t)))
        if create and (force or not os.path.exists(table)):
            weights = self.weights(t)
            makedirs(os.path.dirname(table))
            with open(table, "wt") as f:
                f.write(self.topdir)
                f.write("\n")
                for imgid in self.imgids:
                    if imgid in weights:
                        f.write(os.path.relpath(self.image(imgid), self.topdir))
                        dofs = self.imgdof(imgid, t, step=step, decomposed=True, force=force, create=not batch)
                        if isinstance(dofs, str):
                            dofs = [dofs]
                        for dof in dofs:
                            if dof and dof != "identity":
                                f.write("\t" + os.path.relpath(dof, self.topdir))
                        f.write("\t{0}\n".format(weights[imgid]))
        return table

    def avgimg(self, t, path=None, step=-1, decomposed=False, force=False, create=True, batch=False):
        """Create average image for a given time point."""
        if not path:
            path = os.path.join(self.subdir(step), "img", "{0}.nii.gz".format(self.timename(t)))
        if create and (force or not os.path.exists(path)):
            table = self.imgtable(t, step=step, decomposed=decomposed, force=force, batch=batch)
            makedirs(os.path.dirname(path))
            self._run("average-images", args=[path], opts={
                "images": table,
                "padding": 0,
                "threshold": .5,
                "reference": self.refimg,
                "normalization": "zscore",
                "rescaling": (0, 100),
                "sharpen": True,
                "dtype": "uchar"
            })
        return path

    def avgimgs(self, step=-1, ages=[], outdir=None, decomposed=False, force=False, create=True, queue=None):
        """Create all average images required for (parallel) atlas construction."""
        if step < 0:
            step = self.step
        if ages:
            self.info("Create template for each specified time point", step=step)
        else:
            self.info("Create template for each represented time point", step=step)
            ages = self.ages()
        script = ""
        tasks = 0
        imgs = {}            
        for t in ages:
            if outdir:
                img = os.path.join(outdir, self.timename(t) + ".nii.gz")
            else:
                img = None
            img = self.avgimg(t, path=img, step=step, decomposed=decomposed, force=force, create=create and not queue)
            if queue and (force or not os.path.exists(img)):
                script += 'elif taskid == {taskid}: atlas.avgimg({t}, path="{path}", decomposed={decomposed}, batch=True)\n'.format(
                    taskid=tasks, t=t, path=img, decomposed=decomposed
                )
                tasks += 1
            imgs[t] = img
        if create and queue:
            return (imgs, self._submit("avgimgs", script=script, tasks=tasks, step=step, queue=queue))
        return imgs

    def construct(self, start=-1, niter=10, outdir=None, force=False, queue=None):
        """Perform atlas construction.

        Args:
            start (int): Last completed iteration. (default: step)
            niter (int, optional): Number of atlas construction iterations. (default: 10)
            outdir (str, optional): Directory for final templates. (default: subdir(step))
            force (bool, optional): Force re-creation of already existing files. (default: False)
            queue (str, tuple, optional): Name of queues of batch queuing system.
                When not specified, the atlas construction runs on the local machine
                using the number of threads specified during construction. When a single
                `str` is given, both short and long running jobs are submitted to the same
                queue. Otherwise, the first queue is used for shorter running jobs, while
                the second queue is used for longer running jobs.

        """
        if start < 0:
            start = self.step
        if start < 0:
            raise ValueError("Atlas to be constructed must have step index >= 0!")
        self.info("Performing {0} iterations starting with step {1}".format(niter, start))
        self.info("Average age-dependent deformations: {}".format(self.age_specific_svffds))
        self.info("Age-dependent image deformations:   {}".format(self.age_specific_imgdof))
        if not queue:
            queue = (None, None)
        elif not isinstance(queue, (list, tuple)):
            queue = (queue, queue)  # queues for (short, long) running jobs
        for step in range(start + 1, start + niter + 1):
            # Create template images using the current/initial transformations
            result = self.avgimgs(step=step - 1, decomposed=True, force=force, queue=queue[0])
            if queue[0]:
                self.wait(result[-1], interval=60, verbose=2)
            # Register all images to the current template images
            result = self.svffds(step=step, force=force, queue=queue[1])
            if queue[1]:
                self.wait(result[-1], interval=60, verbose=5)
            # Compute all required average deformations
            if self.residual_svffds:
                result = self.avgdofs(step=step, force=force, queue=queue[0])
                if queue[0]:
                    self.wait(result[-1], interval=60, verbose=2)
            # Compute all required longitudinal deformations
            if self.age_specific_svffds or self.age_specific_imgdof:
                result = self.compose(step=step, force=force, queue=queue[0])
                if queue[0]:
                    self.wait(result[-1], interval=30, verbose=1)
        # Write final template images to specified directory
        if not outdir:
            outdir = self.outdir
        self.info("Creating final mean shape and intensity templates")
        result = self.avgimgs(outdir=outdir, ages=self.means, step=start + niter,
                              decomposed=True, force=force, queue=queue[0])
        if queue[0]:
            self.wait(result[-1], interval=60, verbose=2)
        self.info("Finished atlas construction!")
        self.step = start + niter

    def template(self, i):
        """Get absolute path of i-th template image."""
        if i < 0 or i >= len(self.means):
            raise IndexError()
        img = os.path.join(self.outdir, self.timename(self.means[i]) + ".nii.gz")
        if self.step >= 0 and not os.path.exists(img):
            avg = self.avgimg(self.means[i], step=self.step, create=False)
            if os.path.exists(avg):
                img = avg
        return img

    def __len__(self):
        """Length of the atlas is the number of templates at discrete time points."""
        return len(self.means)

    def __getitem__(self, i):
        """Absolute file path of i-th atlas template."""
        return self.template(i)

    def deformation(self, i, t=None, force=False, create=True):
        """Get absolute path of longitudinal deformation from template i."""
        if i < 0 or i >= len(self.means):
            return "identity"
        if t is None:
            t = self.means[i + 1]
        return self.growth(self.means[i], t, step=self.step, force=force, create=create)

    def deform(self, i, t, path=None, force=False, create=True):
        """Deform i-th template using longitudinal deformations to time point t."""
        source = self.template(i)
        t = self.normtime(t)
        if t == self.means[i]:
            if path and os.path.realpath(os.path.abspath(path)) != os.path.realpath(source):
                if force and (create or not os.path.exists(path)):
                    shutil.copyfile(source, path)
                return path
            else:
                return source
        if not path:
            path = os.path.join(self.subdir(self.step), "img", "{}-{}.nii.gz".format(self.timename(self.means[i]), self.timename(t)))
        if create and (force or not os.path.exists(path)):
            dof = self.deformation(i, t)
            makedirs(os.path.dirname(path))
            self._run("transform-image", args=[source, path], opts={"dofin": dof, "target": self.refimg})
        return path

    def interpolate(self, t, path=None, interp="default", deform=False, sigma=0, force=False, create=True):
        """Interpolate atlas volume for specified time from finite set of templates.

        Unlike avgimg, this function does not evaluate the continuous spatio-temporal function
        to construct a template image for the given age. Instead, it uses the specified
        interpolation kernel and computes the corresponding weighted average of previously
        constructed template images.

        Args:
            interp (str): Temporal interpolation kernel (weights). (default: 'gaussian')
            sigma (float): Standard deviation used for Gaussian interpolation. (default: adaptive)
            deform (bool): Whether to deform each template using the longitudinal deformation.

        """
        interp = interp.lower()
        if interp in ("kernel", "default"):
            interp = "gaussian"
        t = self.normtime(t)
        i = self.timeindex(t)
        if t == self.means[i] and interp.startswith("linear"):
            return self.template(i)
        if not path:
            path = os.path.join(self.outdir, self.timename(t) + ".nii.gz")
        if create and (force or not os.path.exists(path)):
            args = [path]
            # Linear interpolation
            if interp == "linear":
                w = (self.means[i + 1] - t) / (self.means[i + 1] - self.means[i])
                args.extend(["-image", self.template(i), w])
                if deform:
                    args.extend(["-dof", self.deformation(i, t)])
                w = (t - self.means[i]) / (self.means[i + 1] - self.means[i])
                args.extend(["-image", self.template(i + 1), w])
                if deform:
                    args.extend(["-dof", self.deformation(i + 1, t)])
            # Gaussian interpolation
            elif interp == "gaussian":
                for i in range(len(self.means)):
                    w = self.weight(self.means[i], mean=t, sigma=sigma)
                    if w > 0.:
                        args.extend(["-image", self.template(i), w])
                        if deform:
                            dof = self.deformation(self.means[i], t)
                            if dof != "identity":
                                args.extend(["-dof", dof])
            if create:
                self._run("average-images", args=args, opts={"reference": self.refimg, "dtype": "uchar", "padding": 0})
        return path

    def view(self, i=None):
        """View template image at specified time point(s)."""
        if i is None:
            i = range(len(self.means))
        elif not isinstance(i, (list, tuple)):
            i = [i]
        imgs = [self.template(int(idx)) for idx in i]
        mirtk.run("view", opts={"target": imgs})

    def rmtemp(self):
        """Delete temporary files.

        This function removes all temporary files which can be recomputed
        without the need for performing the more costly registrations.
        Intermediate template images, auxiliary CSV files, and composite
        transformations are among these temporary files to be deleted.

        """
        raise NotImplementedError()


##############################################################################
# Batch execution using SLURM

# ----------------------------------------------------------------------------
def sbatch(name, args=[], opts={}, script=None, tasks=0, deps=[], logdir=None, log=None, queue='long', threads=8, verbose=1):
    if threads <= 0:
        raise ValueError("Must specify number of threads when executing as SLURM batch job!")
    if script:
        shexec = "#!/bin/bash\nexec {0} <(cat <<END_OF_SCRIPT\n".format(sys.executable)
        shexec += script.format(**opts)
        shexec += "\nEND_OF_SCRIPT)"
        if tasks > 0:
            shexec += " $SLURM_ARRAY_TASK_ID"
        shexec += "\n"
        script = shexec
    else:
        script = "#!/bin/bash\n"
        script += "mirtk {}".format(name)
        for arg in args:
            if ' ' in arg:
                arg = '"' + arg + '"'
            script += ' ' + str(arg)
        for opt in opts:
            arg = opts[opt]
            if isinstance(arg, (list, tuple)):
                arg = ' '.join(arg)
            if opt[0] != '-':
                opt = '-' + opt
            script += ' ' + opt + ' ' + str(arg)
        script += " -threads {0}".format(threads)
    argv = [
        'sbatch',
        '-J', name,
        '-n', '1',
        '-c', str(threads),
        '-p', queue,
        '--mem=8G'
    ]
    if tasks > 0:
        argv.append('--array=0-{}'.format(tasks - 1))
    if logdir or log:
        if not logdir:
            logdir = os.path.dirname(log)
        elif not log:
            log = os.path.join(logdir, name)
            if tasks > 0:
                log += "_%A.%a.log"
            log += "_%j.log"
        if not os.path.exists(logdir):
            os.makedirs(logdir)
        argv.extend(['-o', log, '-e', log])
    if deps:
        if isinstance(deps, int):
            deps = [deps]
        deps = [str(dep) for dep in deps if dep > 0]
        if deps:
            argv.append('--dependency=afterok:' + ',afterok:'.join(deps))
    proc = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    (out, err) = proc.communicate(input=script.encode('utf-8'))
    if proc.returncode != 0:
        raise Exception(err)
    match = re.match('Submitted batch job ([0-9]+)', out)
    if not match:
        raise Exception("Failed to determine job ID from sbatch output:\n" + out)
    jobid = int(match.group(1))
    if verbose > 0:
        if tasks > 0:
            print("  Submitted job {} (JobId={}, Tasks={})".format(name, jobid, tasks))
        else:
            print("  Submitted job {} (JobId={})".format(name, jobid))
    return jobid


# ----------------------------------------------------------------------------
def swait(jobs, max_time=0, interval=60, verbose=0):
    """Wait for SLURM jobs to finish."""
    if not isinstance(jobs, (list, tuple)):
        jobs = [jobs]
    jobs = [job for job in jobs if job > 0]
    num_wait = len(jobs)
    if num_wait == 0:
        return True
    num_fail = 0
    total_time = 0
    re_state = re.compile("JobState=([A-Z]+)")
    jobs = [str(job) for job in jobs]
    iterations = 0
    while num_wait > 0 and (max_time <= 0 or total_time < max_time):
        time.sleep(interval)
        total_time += interval
        done = []
        fail = []
        batch = []
        table = subprocess.check_output(["sacct", "--brief", "--parsable", "--jobs={}".format(','.join(jobs))])
        for line in table.splitlines():
            cols = line.split("|")
            if cols[0] in jobs:
                if cols[1] == "COMPLETED":
                    done.append(cols[0])
                elif cols[1] not in ("PENDING", "SUSPENDED", "RUNNING"):
                    fail.append(cols[0])
            elif cols[0].endswith(".batch"):
                batch.append(cols[0][0:-6])
        num_jobs = len(jobs)
        num_done = 0
        num_fail = 0
        for job in batch:
            try:
                info = subprocess.check_output(["scontrol", "show", "job", job])
                for line in info.splitlines():
                    match = re_state.search(line)
                    if match:
                        num_jobs += 1
                        if match.group(1) == "COMPLETED":
                            num_done += 1
                        elif match.group(1) not in ("PENDING", "SUSPENDED", "RUNNING"):
                            num_fail += 1
                num_jobs -= 1
                try:
                    done.remove(job)
                except ValueError:
                    pass
                try:
                    fail.remove(job)
                except ValueError:
                    pass
            except subprocess.CalledProcessError:
                pass  # scontrol forgets about no longer queued/running jobs
        num_done += len(done)
        num_fail += len(fail)
        num_wait = num_jobs - num_done - num_fail
        if verbose > 0 and (num_wait <= 0 or iterations % verbose == 0):
            sys.stdout.write("{:%Y-%b-%d %H:%M:%S}".format(datetime.now()))
            sys.stdout.write(" WAIT {} job(s) running/suspended/pending".format(num_wait))
            if num_fail > 0:
                sys.stdout.write(", {} failed".format(num_fail))
            sys.stdout.write("\n")
            sys.stdout.flush()
        ++iterations
    if num_wait > 0 and max_time > 0 and total_time >= max_time:
        raise Exception("Exceeded maximum time waiting for jobs to complete!")
    if total_time > 0:
        time.sleep(10)  # wait a bit for files to be available from all NFS clients
    return num_fail == 0


##############################################################################
# Batch execution using HTCondor

# ----------------------------------------------------------------------------
def cbatch(name, args=[], opts={}, script=None, tasks=0, deps=[], requirements=[], logdir=None, log=None, verbose=1):
    if deps:
        raise NotImplementedError("Cannot submit individual HTCondor jobs with dependencies, this requires use of DAGMan")
    if logdir or log:
        if not logdir:
            log = os.path.abspath(log)
            logdir = os.path.dirname(log)
        elif not log:
            logdir = os.path.abspath(logdir)
            if tasks > 0:
                log = os.path.join(logdir, name + "_$(Cluster).$(Process).log")
            else:
                log = os.path.join(logdir, name + "_$(Cluster).log")
        makedirs(logdir)
    jobdesc = "universe = vanilla\n"
    if requirements:
        jobdesc += "requirements = " + " && ".join(requirements) + "\n"
    if script:
        if not log:
            raise ValueError("Script submission of batch to HTCondor requires log path for script file!")
        script_path = os.path.join(logdir, name + ".py")
        with open(script_path, "wt") as f:
            f.write(script.format(**opts))
        jobdesc += "executable = {0}\n".format(sys.executable)
        jobdesc += "arguments = \"'{0}'".format(script_path)
        if tasks > 0:
            jobdesc += " $(Process)"
        jobdesc += "\"\n"
        if log:
            jobdesc += "output = {0}\n".format(log)
            jobdesc += "error = {0}\n".format(log)
        jobdesc += "queue"
        if tasks > 0:
            jobdesc += " {}".format(tasks)
        jobdesc += "\n"
    else:
        jobdesc += "executable = mirtk\n"
        jobdesc += "arguments = \"" + name
        for arg in args:
            if ' ' in arg:
                arg = "'" + arg + "'"
            jobdesc += ' ' + str(arg)
        for opt in opts:
            arg = opts[opt]
            if isinstance(arg, (list, tuple)):
                arg = ' '.join(arg)
            if opt[0] != '-':
                opt = '-' + opt
            jobdesc += ' ' + opt + ' ' + str(arg)
        jobdesc += "\"\n"
        if log:
            jobdesc += "output = {0}\n".format(log)
            jobdesc += "error = {0}\n".format(log)
        jobdesc += "queue\n"
    proc = subprocess.Popen(["condor_submit", "-"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    (out, err) = proc.communicate(input=jobdesc.encode('utf-8'))
    if proc.returncode != 0:
        raise Exception(err)
    match = re.search('[0-9]+ job(s) submitted to cluster ([0-9]+)\.', out)
    if not match:
        match = re.search('\*\* Proc ([0-9]+)(\.[0-9]+)?:', out)
        if not match:
            raise Exception("Failed to determine job ID from condor_submit output:\n" + out)
    jobid = int(match.group(1))
    if verbose > 0:
        if tasks > 0:
            print("  Submitted job {} (JobId={}, Tasks={})".format(name, jobid, tasks))
        else:
            print("  Submitted job {} (JobId={})".format(name, jobid))
    return jobid if tasks == 0 else (jobid, tasks)


# ----------------------------------------------------------------------------
def parse_condor_xml(s):
    """Parse XML ClassAds returned by condor_q or condor_history with -xml option."""
    # Note: condor_history may return multiple (only one non-empty) <classads> tags
    from xml.etree import ElementTree
    end = 0
    classads = None
    starttag = "<classads>"
    endtag = "</classads>"
    while True:
        start = s.find(starttag, end)
        if start == -1:
            break
        end = s.find(endtag, start)
        if end == -1:
            raise ValueError("Malformed <classads> XML, could not find matching </classads>!")
        end += len(endtag)
        elem = ElementTree.fromstring(s[start:end])
        if classads is None:
            classads = elem
        else:
            classads.extend(elem)
    return classads


# ----------------------------------------------------------------------------
def cwait(jobs, max_time=0, max_error=5, interval=60, verbose=0):
    if not isinstance(jobs, list):
        jobs = [jobs]
    clusters = []
    for job in jobs:
        if isinstance(job, tuple):
            cluid = job[0]
            tasks = job[1]
        else:
            cluid = job
            tasks = 1
        if cluid > 0:
            clusters.append((cluid, tasks))
    num_wait = len(clusters)
    num_error = 0
    total_time = 0
    iterations = 0
    while num_wait > 0 and (max_time <= 0 or total_time < max_time or total_time == 0) and num_error < max_error:
        time.sleep(interval)
        total_time += interval
        try:
            num_pending = 0
            num_running = 0
            num_suspended = 0
            num_held = 0
            num_done = 0
            for cluster, tasks in clusters:
                classads = parse_condor_xml(subprocess.check_output(["condor_q", "-xml", str(cluster)]))
                for process in range(tasks):
                    classad = classads.find(".c/a[@n='ClusterId'][i='{0}']/../a[@n='ProcId'][i='{1}']/..".format(cluster, process))
                    if classad is None:
                        num_done += 1
                    else:
                        status = int(classad.find("a[@n='JobStatus']/i").text)
                        # 1) Idle
                        # 2) Running
                        # 3) Removed
                        # 4) Completed (also when failed, check ExitCode afterwards using condor_history)
                        # 5) Held
                        # 6) Transferring Output
                        # 7) Suspended
                        if status == 1:
                            num_pending += 1
                        elif status == 2 or status == 6:
                            num_running += 1
                        elif status == 3:
                            num_done += 1
                        elif status == 4:
                            num_done += 1
                        elif status == 5:
                            num_held += 1
                        elif status == 7:
                            num_suspended += 1
                        else:
                            raise Exception("Unknown job status: {}".format(status))
            num_wait = num_running + num_pending + num_held + num_suspended
            if verbose > 0 and (num_wait <= 0 or iterations % verbose == 0):
                sys.stdout.write("{:%Y-%b-%d %H:%M:%S}".format(datetime.now()))
                sys.stdout.write(" WAIT {p} pending, {r} running, {s} suspended, {h} held, {d} completed\n".format(
                    p=num_pending, r=num_running, s=num_suspended, h=num_held, d=num_done
                ))
                sys.stdout.flush()
        except subprocess.CalledProcessError:
            sys.stdout.write("{:%Y-%b-%d %H:%M:%S}".format(datetime.now()))
            sys.stdout.write(" WAIT Failed to retrieve job status, will retry {0} more times!\n".format(max_error - num_error))
            sys.stdout.flush()
            num_error += 1
        iterations += 1
    if num_error >= max_error:
        raise Exception("Exceeded maximum number of retries to query status of jobs!")
    if num_wait > 0 and max_time > 0 and total_time >= max_time:
        raise Exception("Exceeded maximum time waiting for jobs to complete!")
    num_fail = 0
    if total_time > 0:
        time.sleep(10)
        num_jobs = 0
        unknown = {}
        for cluster, tasks in clusters:
            num_jobs += tasks
            unknown[cluster] = [0] * tasks
        num_error = 0
        num_unknown = num_jobs
        while num_unknown > 0 and num_error <= max_error:
            try:
                num_fail = 0
                for cluster, tasks in clusters:
                    classads = parse_condor_xml(subprocess.check_output(["condor_history", "-xml", str(cluster)]))
                    for process in range(tasks):
                        classad = classads.find(".c/a[@n='ClusterId'][i='{0}']/../a[@n='ProcId'][i='{1}']/..".format(cluster, process))
                        if classad is None:
                            unknown[cluster][process] += 1
                        else:
                            unknown[cluster][process] = 0
                            status = int(classad.find("a[@n='JobStatus']/i").text)
                            exit_code = int(classad.find("a[@n='ExitCode']/i").text)
                            if status != 4 or exit_code != 0:
                                num_fail += 1
                num_unknown = 0
                for cluster, tasks in clusters:
                    for process in range(tasks):
                        if unknown[cluster][process] > 0:
                            if unknown[cluster][process] > max_error:
                                raise Exception("Could not retrieve exit code of job {}.{} for the past {} attempts!".format(cluster, process, unknown[cluster][process]))
                            num_unknown += 1
                if verbose > 0:
                    sys.stdout.write("{:%Y-%b-%d %H:%M:%S}".format(datetime.now()))
                    sys.stdout.write(" DONE {0} succeeded, {1} failed".format(num_jobs - num_fail - num_unknown, num_fail))
                    if num_unknown > 0:
                        sys.stdout.write(", {} unknown, will retry".format(num_unknown))
                    sys.stdout.write("\n")
                    sys.stdout.flush()
            except subprocess.CalledProcessError:
                sys.stdout.write("{:%Y-%b-%d %H:%M:%S}".format(datetime.now()))
                sys.stdout.write(" WAIT Failed to retrieve exit codes, will retry {0} more times!\n".format(max_error - num_error))
                sys.stdout.flush()
                num_error += 1
        if num_error >= max_error:
            raise Exception("Exceeded maximum number of retries to query exit codes of jobs!")
    return num_fail == 0


##############################################################################
# Auxiliaries

# ----------------------------------------------------------------------------
def read_imgids(path):
    """Read subject/image IDs from text file."""
    imgids = []
    with open(path, "rt") as f:
        for line in f.readlines():
            line = line.strip()
            if len(line) > 0 and line[0] != '#':
                imgids.append(re.split("^[ \t,]+", line)[0])
    return imgids


# ----------------------------------------------------------------------------
def read_ages(path, delimiter=None):
    """Read subject/image age from CSV file."""
    if not delimiter:
        ext = os.path.splitext(path)[1].lower()
        if ext == ".csv":
            delimiter = ","
        elif ext == ".tsv":
            delimiter = "\t"
        elif ext == ".txt":
            delimiter = " "
        else:
            raise ValueError("Cannot determine delimiter of {} file! Use file extension .csv, .tsv, or .txt.".format(path))
    ages = {}
    with open(path, "rt") as csvfile:
        reader = csv.reader(csvfile, delimiter=delimiter, quotechar='"')
        for line in reader:
            ages[line[0]] = float(line[1])
    return ages


# ----------------------------------------------------------------------------
def basename_without_ext(path):
    name = os.path.basename(path)
    name, ext = os.path.splitext(name)
    if ext.lower() == '.gz':
        name = os.path.splitext(name)[0]
    return name


# ----------------------------------------------------------------------------
def makedirs(path):
    """Make directories, throws no error when it already exists."""
    if not path:
        raise ValueError("Path argument is empty or None!")
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

# ----------------------------------------------------------------------------
def remove_or_makedirs(path):
    """Remove file when it exists or make directories otherwise."""
    if os.path.exists(path):
        os.remove(path)
    else:
        makedirs(os.path.dirname(path))

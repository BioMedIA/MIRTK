##############################################################################
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright 2017 Imperial College London
# Copyright 2017 Andreas Schuh
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##############################################################################

import os
import sys
import math
import csv
import json
import re
import shutil
import subprocess
import time
from datetime import datetime

import mirtk
from mirtk.utils import makedirs
from mirtk.batch.slurm  import submit as sbatch, wait as swait
from mirtk.batch.condor import submit as cbatch, wait as cwait


##############################################################################
# Spatio-temporal atlas

class SpatioTemporalAtlas(object):
    """Spatio-temporal image and deformation atlas

    Configuration entries such as file and directory paths, final atlas time
    points (means) and corresponding temporal kernel width/standard deviation
    (sigma) are read from a JSON configuration file.

    The `construct` function performs the group-wise atlas construction.
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

    def __init__(self, config, root=None, step=-1, verbose=1, threads=-1, exit_on_error=False):
        """Load spatio-temporal atlas configuration."""
        self.step = step
        if isinstance(config, str):
            config = os.path.abspath(config)
            with open(config, "rt") as f:
                self.config = json.load(f)
            self.root = os.path.dirname(config)
        else:
            self.config = config
            if root:
                self.root = os.path.abspath(root)
            else:
                self.root = os.getcwd()
        # Paths of image file and output of prior global normalization step
        paths = self.config.get("paths", {})
        images = self.config.get("images", {})
        self.topdir = os.path.normpath(os.path.join(self.root, paths.get("topdir", ".")))
        self.agecsv = self._path(paths, "agecsv", os.path.join(self.topdir, "config", "ages.csv"))
        self.imgcsv = self._path(paths, "imgcsv", os.path.join(self.topdir, "config", "subjects.csv"))
        self.imgage = read_ages(self.agecsv)
        self.imgids = read_imgids(self.imgcsv)
        self.tmpdir = self._path(paths, "tmpdir", "cache")
        self.outdir = self._path(paths, "outdir", os.path.join(self.topdir, "templates"))
        self.refimg = self._path(images, "ref", os.path.join(self.topdir, "global", "average.nii.gz"))
        self.channel = images.get("default", "t2w")
        if len(self.config["images"]) == 1:
            self.channel = self.config["images"].keys()[0]
        # Parameters of temporal kernel regressions
        regression = self.config.get("regression", {})
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
        regcfg = self.config.get("registration", {})
        growth = regcfg.get("growth", {})
        self.num_bch_terms = growth.get("bchterms", 3)  # composition should be symmetric, e.g., 2, 3, or 5 terms
        self.age_specific_imgdof = growth.get("enabled", True)
        if self.age_specific_imgdof:
            self.age_specific_regdof = not growth.get("exclavg", False)
        else:
            self.age_specific_regdof = False
        # Other workflow execution settings
        envcfg = self.config.get("environment", {})
        self.verbose = verbose
        self.queue = {
            "short": envcfg.get("queue", {}).get("short", "local").lower(),
            "long": envcfg.get("queue", {}).get("long", "local").lower()
        }
        self.threads = envcfg.get("threads", 8) if threads < 0 else threads
        self.mintasks = envcfg.get("mintasks", 1)
        self.maxtasks = envcfg.get("maxtasks", 1000)
        self.exit_on_error = exit_on_error
        # Discard images not required for any final time point
        imgids = set()
        for t in self.means:
            imgids |= set(self.weights(mean=t).keys())
        self.imgids = list(imgids)
        self.imgids.sort()

    def _run(self, command, args=[], opts={}, step=-1, workdir=None, queue=None, name=None, submit_kwargs={}, wait_kwargs={}):
        """Execute single MIRTK command."""
        if step < 0:
            step = self.step
        if not name:
            name = command
        if queue and queue.lower() == "local":
            queue = None
        if command in ("edit-dofs", "em-hard-segmentation", "average-measure"):
            # These commands do not support the -threads option (yet)
            threads = 0
        else:
            threads = self.threads
            if "verbose" not in opts:
                if isinstance(opts, list):
                    opts.append(("verbose", self.verbose - 2))
                else:
                    opts["verbose"] = self.verbose - 2
        if "verbose" in opts and command in ("average-measure"):
            # These commands do not support the -verbose option (yet)
            del opts["verbose"]
        if queue:
            job = self._submit(name, command=command, args=args, opts=opts, step=step, workdir=workdir, script=None, tasks=-1, group=1, queue=queue, **submit_kwargs)
            self.wait(job, **wait_kwargs)
        else:
            prevdir = os.getcwd()
            if workdir:
                os.chdir(workdir)
            try:
                if self.verbose > 1:
                    sys.stdout.write("\n\n")
                mirtk.run(command, args=args, opts=opts, showcmd=(self.verbose > 1), threads=threads, onerror='exit' if self.exit_on_error else 'throw')
            finally:
                os.chdir(prevdir)

    def _submit(self, name, command=None, args=[], opts={}, script=None, tasks=-1, group=1, step=-1, queue=None, memory=8 * 1024, workdir=None):
        """Submit batch script."""
        if step < 0:
            step = self.step
        if not queue:
            queue = self.queue
        if tasks == 0:
            return (queue, 0)
        groups = tasks
        threads = self.threads if self.threads > 0 else 1
        if script:
            source = 'import sys\nimport socket\n'
            source += 'sys.stdout.write("Host: " + socket.gethostname() + "\\n\\n")\n'
            source += 'sys.path.insert(0, "{0}")\n'.format(os.path.dirname(os.path.dirname(mirtk.__file__)))
            source += 'sys.path.insert(0, "{0}")\n'.format(os.path.dirname(__file__))
            source += 'from {0} import SpatioTemporalAtlas\n'.format(os.path.splitext(os.path.basename(__file__))[0])
            source += 'atlas = SpatioTemporalAtlas(root="{root}", config={config}, step={step}, threads={threads}, verbose=3, exit_on_error=True)\n'
            if tasks > 0:
                tasks_per_group = max(group, self.mintasks, (tasks + self.maxtasks - 1) // self.maxtasks)
                if tasks_per_group > 1:
                    groups = range(0, tasks, tasks_per_group)
                    source += 'groupid = int(sys.argv[1])\n'
                    source += 'if groupid < 0 or groupid >= {0}:\n'.format(len(groups))
                    source += '    sys.stderr.write("Invalid group ID\\n")\n'
                    source += 'tasks_per_group = {0}\n'.format(tasks_per_group)
                    source += 'for taskid in range(groupid * tasks_per_group, (groupid + 1) * tasks_per_group):\n'
                    source += '    if taskid >= {0}: break\n'.format(tasks)
                    source += '    ' + '\n    '.join(script.splitlines()) + '\n'
                    source += '    else: sys.stderr.write("Invalid task ID\\n")\n'
                    groups = len(groups)
                else:
                    source += 'taskid = int(sys.argv[1])\n'
                    source += 'if taskid < 0: sys.stderr.write("Invalid task ID\\n")\n'
                    source += script
                    source += 'else: sys.stderr.write("Invalid task ID\\n")\n'
            else:
                source += script
            opts.update({
                "root": self.root,
                "config": repr(self.config),
                "step": step,
                "threads": threads
            })
            script = source
        elif not command:
            command = name
        jobname = "i{0:02d}_{1}".format(step, name) if step >= 0 else name
        if queue.lower() in ("condor", "htcondor"):
            if tasks > 0:
                log = os.path.join(self.subdir(step), "log", name + "_$(Cluster).$(Process).log")
            else:
                log = os.path.join(self.subdir(step), "log", name + "_$(Cluster).log")
            condor_config = self.config.get("environment", {}).get("condor", {})
            requirements = condor_config.get("requirements", [])
            environment = condor_config.get("environment", {})
            jobid = cbatch(name=jobname, command=command, args=args, opts=opts, script=script, tasks=groups,
                           log=log, threads=threads, memory=memory, requirements=requirements, environment=environment,
                           workdir=workdir, verbose=0)
        else:
            if tasks > 0:
                log = os.path.join(self.subdir(step), "log", name + "_%A.%a.log")
            else:
                log = os.path.join(self.subdir(step), "log", name + "_%j.log")
            jobid = sbatch(name=jobname, command=command, args=args, opts=opts, script=script, tasks=groups,
                           log=log, threads=threads, memory=memory, queue=queue,
                           workdir=workdir, verbose=0)
        if tasks > 0:
            self.info("Submitted batch '{}' (id={}, #jobs={}, #tasks={})".format(
                name, jobid[0] if isinstance(jobid, tuple) else jobid, groups, tasks)
            )
        else:
            self.info("Submitted job '{}' (id={})".format(name, jobid[0] if isinstance(jobid, tuple) else jobid))
        return (queue, jobid)

    def wait(self, jobs, interval=60, verbose=5):
        """Wait for batch jobs to complete."""
        if not isinstance(jobs, list):
            jobs = [jobs]
        condor_jobs = []
        slurm_jobs = []
        for queue, jobid in jobs:
            if queue and queue.lower() != "local":
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
        if isclose(t1, t2):
            w = 2
            if self.precision > 0:
                w += self.precision + 1
            return "t{0:0{1}.{2}f}".format(self.normtime(t1), w, self.precision)
        return "{0:s}-{1:s}".format(self.timename(t1), self.timename(t2), self.precision)

    def age(self, imgid):
        """Get time associated with the specified image, i.e., mean of temporal Gaussian kernel."""
        return self.normtime(self.imgage[imgid])

    def ages(self):
        """Get set of all ages associated with the images from which atlas is constructed."""
        ages = set()
        for imgid in self.imgids:
            ages.add(self.age(imgid))
        return list(ages)

    def stdev(self, t):
        """Get standard deviation of temporal Gaussian kernel centered at time t."""
        i = self.timeindex(t)
        if i == len(self.sigma) - 1:
            return self.sigma[i]
        else:
            alpha = (t - self.means[i]) / (self.means[i + 1] - self.means[i])
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

    def splitchannel(self, channel=None):
        """Split channel specification into name and label specification."""
        if not channel:
            channel = self.channel
        if "=" in channel:
            channel, label = channel.split("=")
            try:
                l = int(label)
                label = l
            except ValueError:
                label = label.split(",")
        else:
            label = 0
        return (channel, label)

    def image(self, imgid, channel=None, label=0, force=False, create=True):
        """Get absolute path of requested input image."""
        if not channel:
            channel = self.channel
        cfg = self.config["images"][channel]
        prefix = cfg["prefix"].replace("/", os.path.sep)
        suffix = cfg.get("suffix", ".nii.gz")
        img = os.path.normpath(os.path.join(self.topdir, prefix + imgid + suffix))
        if label:
            if isinstance(label, (tuple, list)):
                lblstr = ','.join([str(l).strip() for l in label])
            else:
                lblstr = str(label).strip()
            lblstr = lblstr.replace("..", "-")
            msk = os.path.join(self.topdir, "masks", "{}_{}".format(channel, lblstr), imgid + ".nii.gz")
            if create and (force or not os.path.exists(msk)):
                makedirs(os.path.dirname(msk))
                if not isinstance(label, list):
                    label = [label]
                self._run("calculate-element-wise", args=[img, "-label"] + label + ["-set", 1, "-pad", 0, "-out", msk, "binary"])
            img = msk
        return img

    def affdof(self, imgid):
        """Get absolute path of affine transformation of global normalization."""
        cfg = self.config.get("registration", {}).get("affine", None)
        if cfg is None:
            return "identity"
        prefix = cfg["prefix"].replace("/", os.path.sep)
        suffix = cfg.get("suffix", ".dof")
        return os.path.normpath(os.path.join(self.topdir, prefix + imgid + suffix))

    def regcfg(self, step=-1):
        """Get registration configuration entries."""
        configs = self.config.get("registration", {}).get("config", {})
        if not isinstance(configs, list):
            configs = [configs]
        cfg = {}
        for i in range(min(step, len(configs)) if step > 0 else len(configs)):
            cfg.update(configs[i])
        return cfg

    def parin(self, step=-1, force=False, create=True):
        """Write registration configuration file and return path."""
        if step < 0:
            step = self.step
        parin = os.path.join(self.subdir(step), "config", "register.cfg")
        if create and (force or not os.path.exists(parin)):
            cfg = self.regcfg(step)
            images = self.config["images"]
            channels = cfg.get("channels", ["t2w"])
            if not isinstance(channels, list):
                channels = [channels]
            if len(channels) == 0:
                raise ValueError("Empty list of images/channels specified for registration!")
            params = ["[default]"]
            # transformation model
            model = cfg.get("model", "SVFFD")
            mffd, model = model.split(":") if ":" in model else "None", model
            params.append("Transformation model = {}".format(model))
            params.append("Multi-level transformation = {}".format(mffd))
            # energy function
            energy = cfg.get("energy", "sym" if "svffd" in model.lower() else "asym")
            if energy.lower() in ("asym", "asymmetric"):
                sim_term_type = 0
            elif energy.lower() in ("ic", "inverseconsistent", "inverse-consistent"):
                sim_term_type = 1
            elif energy.lower() in ("sym", "symmetric"):
                sim_term_type = 2
            else:
                sim_term_type = -1
            if sim_term_type < 0:
                formula = energy
            else:
                formula = ""
                measures = cfg.get("measures", {})
                for c in range(len(channels)):
                    target = 2 * c + 1
                    source = 2 * c + 2
                    if isinstance(measures, list):
                        measure = measures[c] if c < len(measures) else measures[-1]
                    elif isinstance(measures, dict):
                        channel = self.splitchannel(channels[c])[0]
                        measure = measures.get(channel, "NMI")
                    else:
                        measure = measures
                    if sim_term_type == 2:
                        term = "{sim}[{channel} sim](I({tgt}) o T^-0.5, I({src}) o T^0.5)"
                    elif sim_term_type == 1:
                        term = "{sim}[{channel} fwd-sim](I({tgt}) o T^-1, I({src})) + {sim}[{channel} bwd-sim](I({tgt}), I({src}) o T)"
                    else:
                        term = "{sim}[{channel} sim](I({tgt}), I({src}) o T)"
                    if c > 0:
                        formula += " + "
                    formula += term.format(sim=measure, channel=channels[c].capitalize().replace("=", " "), tgt=target, src=source)
            if "bending" in cfg:
                formula += " + 0 BE[Bending energy](T)"
            if "jacobian" in cfg:
                formula += " + 0 JAC[Jacobian penalty](T)"
            params.append("Energy function = " + formula)
            params.append("No. of bins = {}".format(cfg.get("bins", 64)))
            params.append("Local window size [box] = {}".format(cfg.get("window", "5 vox")))
            params.append("No. of last function values = 10")
            if "svffd" in model.lower():
                params.append("Integration method = {}".format(cfg.get("ffdim", "FastSS")))
                params.append("No. of integration steps = {}".format(cfg.get("intsteps", 64)))
                params.append("No. of BCH terms = {}".format(cfg.get("bchterms", 4)))
                params.append("Use Lie derivative = {}".format(cfg.get("liederiv", "No")))
            # resolution pyramid
            spacings = cfg.get("spacing", [])
            if not isinstance(spacings, list):
                spacings = [spacings]
            resolutions = cfg.get("resolution", [])
            if not isinstance(resolutions, list):
                resolutions = [resolutions]
            blurring = cfg.get("blurring", {})
            if isinstance(blurring, list):
                blurring = {self.splitchannel(channels[0])[0]: blurring}
            be_weights = cfg.get("bending", [0])
            if not isinstance(be_weights, list):
                be_weights = [be_weights]
            lj_weights = cfg.get("jacobian", [0])
            if not isinstance(lj_weights, list):
                lj_weights = [lj_weights]
            levels = cfg.get("levels", 0)
            if levels <= 0:
                levels = len(resolutions) if isinstance(resolutions, list) else 4
            resolution = 0.
            spacing = 0.
            params.append("No. of resolution levels = {0}".format(levels))
            params.append("Image interpolation mode = {0}".format(cfg.get("interpolation", "Linear with padding")))
            params.append("Downsample images with padding = {0}".format("Yes" if cfg.get("padding", True) else "No"))
            params.append("Image similarity foreground = {0}".format(cfg.get("foreground", "Overlap")))
            params.append("Strict step length range = No")
            params.append("Maximum streak of rejected steps = 2")
            for level in range(1, levels + 1):
                params.append("")
                params.append("[level {}]".format(level))
                resolution = float(resolutions[level - 1]) if level <= len(resolutions) else 2. * resolution
                if resolution > 0.:
                    params.append("Resolution [mm] = {}".format(resolution))
                for c in range(len(channels)):
                    target = 2 * c + 1
                    source = 2 * c + 2
                    image = images[self.splitchannel(channels[c])[0]]
                    bkgrnd = float(image.get("bkgrnd", -1))
                    params.append("Background value of image {0} = {1}".format(target, bkgrnd))
                    bkgrnd = -1. if "labels" in image else 0.
                    params.append("Background value of image {0} = {1}".format(source, bkgrnd))
                for c in range(len(channels)):
                    target = 2 * c + 1
                    source = 2 * c + 2
                    channel = self.splitchannel(channels[c])[0]
                    sigmas = blurring.get(channel, [])
                    if not isinstance(sigmas, list):
                        sigmas = [sigmas]
                    if len(sigmas) == 0 and "labels" in images[channel]:
                        sigmas = [2]
                    if len(sigmas) > 0:
                        sigma = float(sigmas[level - 1] if level <= len(sigmas) else sigmas[-1])
                        params.append("Blurring of image {0} [vox] = {1}".format(target, sigma))
                spacing = float(spacings[level - 1]) if level <= len(spacings) else 2. * spacing
                if spacing > 0.:
                    params.append("Control point spacing = {0}".format(spacing))
                be_weight = float(be_weights[level - 1] if level <= len(be_weights) else be_weights[-1])
                params.append("Bending energy weight = {0}".format(be_weight))
                lj_weight = float(lj_weights[level - 1] if level <= len(lj_weights) else lj_weights[-1])
                params.append("Jacobian penalty weight = {0}".format(lj_weight))
            # write -parin file for "register" command
            makedirs(os.path.dirname(parin))
            with open(parin, "wt") as f:
                f.write("\n".join(params) + "\n")
        return parin

    def regdof(self, imgid, t=None, step=-1, path=None, force=False, create=True, batch=False):
        """Register atlas to specified image and return path of transformation file."""
        if step < 0:
            step = self.step
        age = self.age(imgid)
        if t is None or not self.age_specific_regdof:
            t = age
        if not path:
            path = os.path.join(self.subdir(step), "dof", "deformation", self.timename(t), "{0}.dof.gz".format(imgid))
        dof = path
        if isclose(t, age):
            if force or not os.path.exists(path):
                if step < 1:
                    dof = "identity"
                elif create:
                    cfg = self.regcfg(step)
                    args = []
                    affdof = self.affdof(imgid)
                    if affdof == "identity":
                        affdof = None
                    channels = cfg.get("channels", self.channel)
                    if not isinstance(channels, list):
                        channels = [channels]
                    if len(channels) == 0:
                        raise ValueError("Empty list of images/channels specified for registration!")
                    for channel in channels:
                        channel, label = self.splitchannel(channel)
                        args.append("-image")
                        args.append(self.image(imgid, channel=channel, label=label))
                        if affdof:
                            args.append("-dof")
                            args.append(affdof)
                        args.append("-image")
                        args.append(self.avgimg(t, channel=channel, label=label, step=step - 1, create=not batch))
                    makedirs(os.path.dirname(dof))
                    self._run("register", args=args, opts={
                        "parin": self.parin(step=step, force=force, create=not batch),
                        "mask": self.refimg,
                        "dofin": "identity",
                        "dofout": dof
                    })
        else:
            dof1 = self.regdof(imgid, t=age, step=step, force=force, create=create and not batch)
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

    def regdofs(self, step=-1, force=False, create=True, queue=None, batchname="register"):
        """Register all images to their age-specific template."""
        if step < 0:
            step = self.step
        if not queue:
            queue = self.queue["long"]
        if queue == "local":
            queue = None
        self.info("Register images to template of corresponding age", step=step)
        script = ""
        tasks = 0
        dofs = {}
        for imgid in self.imgids:
            dof = self.regdof(imgid, step=step, force=force, create=create and not queue)
            if dof != "identity" and queue and (force or not os.path.exists(dof)):
                remove_or_makedirs(dof)
                script += 'elif taskid == {taskid}: atlas.regdof("{imgid}", batch=True)\n'.format(taskid=tasks, imgid=imgid)
                tasks += 1
            dofs[imgid] = dof
        if create and queue:
            self.parin(step=step, force=force)  # create -parin file if missing
            job = self._submit(batchname, script=script, tasks=tasks, step=step, queue=queue)
        else:
            job = (None, 0)
        return (job, dofs)

    def doftable(self, t, step=-1, force=False, create=True, batch=False):
        """Write table with image to atlas deformations of images with non-zero weight."""
        dofs = []
        t = self.normtime(t)
        weights = self.weights(t)
        all_dofs_are_identity = True
        for imgid in self.imgids:
            if imgid in weights:
                if not self.age_specific_regdof or isclose(t, self.age(imgid)):
                    dof = self.regdof(imgid, t=t, step=step, force=force, create=create and not batch)
                else:
                    dof = self.regdof(imgid, t=t, step=step, force=force, create=create, batch=batch)
                dofs.append((dof, weights[imgid]))
                if dof != "identity":
                    all_dofs_are_identity = False
        if all_dofs_are_identity:
            dofdir = None
            dofnames = "identity"
        elif len(dofs) > 0:
            dofdir = self.topdir
            dofnames = os.path.join(self.subdir(step), "config", "{}-dofs.tsv".format(self.timename(t)))
            if create and (force or not os.path.exists(dofnames)):
                makedirs(os.path.dirname(dofnames))
                with open(dofnames, "wt") as table:
                    for dof, w in dofs:
                        if dofdir:
                            dof = os.path.relpath(dof, dofdir)
                        table.write("{}\t{}\n".format(dof, w))
        else:
            raise ValueError("No image has non-zero weight for time {0}!".format(t))
        return (dofdir, dofnames)

    def avgdof(self, t, path=None, step=-1, force=False, create=True, batch=False):
        """Get mean cross-sectional SV FFD transformation at given time."""
        t = self.normtime(t)
        if not path:
            path = os.path.join(self.subdir(step), "dof", "average", "{0}.dof.gz".format(self.timename(t)))
        if create and (force or not os.path.exists(path)):
            dofdir, dofnames = self.doftable(t, step=step, force=force, batch=batch)
            if dofnames == "identity":
                path = "identity"
            else:
                makedirs(os.path.dirname(path))
                self._run("average-dofs", args=[path], opts={
                    "dofdir": dofdir,
                    "dofnames": dofnames,
                    "type": "SVFFD",
                    "target": "common",
                    "invert": True,
                    "global": False
                })
        return path

    def avgdofs(self, step=-1, force=False, create=True, queue=None, batchname="avgdofs"):
        """Compute all average SV FFDs needed for (parallel) atlas construction."""
        if step < 0:
            step = self.step
        if not queue:
            queue = self.queue["short"]
        if queue == "local":
            queue = None
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
            job = self._submit(batchname, script=script, tasks=tasks, step=step, queue=queue)
        else:
            job = (None, 0)
        return (job, dofs)

    def growth(self, t1, t2, step=-1, force=False, create=True, batch=False):
        """Make composite SV FFD corresponding to longitudinal change from t1 to t2."""
        if step < 0:
            step = self.step
        t1 = self.normtime(t1)
        t2 = self.normtime(t2)
        if isclose(t1, t2):
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

    def compose(self, step=-1, ages=[], allpairs=False, force=False, create=True, queue=None, batchname="compose"):
        """Compose longitudinal deformations with residual average deformations."""
        if step < 0:
            step = self.step
        if not queue:
            queue = self.queue["short"]
        if queue == "local":
            queue = None
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
            job = self._submit(batchname, script=script, tasks=tasks, step=step, queue=queue)
        else:
            job = (None, 0)
        return (job, dofs)

    def imgdof(self, imgid, t, step=-1, decomposed=False, force=False, create=True, batch=False):
        """Compute composite image to atlas transformation."""
        if step < 0:
            step = self.step
        if step > 0:
            dof = os.path.join(self.subdir(step), "dof", "composite", self.timename(t), "{0}.dof.gz".format(imgid))
            if decomposed or (create and (force or not os.path.exists(dof))):
                dofs = [
                    self.affdof(imgid),
                    self.regdof(imgid, step=step, force=force, create=create and not batch)
                ]
                if self.age_specific_imgdof:
                    growth = self.growth(self.age(imgid), t, step=step - 1, force=force, create=create and not batch)
                    if growth != "identity":
                        dofs.append(growth)
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

    def imgdofs(self, ages=[], step=-1, force=False, create=True, queue=None, batchname="imgdofs"):
        """Compute all composite image to atlas transformations."""
        if step < 0:
            step = self.step
        if not queue:
            queue = self.queue["short"]
        if queue == "local":
            queue = None
        self.info("Update composite image to atlas transformations", step=step)
        if ages:
            if not isinstance(ages, (tuple, list)):
                ages = [ages]
            ages = [self.normtime(t) for t in ages]
        else:
            ages = self.ages()
        script = ""
        tasks = 0
        dofs = {}
        for t in ages:
            dofs[t] = {}
            weights = self.weights(t)
            for imgid in self.imgids:
                if imgid in weights:
                    dof = self.imgdof(imgid=imgid, t=t, step=step, force=force, create=create and not queue)
                    if dof != "identity" and queue and (force or not os.path.exists(dof)):
                        remove_or_makedirs(dof)
                        script += 'elif taskid == {taskid}: atlas.imgdof(imgid="{imgid}", t={t}, batch=True)\n'.format(
                            taskid=tasks, imgid=imgid, t=t
                        )
                        tasks += 1
                    dofs[t][imgid] = dof
        if create and queue:
            job = self._submit(batchname, script=script, tasks=tasks, step=step, queue=queue)
        else:
            job = (None, 0)
        return (job, dofs)

    def defimg(self, imgid, t, channel=None, path=None, step=-1, decomposed=True, force=False, create=True, batch=False):
        """Transform sample image to atlas space at given time point."""
        if not channel:
            channel = self.channel
        if not path:
            path = os.path.join(self.subdir(step), channel, self.timename(t), imgid + ".nii.gz")
        if create and (force or not os.path.exists(path)):
            cfg = self.config["images"][channel]
            img = self.image(imgid, channel=channel)
            dof = self.imgdof(imgid=imgid, t=t, step=step, decomposed=decomposed, force=force, create=not batch)
            opts = {
                "interp": cfg.get("interp", "linear"),
                "target": self.refimg,
                "dofin": dof,
                "invert": True
            }
            if "bkgrnd" in cfg:
                opts["source-padding"] = cfg["bkgrnd"]
            if "labels" in cfg:
                opts["labels"] = cfg["labels"].split(",")
            if "datatype" in cfg:
                opts["datatype"] = cfg["datatype"]
            makedirs(os.path.dirname(path))
            self._run("transform-image", args=[img, path], opts=opts)
        return path

    def defimgs(self, ages=[], channels=[], step=-1, decomposed=True, force=False, create=True, queue=None, batchname="defimgs"):
        """Transform all images to discrete set of atlas time points."""
        if step < 0:
            step = self.step
        if not queue:
            queue = self.queue["short"]
        if queue == "local":
            queue = None
        single_channel = channels and isinstance(channels, basestring)
        if not channels:
            channels = self.regcfg(step).get("channels", self.channel)
        if not isinstance(channels, list):
            channels = [channels]
        if ages:
            self.info("Deform images to discrete time points", step=step)
            if not isinstance(ages, (tuple, list)):
                ages = [ages]
            ages = [self.normtime(t) for t in ages]
        else:
            self.info("Deform images to observed time points", step=step)
            ages = self.ages()
        script = ""
        tasks = 0
        imgs = {}
        for channel in channels:
            imgs[channel] = {}
            for t in ages:
                weights = self.weights(t)
                if len(weights) == 0:
                    raise Exception("No image has non-zero weight for t={}".format(t))
                imgs[channel][t] = []
                for imgid in self.imgids:
                    if imgid in weights:
                        img = self.defimg(imgid=imgid, t=t, channel=channel, step=step, decomposed=decomposed, force=force, create=create and not queue)
                        if queue and (force or not os.path.exists(img)):
                            remove_or_makedirs(img)
                            script += 'elif taskid == {taskid}: atlas.defimg(imgid="{imgid}", t={t}, channel="{channel}", path="{path}", decomposed={decomposed}, batch=True)\n'.format(
                                taskid=tasks, t=t, imgid=imgid, channel=channel, path=img, decomposed=decomposed
                            )
                            tasks += 1
                        imgs[channel][t].append(img)
        if create and queue:
            job = self._submit(batchname, script=script, tasks=tasks, step=step, queue=queue)
        else:
            job = (None, 0)
        if single_channel:
            imgs = imgs[channel]
        return (job, imgs)

    def imgtable(self, t, channel=None, step=-1, decomposed=True, force=False, create=True, batch=False):
        """Write image table with weights for average-images, and deform images to given time point."""
        t = self.normtime(t)
        if not channel:
            channel = self.channel
        image = self.config["images"][channel]
        table = os.path.join(self.subdir(step), "config", "{t}-{channel}.tsv".format(t=self.timename(t), channel=channel))
        if create and (force or not os.path.exists(table)):
            weights = self.weights(t)
            if len(weights) == 0:
                raise Exception("No image has non-zero weight for t={}".format(t))
            makedirs(os.path.dirname(table))
            with open(table, "wt") as f:
                f.write(self.topdir)
                f.write("\n")
                for imgid in self.imgids:
                    if imgid in weights:
                        img = self.defimg(t=t, imgid=imgid, channel=channel, step=step, decomposed=decomposed, force=force, create=not batch)
                        f.write(os.path.relpath(img, self.topdir))
                        f.write("\t{0}\n".format(weights[imgid]))
        return table

    def avgimg(self, t, channel=None, label=0, path=None, sharpen=True, outdir=None, step=-1, decomposed=True, force=False, create=True, batch=False):
        """Create average image for a given time point."""
        if not channel:
            channel = self.channel
        cfg = self.config["images"][channel]
        if isinstance(label, basestring):
            label = label.split(",")
        if isinstance(label, list) and len(label) == 1:
            label = label[0]
        if not path:
            if label:
                if isinstance(label, (tuple, list)):
                    lblstr = ','.join([str(l).strip() for l in label])
                elif isinstance(label, int):
                    max_label = max(parselabels(cfg["labels"])) if "labels" in cfg else 9
                    lblstr = "{0:0{1}d}".format(label, len(str(max_label)))
                else:
                    lblstr = str(label).strip()
                lblstr = lblstr.replace("..", "-")
            else:
                lblstr = None
            if outdir:
                outdir = os.path.abspath(outdir)
                if "labels" in cfg:
                    if lblstr:
                        path = os.path.join(outdir, "pbmaps", "_".join([channel, lblstr]))
                    else:
                        path = os.path.join(outdir, "labels", channel)
                else:
                    path = os.path.join(outdir, "templates", channel)
            else:
                path = os.path.join(self.subdir(step), channel)
                if lblstr:
                    path = os.path.join(path, "prob_" + lblstr)
                elif sharpen:
                    path = os.path.join(path, "templates")
                else:
                    path = os.path.join(path, "mean")
            path = os.path.join(path, self.timename(t) + ".nii.gz")
        if create and (force or not os.path.exists(path)):
            makedirs(os.path.dirname(path))
            if "labels" in cfg and not label:
                labels = parselabels(cfg["labels"])
                args = [self.avgimg(t, channel=channel, label=label, step=step, decomposed=decomposed, force=force, batch=batch) for label in labels]
                self._run("em-hard-segmentation", args=[len(args)] + args + [path])
            else:
                table = self.imgtable(t, step=step, channel=channel, decomposed=decomposed, force=force, batch=batch)
                opts = {
                    "images": table,
                    "reference": self.refimg
                }
                if "bkgrnd" in cfg:
                    opts["padding"] = float(cfg["bkgrnd"])
                opts["datatype"] = cfg.get("datatype", "float")
                if "labels" in cfg:
                    opts["label"] = label
                    if opts["datatype"] not in ["float", "double"]:
                        opts["rescaling"] = cfg.get("rescaling", [0, 100])
                    elif "rescaling" in cfg:
                        opts["rescaling"] = cfg["rescaling"]
                else:
                    opts["threshold"] = .5
                    opts["normalization"] = cfg.get("normalization", "zscore")
                    opts["rescaling"] = cfg.get("rescaling", [0, 100])
                    if sharpen:
                        opts["sharpen"] = cfg.get("sharpen", True)
                self._run("average-images", args=[path], opts=opts)
        return path

    def avgimgs(self, step=-1, ages=[], channels=[], labels={}, sharpen=True, outdir=None, decomposed=True, force=False, create=True, queue=None, batchname="avgimgs"):
        """Create all average images required for (parallel) atlas construction."""
        if step < 0:
            step = self.step
        if not queue:
            queue = self.queue["short"]
        if queue == "local":
            queue = None
        if ages:
            self.info("Average images at discrete time points", step=step)
            ages = [self.normtime(t) for t in ages]
        else:
            self.info("Average images at observed time points", step=step)
            ages = self.ages()
        single_channel = channels and isinstance(channels, basestring)
        if not channels:
            channels = self.regcfg(step).get("channels", self.channel)
        if not isinstance(channels, list):
            channels = [channels]
        if not isinstance(labels, dict):
            dlabels = {}
            for channel in channels:
                dlabels[channel] = labels
            labels = dlabels
        script = ""
        tasks = 0
        imgs = {}
        for channel in channels:
            imgs[channel] = {}
            for t in ages:
                imgs[channel][t] = []
        for channel in channels:
            segments = [0]
            if "labels" in self.config["images"][channel]:
                lbls = labels.get(channel, [])
                if lbls:
                    if isinstance(lbls, basestring):
                        if lbls.lower() == "all":
                            segments = parselabels(self.config["images"][channel]["labels"])
                        else:
                            segments = parselabels(lbls)
                    else:
                        segments = lbls
            for segment in segments:
                for t in ages:
                    img = self.avgimg(t, channel=channel, label=segment, sharpen=sharpen, outdir=outdir, step=step,
                                      decomposed=decomposed, force=force, create=create and not queue)
                    if queue and (force or not os.path.exists(img)):
                        remove_or_makedirs(img)
                        script += 'elif taskid == {taskid}: atlas.avgimg(t={t}, channel="{channel}", label={segment}, sharpen={sharpen}, outdir={outdir}, decomposed={decomposed}, batch=True)\n'.format(
                            taskid=tasks, t=t, channel=channel, segment=repr(segment), sharpen=sharpen, outdir=repr(outdir), decomposed=decomposed
                        )
                        tasks += 1
                    if len(segments) == 1:
                        imgs[channel][t] = img
                    else:
                        imgs[channel][t].append(img)
        if create and queue:
            job = self._submit(batchname, script=script, tasks=tasks, step=step, queue=queue)
        else:
            job = (None, 0)
        if single_channel:
            imgs = imgs[channels[0]]
        return (job, imgs)

    def construct(self, start=-1, niter=10, outdir=None, force=False, queue=None):
        """Perform atlas construction.

        Args:
            start (int): Last completed iteration. (default: step)
            niter (int, optional): Number of atlas construction iterations. (default: 10)
            outdir (str, optional): Directory for final templates. (default: config.paths.outdir)
            force (bool, optional): Force re-creation of already existing files. (default: False)
            queue (str, dict, optional): Name of queues of batch queuing system.
                When not specified, the atlas construction runs on the local machine
                using the number of threads specified during construction. When a single `str`
                is given, both short and long running jobs are submitted to the same queue.
                Otherwise, separate environments can be specified for "short" or "long" running
                jobs using the respective dictionary keys. The supported environments are:
                - "local": Multi-threaded execution on host machine
                - "condor": Batch execution using HTCondor
                - "<other>": Batch execution using named SLURM partition
                (default: config.environment.queue)

        """
        if start < 0:
            start = self.step
        if start < 0:
            raise ValueError("Atlas to be constructed must have step index >= 0!")
        if start > 0:
            self.info("Performing {0} iterations starting with step {1}".format(niter, start))
        else:
            self.info("Performing {0} iterations".format(niter))
        self.info("Age-dependent image deformations = {}".format(self.age_specific_imgdof))
        self.info("Average age-dependent deformations = {}".format(self.age_specific_regdof))
        # Initialize dict of queues used for batch execution
        if not queue:
            queue = self.queue
        if not isinstance(queue, dict):
            queue = {"short": queue, "long": queue}
        else:
            if "short" not in queue:
                queue["short"] = "local"
            if "long" not in queue:
                queue["long"] = "local"
        # Save considerable amount of disk memory by not explicitly storing the
        # composite image to atlas transformations of type FluidFreeFormTransformation
        # (alternatively, approximate composition). The transform-image and/or
        # average-images commands can apply a sequence of transforms in order to
        # perform the composition quasi on-the-fly.
        decomposed_imgdofs = True
        rmtemp_regdofs = True
        # Iterate atlas construction steps
        weights = {}
        for t in self.ages():
            weights[t] = self.weights(t)
        for step in range(start + 1, start + niter + 1):
            # Deform images to atlas space
            if not decomposed_imgdofs:
                job = self.imgdofs(step=step - 1, force=force, queue=queue["short"])[0]
                self.wait(job, interval=30, verbose=1)
            job = self.defimgs(step=step - 1, decomposed=decomposed_imgdofs, force=force, queue=queue["short"])[0]
            self.wait(job, interval=30, verbose=1)
            # Average images in atlas space
            job = self.avgimgs(step=step - 1, force=force, queue=queue["short"])[0]
            self.wait(job, interval=60, verbose=2)
            # Register all images to the current template images
            job = self.regdofs(step=step, force=force, queue=queue["long"])[0]
            self.wait(job, interval=60, verbose=5)
            # Compute all required average deformations
            job = self.avgdofs(step=step, force=force, queue=queue["short"])[0]
            self.wait(job, interval=60, verbose=2)
            if rmtemp_regdofs and self.age_specific_regdof:
                self.info("Deleting temporary deformation files", step=step)
                for t in weights:
                    for imgid in weights[t]:
                        dof1 = self.regdof(imgid, step=step, create=False)
                        dof2 = self.regdof(imgid, t=t, step=step, create=False)
                        if dof2 != "identity" and dof1 != dof2 and os.path.exists(dof2):
                            os.remove(dof2)
            # Compute all required longitudinal deformations
            if self.age_specific_regdof or self.age_specific_imgdof:
                job = self.compose(step=step, force=force, queue=queue["short"])[0]
                self.wait(job, interval=30, verbose=1)
        # Write final template images to specified directory
        self.step = start + niter
        self.info("Creating final mean shape templates")
        if outdir is None:
            outdir = self.outdir
        if outdir:
            ages = self.means
        else:
            ages = self.ages()
            outdir = None
        if not decomposed_imgdofs:
            job = self.imgdofs(ages=ages, force=force, queue=queue["short"])[0]
            self.wait(job, interval=30, verbose=1)
        channels = [channel for channel in self.config["images"].keys() if channel not in ("default", "ref")]
        job = self.defimgs(channels=channels, ages=ages, decomposed=decomposed_imgdofs, force=force, queue=queue["short"])[0]
        self.wait(job, interval=30, verbose=1)
        job = self.avgimgs(channels=channels, ages=ages, force=force, queue=queue["short"], outdir=outdir)[0]
        self.wait(job, interval=60, verbose=2)
        self.info("Finished atlas construction!")

    def evaluate(self, ages=[], step=-1, force=False, queue=None):
        """Evaluate atlas sharpness measures."""
        if not ages:
            ages = self.means
        if isinstance(ages, int):
            ages = [float(ages)]
        elif isinstance(ages, float):
            ages = [ages]
        if isinstance(step, int):
            if step < 0:
                if self.step < 0:
                    raise ValueError("Need to specify which step/iteration of atlas construction to evaluate!")
                step = self.step
            steps = [step]
        else:
            steps = step
        measures = self.config["evaluation"]["measures"]
        rois_spec = self.config["evaluation"].get("rois", {})
        roi_paths = {}
        roi_label = {}
        roi_labels = {}
        roi_channels = set()
        for roi_name, roi_path in rois_spec.items():
            labels = []
            if isinstance(roi_path, list):
                if len(roi_path) != 2:
                    raise ValueError("Invalid evaluation ROI value, must be either path (format) string of individual ROI or [<path_format>, <range>]")
                labels = roi_path[1]
                roi_path = roi_path[0]
                if roi_path in self.config["images"] and isinstance(labels, basestring) and labels.lower() == "all":
                    labels = parselabels(self.config["images"][roi_path].get("labels", ""))
                    if not labels:
                        raise ValueError("ROI channel {} must have 'labels' specified to use for 'all'!".format(roi_path))
                else:
                    labels = parselabels(labels)
                roi_name_format = roi_name
                if roi_name_format.format(l=0) == roi_name_format:
                    raise ValueError("Invalid evaluation ROI key name, name must include '{l}' format string!")
                for label in labels:
                    roi_name = roi_name_format.format(l=label)
                    roi_paths[roi_name] = roi_path
                    roi_label[roi_name] = label
            else:
                roi_paths[roi_name] = roi_path
            if roi_path in self.config["images"]:
                roi_channels.add(roi_path)
                roi_labels[roi_path] = labels
        roi_names = roi_paths.keys()
        roi_names.sort()
        roi_channels = list(roi_channels)
        re_measure = re.compile("^\s*(\w+)\s*\((.*)\)\s*$")
        # Evaluate voxel-wise measures
        for step in steps:
            voxelwise_measures = {}
            for t in ages:
                voxelwise_measures[t] = []
            for channel in measures:
                channel_info = self.config["images"][channel]
                channel_measures = measures[channel]
                if isinstance(channel_measures, basestring):
                    channel_measures = [channel_measures]
                if channel_measures:
                    # Deform individual images to atlas time point
                    name = "eval_defimgs_{channel}".format(channel=channel)
                    job, imgs = self.defimgs(channels=channel, ages=ages, step=step, queue=queue, batchname=name)
                    self.wait(job, interval=30, verbose=1)
                    # Evaluate measures for this image channel/modality
                    for measure in channel_measures:
                        measure = measure.lower().strip()
                        match = re_measure.match(measure)
                        if match:
                            measure = match.group(1)
                            args = match.group(2).strip()
                        else:
                            args = ""
                        # Evaluate gradient magnitude of average image
                        if measure == "grad":
                            sharpen = False
                            if "labels" in channel_info:
                                if not args:
                                    raise ValueError("Gradient magnitude of segmentation can only be computed for probability map of one or more label(s)!")
                                labels = [args]
                            else:
                                if args and args not in ("mean", "avg", "average", "template"):
                                    raise ValueError("Gradient magnitude of intensity images can only be computed from 'mean'/'avg'/'average'/'template' image!")
                                if args and args == "template":
                                    sharpen = True
                                labels = []
                            name = "eval_avgimgs_{channel}".format(channel=channel)
                            job, avgs = self.avgimgs(channels=channel, labels={channel: labels}, ages=ages, sharpen=sharpen, step=step, queue=queue)
                            self.wait(job, interval=60, verbose=2)
                            if args:
                                measure += "_" + args
                            for t in ages:
                                path = os.path.join(self.subdir(step), channel, measure, self.timename(t) + ".nii.gz")
                                if force or not os.path.exists(path):
                                    makedirs(os.path.dirname(path))
                                    self._run("detect-edges", step=step, queue=queue, wait_kwargs={"interval": 10, "verbose": 3},
                                              name="eval_{channel}_{measure}_{age}".format(channel=channel, measure=measure, age=self.timename(t)),
                                              args=[avgs[t], path], opts={"padding": channel_info.get("bkgrnd", -1), "central": None})
                                voxelwise_measures[t].append((channel, measure, path))
                        else:
                            if args:
                                raise ValueError("Measure '{}' has no arguments!".format(measure))
                            opts = {
                                "bins": channel_info.get("bins", 0),
                                "normalization": channel_info.get("normalization", "none")
                            }
                            if "rescaling" in channel_info:
                                opts["rescale"] = channel_info["rescaling"]
                            if "bkgrnd" in channel_info:
                                opts["padding"] = channel_info["bkgrnd"]
                            for t in ages:
                                path = os.path.join(self.subdir(step), channel, measure, self.timename(t) + ".nii.gz")
                                if force or not os.path.exists(path):
                                    makedirs(os.path.dirname(path))
                                    opts["output"] = os.path.relpath(path, self.topdir)
                                    mask = self.config["evaluation"].get("mask", "").format(t=t)
                                    if mask:
                                        opts["mask"] = mask
                                    self._run("aggregate-images", step=step, queue=queue, workdir=self.topdir, wait_kwargs={"interval": 30, "verbose": 1},
                                              name="eval_{channel}_{measure}_{age}".format(channel=channel, measure=measure, age=self.timename(t)),
                                              args=[measure] + [os.path.relpath(img, self.topdir) for img in imgs[t]], opts=opts)
                                voxelwise_measures[t].append((channel, measure, path))
            # Average voxel-wise measures within each ROI
            for channel in roi_channels:
                # Deform individual images to atlas time point
                name = "eval_defimgs_{channel}".format(channel=channel)
                job, imgs = self.defimgs(channels=channel, ages=ages, step=step, queue=queue, batchname=name)
                self.wait(job, interval=30, verbose=1)
                name = "eval_avgrois_{channel}".format(channel=channel)
                job, rois = self.avgimgs(channels=[channel], labels=roi_labels, ages=ages, step=step, queue=queue, batchname=name)
                self.wait(job, interval=60, verbose=2)
            for t in ages:
                if voxelwise_measures[t] and roi_paths:
                    subdir = self.subdir(step)
                    mean_table = os.path.join(subdir, "qc-measures", self.timename(t) + "-mean.csv")
                    sdev_table = os.path.join(subdir, "qc-measures", self.timename(t) + "-sdev.csv")
                    wsum_table = os.path.join(subdir, "qc-measures", self.timename(t) + "-wsum.csv")
                    if force or not os.path.exists(mean_table) or not os.path.exists(sdev_table) or not os.path.exists(wsum_table):
                        outdir = os.path.dirname(mean_table)
                        tmp_mean_table = os.path.join(outdir, "." + os.path.basename(mean_table))
                        tmp_sdev_table = os.path.join(outdir, "." + os.path.basename(sdev_table))
                        tmp_wsum_table = os.path.join(outdir, "." + os.path.basename(wsum_table))
                        args = [x[2] for x in voxelwise_measures[t]]
                        opts = {
                            "name": [],
                            "roi-name": [],
                            "roi-path": [],
                            "header": None,
                            "preload": None,
                            "mean": tmp_mean_table,
                            "sdev": tmp_sdev_table,
                            "size": tmp_wsum_table,
                            "digits": self.config["evaluation"].get("digits", 9)
                        }
                        for voxelwise_measure in voxelwise_measures[t]:
                            channel = voxelwise_measure[0]
                            measure = voxelwise_measure[1]
                            opts["name"].append("{}/{}".format(channel, measure))
                        for roi_name in roi_names:
                            roi_path = roi_paths[roi_name]
                            if roi_path in roi_channels:
                                roi_path = self.avgimg(t, channel=roi_path, label=roi_label.get(roi_name, 0), step=step, create=False)
                            else:
                                roi_path = roi_path.format(
                                    subdir=subdir, tmpdir=self.tmpdir, topdir=self.topdir,
                                    i=step, t=t, l=roi_label.get(roi_name, 0)
                                )
                            opts["roi-name"].append(roi_name)
                            opts["roi-path"].append(roi_path)
                        makedirs(outdir)
                        try:
                            self._run("average-measure", step=step, queue=queue,
                                      name="eval_average_{age}".format(age=self.timename(t)),
                                      args=args, opts=opts)
                        except Exception as e:
                            for tmp_table in [tmp_mean_table, tmp_sdev_table, tmp_wsum_table]:
                                if os.path.exists(tmp_table):
                                    os.remove(tmp_table)
                            raise e
                        cur_wait_time = 0
                        inc_wait_time = 10
                        max_wait_time = 6 * inc_wait_time
                        if queue and queue.lower() != "local":
                            while True:
                                missing = []
                                for tmp_table in [tmp_mean_table, tmp_sdev_table, tmp_wsum_table]:
                                    if not os.path.exists(tmp_table):
                                        missing.append(tmp_table)
                                if not missing:
                                    break
                                if cur_wait_time < max_wait_time:
                                    time.sleep(inc_wait_time)
                                    cur_wait_time += inc_wait_time
                                else:
                                    raise Exception("Job average-measure finished, but output files still missing after {}s: {}".format(cur_wait_time, missing))
                        for src, dst in zip([tmp_mean_table, tmp_sdev_table, tmp_wsum_table], [mean_table, sdev_table, wsum_table]):
                            try:
                                os.rename(src, dst)
                            except OSError as e:
                                sys.stderr.write("Failed to rename '{}' to '{}'".format(src, dst))
                                raise e

    def template(self, i, channel=None):
        """Get absolute path of i-th template image."""
        if i < 0 or i >= len(self.means):
            raise IndexError()
        if not channel:
            channel = self.channel
        img = self.avgimg(self.means[i], channel=channel, step=self.step, create=False, outdir=self.outdir)
        if self.step >= 0 and not os.path.exists(img):
            avg = self.avgimg(self.means[i], channel=channel, step=self.step, create=False)
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

    def deform(self, i, t, path=None, channel=None, force=False, create=True):
        """Deform i-th template using longitudinal deformations to time point t."""
        if not channel:
            channel = self.channel
        source = self.template(i, channel=channel)
        t = self.normtime(t)
        if t == self.means[i]:
            if path and os.path.realpath(os.path.abspath(path)) != os.path.realpath(source):
                if force and (create or not os.path.exists(path)):
                    shutil.copyfile(source, path)
                return path
            else:
                return source
        if not path:
            path = os.path.join(self.subdir(self.step), channel, "temp", "{}-{}.nii.gz".format(self.timename(self.means[i]), self.timename(t)))
        if create and (force or not os.path.exists(path)):
            dof = self.deformation(i, t)
            makedirs(os.path.dirname(path))
            self._run("transform-image", args=[source, path], opts={"dofin": dof, "target": self.refimg})
        return path

    def interpolate(self, t, path=None, channel=None, interp="default", deform=False, sigma=0, force=False, create=True):
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
        if not channel:
            channel = self.channel
        t = self.normtime(t)
        i = self.timeindex(t)
        if t == self.means[i] and interp.startswith("linear"):
            return self.template(i, channel=channel)
        if not path:
            path = os.path.join(self.outdir, self.timename(t) + ".nii.gz")
        if create and (force or not os.path.exists(path)):
            args = [path]
            # Linear interpolation
            if interp == "linear":
                w = (self.means[i + 1] - t) / (self.means[i + 1] - self.means[i])
                args.extend(["-image", self.template(i, channel=channel), w])
                if deform:
                    args.extend(["-dof", self.deformation(i, t)])
                w = (t - self.means[i]) / (self.means[i + 1] - self.means[i])
                args.extend(["-image", self.template(i + 1, channel=channel), w])
                if deform:
                    args.extend(["-dof", self.deformation(i + 1, t)])
            # Gaussian interpolation
            elif interp == "gaussian":
                for i in range(len(self.means)):
                    w = self.weight(self.means[i], mean=t, sigma=sigma)
                    if w > 0.:
                        args.extend(["-image", self.template(i, channel=channel), w])
                        if deform:
                            dof = self.deformation(self.means[i], t)
                            if dof != "identity":
                                args.extend(["-dof", dof])
            if create:
                bkgrnd = self.config["images"][channel].get("bkgrnd", -1)
                self._run("average-images", args=args, opts={"reference": self.refimg, "datatype": "uchar", "padding": bkgrnd})
        return path

    def view(self, i=None, channel=None):
        """View template image at specified time point(s)."""
        if i is None:
            i = range(len(self.means))
        elif not isinstance(i, (list, tuple)):
            i = [i]
        if not channel:
            channel = self.channel
        imgs = [self.template(int(idx), channel=channel) for idx in i]
        mirtk.run("view", target=imgs)

    def rmtemp(self):
        """Delete temporary files.

        This function removes all temporary files which can be recomputed
        without the need for performing the more costly registrations.
        Intermediate template images, auxiliary CSV files, and composite
        transformations are among these temporary files to be deleted.

        """
        raise NotImplementedError()


##############################################################################
# Auxiliaries

# ----------------------------------------------------------------------------
def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    """Compare too floating point numbers."""
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


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
def remove_or_makedirs(path):
    """Remove file when it exists or make directories otherwise."""
    if os.path.exists(path):
        os.remove(path)
    else:
        makedirs(os.path.dirname(path))


# ----------------------------------------------------------------------------
def parselabels(labels):
    """Parse labels specification."""
    values = []
    isstr = isinstance(labels, basestring)
    if isstr and "," in labels:
        labels = [arg.trim() for arg in labels.split(",")]
    if isinstance(labels, (tuple, list)):
        for label in labels:
            values.extend(parselabels(label))
    elif isstr:
        m = re.match("([0-9]+)..([0-9]+)", labels)
        if m:
            values.extend(range(int(m.group(1)), int(m.group(2)) + 1))
        elif ":" in labels:
            range_spec = [int(x) for x in labels.split(":")]
            if len(range_spec) == 2:
                values.extend(list(range(range_spec[0], range_spec[2] + 1, range_spec[1])))
            else:
                values.extend(list(range(range_spec[0], range_spec[1] + 1)))
        elif labels != "":
            values.append(int(labels))
    else:
        values.append(int(labels))
    return values


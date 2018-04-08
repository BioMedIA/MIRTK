#!/usr/bin/python

import os
import sys
import json
import argparse

from mirtk.atlas.spatiotemporal import SpatioTemporalAtlas


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
        Construct (spatio-temporal) atlas from images of the same anatomy of different subjects.

        To construct a spatio-temporal atlas, images of subjects at different ages spread over
        the atlas time range are required. Please cite the following preprint when you use this
        command in your research (or the respective peer-reviewed article when accepted):

        Schuh et al., "Unbiased construction of a temporally consistent morphological
        atlas of neonatal brain development", bioRxiv, 2018. doi:10.1101/251512""")
    parser.add_argument("config", help="JSON file with atlas configuration.")
    parser.add_argument("-a", "--ages", "--means", dest="means", type=float, nargs="+",
                        help="Discrete time points for which to construct atlas.")
    parser.add_argument("-s", "--sigma", dest="sigma", type=float, nargs="+", default=1.,
                        help="Standard deviation of temporal kernel(s).")
    parser.add_argument("-o", "--output", "--outdir", dest="outdir",
                        help="Output directory of final atlas.")
    parser.add_argument("-w", "--workdir", "--tmpdir", dest="tmpdir",
                        help="Working directory for intermediate files.")
    parser.add_argument("-e", "--energy", type=str,
                        help="Energy function of image registration (e.g., 'asym', sym', 'ic').")
    parser.add_argument("-c", "--channels", nargs="+", type=str,
                        help="Images to use for multi-channel registration.")
    parser.add_argument("-m", "--measures", nargs="+", type=str,
                        help="Image (dis-)similarity measure to use for each channel.")
    parser.add_argument("-b", "--bending", nargs="+", type=float,
                        help="Bending energy weight for each step.")
    parser.add_argument("--bins", type=int,
                        help="No. of bins to use for histogram-based measures.")
    parser.add_argument("--window", type=int,
                        help="Local window size for NCC measure in voxels.")
    parser.add_argument("--spacing", dest="spacing", nargs="+", type=float,
                        help="Control point spacing for each step.")
    parser.add_argument("-j", "--jacobian", nargs="+", type=float,
                        help="Jacobian penalty weight for each step.")
    parser.add_argument("--bch-terms", dest="bchterms", type=int,
                        help="No. of BCH terms to use for composition of SV FFDs.")
    parser.add_argument("-g", "--growth", type=bool,
                        help="Whether to joinedly estimate mean shape and change.")
    parser.add_argument("-i", "--start", dest="start", default=0, type=int,
                        help="Continue construction after the specified step.")
    parser.add_argument("-n", "--steps", dest="steps", default=-1, type=int,
                        help="Number of steps of iterative atlas construction/refinement. Use 'iterations' config entry if not specified (default: 10).")
    parser.add_argument("-q", "--queue", "--long-queue", dest="longqueue",
                        help="Name of batch system queue. Use 'condor' for HTCondor. Otherwise, the argument is assumed to be the name of a SLURM partition.")
    parser.add_argument("--short-queue", dest="shortqueue",
                        help="Name of batch system queue to use for short running jobs (about 1-30 min). Use --long-queue by default.")
    parser.add_argument("-t", "--threads", type=int,
                        help="Maximum number of CPU cores/threads to use.")
    parser.add_argument("-v", "--verbose", default=1, type=int,
                        help="Verbosity level of output messages: 0) no output, 1) report progress, 2) print command arguments.")
    #parser.add_argument("--debug", default=2, type=int,
    #                    help="Set debugging level: 0) delete all temp file, 1) keep atlas at each iteration, 2) keep all.")
    args = parser.parse_args()
    args.config = os.path.abspath(args.config)
    root = os.path.dirname(args.config)
    with open(args.config, "rt") as f:
        config = json.load(f)
    # paths
    if "paths" not in config:
        config["paths"] = {
            "topdir": os.getcwd()
        }
    if args.outdir:
        config["paths"]["outdir"] = os.path.abspath(args.outdir)
    elif args.outdir is not None:
        config["paths"]["outdir"] = ""
    if args.tmpdir:
        config["paths"]["tmpdir"] = os.path.abspath(args.tmpdir)
    # no. of iterations
    if args.steps < 0:
        args.steps = config.get("iterations", 10)
    # registration parameters
    if "registration" not in config:
        config["registration"] = {"config": [{}], "growth": {}}
    elif "config" not in config["registration"]:
        config["registration"]["config"] = [{}]
    elif not isinstance(config["registration"]["config"], list):
        config["registration"]["config"] = [config["registration"]["config"]]
    regcfg = config["registration"]["config"]
    for cfg in regcfg:
        if args.energy:
            cfg["energy"] = args.energy
        if args.channels:
            cfg["channels"] = args.channels
        if args.measures:
            cfg["measures"] = args.measures
        if args.bins:
            cfg["bins"] = args.bins
        if args.window:
            cfg["window"] = args.window
    if args.spacing:
        for i in range(max(len(regcfg), len(args.spacing))):
            if i >= len(regcfg):
                regcfg.append({})
            regcfg[i]["spacing"] = args.spacing[i] if i < len(args.spacing) else args.spacing[-1]
    if args.bending:
        for i in range(max(len(regcfg), len(args.bending))):
            if i >= len(regcfg):
                regcfg.append({})
            regcfg[i]["bending"] = args.bending[i] if i < len(args.bending) else args.bending[-1]
    if args.jacobian:
        for i in range(max(len(regcfg), len(args.jacobian))):
            if i >= len(regcfg):
                regcfg.append({})
            regcfg[i]["jacobian"] = args.jacobian[i] if i < len(args.jacobian) else args.jacobian[-1]
    # longitudinal growth modeling
    if "growth" not in config["registration"]:
        config["registration"]["growth"] = {}
    if args.bchterms:
        config["registration"]["growth"]["bchterms"] = args.bchterms
    if args.growth:
        config["registration"]["growth"]["enabled"] = args.growth
    # regression kernels
    if "regression" not in config:
        config["regression"] = {}
    if args.means:
        config["regression"]["means"] = args.means
        config["regression"]["sigma"] = args.sigma
    if args.verbose > 3:
        json.dump(config, sys.stdout, indent=4, separators=(',', ': '))
        sys.stdout.write("\n")
        sys.stdout.flush()
    if "environment" not in config:
        config["environment"] = {}
    if "queue" not in config["environment"]:
        config["environment"]["queue"] = {"short": "local", "long": "local"}
    if args.threads:
        config["environment"]["threads"] = args.threads
    if not args.shortqueue and args.longqueue:
        args.shortqueue = args.longqueue
    if args.shortqueue:
        config["environment"]["queue"]["short"] = args.shortqueue
    if args.longqueue:
        config["environment"]["queue"]["long"] = args.longqueue
    atlas = SpatioTemporalAtlas(config=config, root=root, verbose=args.verbose, exit_on_error=True)
    atlas.construct(start=args.start, niter=args.steps - args.start)

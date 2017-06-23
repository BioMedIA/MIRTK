#!/usr/bin/python

import os
import sys
import json
import argparse

from mirtk.atlas.spatiotemporal import SpatioTemporalAtlas


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Evaluate (spatio-temporal) atlas sharpness measures.""")
    parser.add_argument("config", help="JSON file with atlas configuration.")
    parser.add_argument("-a", "--ages", "--means", dest="means", type=float, nargs="+",
                        help="Discrete time points for which to perform evaluation.")
    parser.add_argument("-w", "--workdir", "--tmpdir", dest="tmpdir",
                        help="Working directory of intermediate files.")
    parser.add_argument("-i", "--step", dest="step", default=-1, type=int, nargs="+",
                        help="Atlas construction steps to evaluate. Evaluate all 'iterations' specified in config otherwise (default: 10).")
    parser.add_argument("-q", "--queue", dest="queue",
                        help="Name of batch system queue. Use 'condor' for HTCondor. Otherwise, the argument is assumed to be the name of a SLURM partition.")
    parser.add_argument("-t", "--threads", type=int,
                        help="Maximum number of CPU cores/threads to use.")
    parser.add_argument("-v", "--verbose", default=1, type=int,
                        help="Verbosity level of output messages: 0) no output, 1) report progress, 2) print command arguments.")
    args = parser.parse_args()
    args.config = os.path.abspath(args.config)
    root = os.path.dirname(args.config)
    with open(args.config, "rt") as f:
        config = json.load(f)
    if "paths" not in config:
        config["paths"] = {
            "topdir": os.getcwd()
        }
    if args.tmpdir:
        config["paths"]["tmpdir"] = os.path.abspath(args.tmpdir)
    if "environment" not in config:
        config["environment"] = {}
    if args.threads:
        config["environment"]["threads"] = args.threads
    if args.step == -1:
        args.step = list(range(1, config.get("iterations", 10) + 1))
    atlas = SpatioTemporalAtlas(root=root, config=config, verbose=args.verbose, exit_on_error=True)
    atlas.evaluate(step=args.step, ages=args.means, queue=args.queue)

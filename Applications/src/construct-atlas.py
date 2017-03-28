#!/usr/bin/python

import os
import argparse

from mirtk.atlas.spatiotemporal import SpatioTemporalAtlas


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
        Construct (spatio-temporal) atlas from images of the same anatomy of different subjects.
        To construct a spatio-temporal atlas, images of subjects at different ages spread over
        the atlas time range are required.""")
    parser.add_argument("config", help="JSON file with atlas configuration.")
    parser.add_argument("-i", "--start", dest="start", default=0, type=int,
                        help="Continue construction after the specified step.")
    parser.add_argument("-n", "--steps", dest="steps", default=10, type=int,
                        help="Number of steps of iterative atlas construction/refinement.")
    parser.add_argument("-t", "--threads", default=8, type=int,
                        help="Maximum number of CPU cores/threads to use.")
    parser.add_argument("-q", "--queue", default=None,
                        help="Name of batch system queue. Use 'condor' for HTCondor. Otherwise, the argument is assumed to be the name of a SLURM partition.")
    parser.add_argument("--short-queue", default=None,
                        help="Name of batch system queue to use for short running jobs (about 1-30 min).")
    parser.add_argument("-v", "--verbose", default=1, type=int,
                        help="Verbosity level of output messages: 0) no output, 1) report progress, 2) print command arguments.")
    #parser.add_argument("--debug", default=2, type=int,
    #                    help="Set debugging level: 0) delete all temp file, 1) keep atlas at each iteration, 2) keep all.")
    args = parser.parse_args()
    if not args.short_queue:
        args.short_queue = args.queue
    atlas = SpatioTemporalAtlas(args.config, threads=args.threads, verbose=args.verbose)
    atlas.construct(start=args.start, niter=args.steps - args.start, queue=[args.short_queue, args.queue])

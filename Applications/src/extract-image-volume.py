#!/usr/bin/env python

"""
Split image sequence into individual volumes.
"""

import mirtk
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help="Input sequence.")
    parser.add_argument('output', help="Output file path.")
    parser.add_argument('-t', type=int, default=0, help="Index of volume to extract.")
    parser.add_argument('-n', type=int, default=0,
                        help="Number of volumes to extract. Extract all from -t to maximum index when non-positive.")
    args = parser.parse_args()
    opts = {"Rt1": args.t}
    if args.n > 0:
        opts["Rt2"] = args.t + args.n - 1
    mirtk.extract_image_region(args.input, args.output, split='t', **opts)

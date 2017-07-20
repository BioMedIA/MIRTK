#!/usr/bin/env python

"""
Split image volume into individual slices.
"""

import mirtk
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help="Input volume")
    parser.add_argument('output', help="Output file path")
    parser.add_argument('-z', type=int, default=0, help="Index of slice to extract")
    parser.add_argument('-n', type=int, default=0,
                        help="Number of slices to extract. Extract all from -z to maximum index when non-positive")
    args = parser.parse_args()
    opts = {"Rz1": args.z}
    if args.n > 0:
        opts["Rz2"] = args.z + args.n - 1
    mirtk.extract_image_region(args.input, args.output, split='z', **opts)

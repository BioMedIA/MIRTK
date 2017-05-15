#!/usr/bin/python

import os
import sys
import numpy as np
import nibabel as nib
import argparse

def read_image_data(path):
    """Read image data into 1D numpy array of type float64."""
    return nib.load(path).get_data().flatten().astype(np.float64)


def filename(path):
    """Get filename without extension."""
    name = os.path.basename(path)
    name, ext = os.path.splitext(name)
    if ext.lower() == ".gz":
        name = os.path.splitext(name)[0]
    return name


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Average voxel-wise measure within each region of interest (ROI).")
    parser.add_argument("measure", nargs="+",
                        help="Input image(s) of voxel-wise scalar measure.")
    parser.add_argument("-n", "-name", "-names", "--name", "--names", dest="name", nargs="+", default=[],
                        help="Names of input measures.")
    parser.add_argument("-roi-name", "-roi-names", "--roi-name", "--roi-names", dest="roi_name", nargs="+", default=[],
                        help="Names of ROIs.")
    parser.add_argument("-o", "-output", "-table", "--output", "--table", dest="table",
                        help="Output table")
    parser.add_argument("-r", "-region", "-roi", "-roi-path", "--region", "--roi", "--roi-path", dest="roi", nargs="+",
                        help="Binary mask or segmentation probability map of ROI(s).")
    parser.add_argument("-d", "-delim", "-delimiter", "--delim", "--delimiter", dest="delim", default=",",
                        help="String/character used to delimit table entries.")
    parser.add_argument("-s", "-scale", "--scale", dest="scale", type=float, default=1.,
                        help="Common scaling factor by which to multiply input values.")
    parser.add_argument("-preload", "--preload", dest="preload", action="store_true",
                        help="Read each input image only once. Requires enough memory to preload all measures into memory.")
    parser.add_argument("-transpose", "--transpose", dest="transpose", action="store_true",
                        help="Transpose output table s.t. each row corresponds to an input measure")
    parser.add_argument("-header", "--header", dest="header", action="store_true",
                        help="Print table header")
    parser.add_argument("-id", "--id", dest="id_name",
                        help="Header value of ID column")
    parser.add_argument("-noid", "--noid", dest="noid", action="store_true",
                        help="Do not print ID column")
    parser.add_argument("-a", "-append", "--append", action="store_true",
                        help="Append results to existing table.")
    args = parser.parse_args()
    if len(args.name) > len(args.measure):
        raise ValueError("More --names specified than input measures given!")
    for i in range(len(args.name), len(args.measure)):
        args.name.append(filename(args.measure[i]))
    if not args.roi_name:
        args.roi_name = [""] * len(args.roi)
        for i in range(len(args.roi)):
            #args.roi_name[i] = "{0:d}".format(i + 1)
            args.roi_name[i] = filename(args.roi[i])
    out = sys.stdout
    sep = args.delim
    col = {}
    if args.preload:
        if args.transpose:
            for name, path in zip(args.roi_name, args.roi):
                col[name] = read_image_data(path)
        else:
            for name, path in zip(args.name, args.measure):
                col[name] = read_image_data(path)
    # FIXME: Read existing table and insert missing/newly computed values
    #        Consider pandas DataFrame to easily be able to refer to columns by measure name
    if args.table:
        if args.append:
            out = open(args.table, "at")
        else:
            out = open(args.table, "wt")
    if not args.id_name:
        args.id_name = "measure" if args.transpose else "roi"
    try:
        if args.header and (not args.table or not args.append):
            c = 0
            if not args.noid:
                out.write(args.id_name)
                c += 1
            if args.transpose:
                for roi_name in args.roi_name:
                    c += 1
                    if c > 1:
                        out.write(sep)
                    out.write(roi_name)
            else:
                for name in args.name:
                    c += 1
                    if c > 1:
                        out.write(sep)
                    out.write(name)
            out.write("\n")
        if args.transpose:
            for name, path in zip(args.name, args.measure):
                c = 0
                if not args.noid:
                    c += 1
                    out.write(name)
                values = read_image_data(path)
                for roi_name, roi_path in zip(args.roi_name, args.roi):
                    c += 1
                    if c > 1:
                        out.write(sep)
                    if roi_name in col:
                        roi = col[roi_name]
                    else:
                        roi = read_image_data(roi_path)
                    # Note: When using float32, the result differs from MIRTK calculate-element-wise and builtin sum
                    out.write("{0:.5f}".format(args.scale * np.sum(values * (roi / np.sum(roi)))))
                out.write("\n")
        else:
            for roi_name, roi_path in zip(args.roi_name, args.roi):
                c = 0
                if not args.noid:
                    c += 1
                    out.write(roi_name)
                roi = read_image_data(roi_path)
                for name, path in zip(args.name, args.measure):
                    c += 1
                    if c > 1:
                        out.write(sep)
                    if name in col:
                        values = col[name]
                    else:
                        values = read_image_data(path)                    
                    # Note: When using float32, the result differs from MIRTK calculate-element-wise and builtin sum
                    out.write("{0:.5f}".format(args.scale * np.sum(values * (roi / np.sum(roi)))))
                out.write("\n")
    finally:
        out.close()

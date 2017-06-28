#!/usr/bin/python

import os
import sys
import numpy as np
import nibabel as nib
import argparse


def read_image_data(path):
    """Read image data into 1D numpy array of type float64.

    Note: When using float32, the mean values differ between MIRTK calculate-element-wise and Python sum
    """
    return nib.load(path).get_data().flatten().astype(np.float64)


def filename(path):
    """Get filename without extension."""
    name = os.path.basename(path)
    name, ext = os.path.splitext(name)
    if ext.lower() == ".gz":
        name = os.path.splitext(name)[0]
    return name


def open_table(table_name, append=False, default=None):
    """Open output file if file name specified."""
    if not table_name:
        return default
    return open(table_name, "at" if append else "wt")


def close_table(out):
    """Close table stream if open."""
    if out is not None:
        out.close()


def write(out, value, number=False, digits=5):
    """Write string or floating point value to output stream if open."""
    if number:
        value = "{0:.{1}f}".format(value, digits)
    if isinstance(out, (tuple, list)):
        for f in out:
            write(f, value)
    elif out is not None:
        out.write(value)


def write_header(out, args):
    """Write table header to output stream."""
    if isinstance(out, (tuple, list)):
        for f in out:
            write_header(f, args)
    elif out is not None:
        c = 0
        if not args.noid:
            out.write(args.id_name)
            c += 1
        if args.transpose:
            for roi_name in args.roi_name:
                c += 1
                if c > 1:
                    out.write(args.delim)
                out.write(roi_name)
        else:
            for name in args.name:
                c += 1
                if c > 1:
                    out.write(args.delim)
                out.write(name)
        out.write("\n")


def write_stats(out, x, w, sum_w=0., sum_w2=0., unbiased=False, digits=5):
    """Write measurement statistics to output table(s)."""
    if sum_w <= 0.:
        sum_w = np.sum(w)
    mean = np.sum(w * x) / sum_w
    write(out[0], mean, number=True, digits=digits)
    if out[1] is not None or out[2] is not None or out[3] is not None:
        var = np.sum(w * np.square(x - mean))
        if unbiased:
            if sum_w2 <= 0.:
                sum_w2 = np.sum(np.square(w))
            var = var / (sum_w - (sum_w2 / sum_w))
        else:
            var = var / sum_w
        sdev = np.sqrt(var)
        write(out[1], var, number=True, digits=digits)
        write(out[2], sdev, number=True, digits=digits)
        write(out[3], sdev / np.sqrt(sum_w), number=True, digits=digits)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Average voxel-wise measure within each region of interest (ROI).")
    parser.add_argument("measure", nargs="+",
                        help="Input image(s) of voxel-wise scalar measure.")
    parser.add_argument("-n", "-name", "-names", "--name", "--names", dest="name", nargs="+", default=[],
                        help="Names of input measures. Input file name excluding file name extension is used by default.")
    parser.add_argument("-roi-name", "-roi-names", "--roi-name", "--roi-names", dest="roi_name", nargs="+", default=[],
                        help="Names of ROIs. By default, one-based ROI IDs are generated.")
    parser.add_argument("-size", "-roi-size", "--size", "--roi-size", dest="size_table",
                        help="Write sizes (sum of weights) of ROIs to specified file.")
    parser.add_argument("-o", "-mean", "-output", "-table", "--mean", "--output", "--table", dest="mean_table",
                        help="Write mean values to specified table")
    parser.add_argument("-sdev", "-stdev", "-stddev", "-standard-deviation", "--sdev", "--stdev", "--stddev", "--standard-deviation", dest="sdev_table",
                        help="Write standard deviation of weighted measurements to specified file.")
    parser.add_argument("-var", "-variance", "--var", "--variance", dest="var_table",
                        help="Write variance of weighted measurements to specified file.")
    parser.add_argument("-se", "-standard-error", "--se", "--standard-error", dest="se_table",
                        help="Write standard error of weighted mean to specified file.")
    parser.add_argument("-u", "-unbiased", "--unbiased", dest="unbiased", action='store_true',
                        help="When computing second moments, apply Bessel's correction.")
    parser.add_argument("-r", "-region", "-roi", "-roi-path", "--region", "--roi", "--roi-path", dest="roi", nargs="+", required=True,
                        help="Binary mask or segmentation probability map of ROI(s).")
    parser.add_argument("-d", "-delim", "-delimiter", "--delim", "--delimiter", dest="delim", default=",",
                        help="String/character used to delimit table entries.")
    parser.add_argument("-s", "-scale", "--scale", dest="scale", type=float, default=1.,
                        help="Common scaling factor by which to multiply input values.")
    parser.add_argument("-digits", "-precision", "--digits", "--precision", dest="digits", default=5,
                        help="Number of digits after the decimal point.")
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
            args.roi_name[i] = "{0:d}".format(i + 1)
    col = {}
    if args.preload:
        if args.transpose:
            for name, path in zip(args.roi_name, args.roi):
                roi = read_image_data(path)
                roi = roi / np.max(roi)
                col[name] = roi
        else:
            for name, path in zip(args.name, args.measure):
                col[name] = args.scale * read_image_data(path)
    if not args.id_name:
        args.id_name = "measure" if args.transpose else "roi"
    f_mean = None
    f_var = None
    f_sdev = None
    f_se = None
    f_size = None
    f_out = []
    try:
        f_mean = open_table(args.mean_table, args.append, default=sys.stdout)
        f_sdev = open_table(args.sdev_table, args.append)
        f_var = open_table(args.var_table, args.append)
        f_se = open_table(args.se_table, args.append)
        f_size = open_table(args.size_table, args.append)
        f_out = [f_mean, f_var, f_sdev, f_se]  # order must be same as expected by write_stats!
        if args.header:
            write_header(f_out, args)
            if args.transpose:
                write_header(f_size, args)
            elif args.noid:
                if f_var is not None or f_sdev is not None or f_se is not None:
                    write(f_size, "n" + args.delim + "n2" + "\n")
                else:
                    write(f_size, "n\n")
            else:
                if f_var is not None or f_sdev is not None or f_se is not None:
                    write(f_size, args.id_name + args.delim + "n" + args.delim + "n2" + "\n")
                else:
                    write(f_size, args.id_name + args.delim + "n\n")
        if args.transpose:
            if f_size is not None:
                rows = 1
                if f_var is not None or f_sdev is not None or f_se is not None:
                    rows = 2
                for r in range(1, rows + 1):
                    c = 0
                    if not args.noid:
                        c += 1
                        if r == 1:
                            write(f_size, "n")
                        elif r == 2:
                            write(f_size, "n2")
                    for roi_name, roi_path in zip(args.roi_name, args.roi):
                        c += 1
                        if c > 1:
                            write(f_size, args.delim)
                        if roi_name in col:
                            roi = col[roi_name]
                        else:
                            roi = read_image_data(roi_path)
                            roi = roi / np.max(roi)
                        if r == 1:
                            sum_w = np.sum(roi)
                            write(f_size, sum_w, number=True, digits=args.digits)
                        elif r == 2:
                            sum_w2 = np.sum(np.square(roi))
                            write(f_size, sum_w2, number=True, digits=args.digits)
                    write(f_size, "\n")
            for name, path in zip(args.name, args.measure):
                c = 0
                if not args.noid:
                    c += 1
                    write(f_out, name)
                values = args.scale * read_image_data(path)
                for roi_name, roi_path in zip(args.roi_name, args.roi):
                    c += 1
                    if c > 1:
                        write(f_out, args.delim)
                    if roi_name in col:
                        roi = col[roi_name]
                    else:
                        roi = read_image_data(roi_path)
                        roi = roi / np.max(roi)
                    write_stats(f_out, x=values, w=roi, unbiased=args.unbiased, digits=args.digits)
                write(f_out, "\n")
        else:
            for roi_name, roi_path in zip(args.roi_name, args.roi):
                c = 0
                if not args.noid:
                    c += 1
                    write(f_out, roi_name)
                    write(f_size, roi_name + args.delim)
                roi = read_image_data(roi_path)
                roi = roi / np.max(roi)
                sum_w = np.sum(roi)
                sum_w2 = 0.
                if f_size is not None:
                    write(f_size, sum_w, number=True, digits=args.digits)
                    if f_var is not None or f_sdev is not None or f_se is not None:
                        sum_w2 = np.sum(np.square(roi))
                        write(f_size, args.delim)
                        write(f_size, sum_w2, number=True, digits=args.digits)
                for name, path in zip(args.name, args.measure):
                    c += 1
                    if c > 1:
                        write(f_out, args.delim)
                    if name in col:
                        values = col[name]
                    else:
                        values = args.scale * read_image_data(path)
                    write_stats(f_out, x=values, w=roi, sum_w=sum_w, sum_w2=sum_w2, unbiased=args.unbiased, digits=args.digits)
                write(f_out, "\n")
                write(f_size, "\n")
    finally:
        close_table(f_mean)
        close_table(f_var)
        close_table(f_sdev)
        close_table(f_se)
        close_table(f_size)

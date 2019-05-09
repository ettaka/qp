#!env python

import numpy as np
import matplotlib.pyplot as plt
import time
#import pandas as pd
import argparse
import codecs
import scipy.interpolate

from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


fig = plt.figure()
ax = fig.add_subplot(111)
ax.autoscale(enable=True, axis='y', tight=True)

def read_coil_size_data(filepath):
    coil_size_dict = {}
    header_lines = 10
    with codecs.open(filepath) as f:
        head = [next(f).decode('utf-8', 'ignore') for x in xrange(header_lines)]

    for line in head:
        if 'shim' in line: shim = float(line.split()[1])
    print "shim:", shim

    col_names = ['L+R' if 'L+R' in s else s for s in head[-1].strip('\r\n').split('\t')]
    print col_names
    coil_size_dict['column_names'] = col_names
    coil_size_dict['column_indices'] = {}
    for i,col_name in enumerate(col_names):
        coil_size_dict['column_indices'][col_name] = i

    with codecs.open(filepath) as f:
        raw_data = np.loadtxt(f, skiprows=header_lines)

    coil_size_dict['raw_data'] = raw_data
    print raw_data, col_names
    return coil_size_dict

def add_interp_to_coil_dict(coil_size_dict):
    col_inds = coil_size_dict['column_indices']

    row_dict['interp'] = scipy.interpolate.InterpolatedUnivariateSpline(x_vec, row_dict['values'],k=1)

def plot_coil_size(coil_size_dict):
    xdata = coil_size_dict['raw_data'][:,0]
    ydata = coil_size_dict['raw_data'][:,1]
    label = coil_size_dict['column_names'][1]

    ax.plot(xdata,ydata,'--ro',label=label)
    plt.show()

def plot_coil_sizes(args):
    for filepath in args.paths:
        coil_size_dict = read_coil_size_data(filepath)
        plot_coil_size(coil_size_dict)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot coil sizes')
    parser.add_argument('paths', nargs='+', type=str)

    args = parser.parse_args()
    plot_coil_sizes(args)

    
    

 

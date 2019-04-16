#!env python

import numpy as np
import matplotlib.pyplot as plt
import time
#import pandas as pd
import argparse
import codecs

from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


fig = plt.figure()
ax = fig.add_subplot(111)
ax.autoscale(enable=True, axis='y', tight=True)

def make_error_boxes(ax, xdata, ydata, xerror, yerror, facecolor='r',
                     edgecolor='None', alpha=0.5):

    # Create list for all the error patches
    errorboxes = []

    # Loop over data points; create box from errors at each point
    for x, y, xe, ye in zip(xdata, ydata, xerror.T, yerror.T):
        rect = Rectangle((x - xe[0], y - ye[0]), xe.sum(), ye.sum())
        errorboxes.append(rect)

    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes, facecolor=facecolor, alpha=alpha,
                         edgecolor=edgecolor)

    # Add collection to axes
    ax.add_collection(pc)

    # Plot errorbars
    artists = ax.errorbar(xdata, ydata, xerr=xerror, yerr=yerror,
                          fmt='None', ecolor='k')

    return artists

def compute_data_dict_avg_min_max(data_dict):
    use_cols = []
    for i in range(np.shape(data_dict['raw_data'])[1]):
        if np.any(data_dict['raw_data'][:,i]!=0): use_cols.append(i)

    data_dict['average'] = np.average(data_dict['raw_data'][:,use_cols], axis=1)
    data_dict['min'] = np.min(data_dict['raw_data'][:,use_cols], axis=1)
    data_dict['max'] = np.max(data_dict['raw_data'][:,use_cols], axis=1)
    data_dict['error'] = np.array([data_dict['average'] - data_dict['min'], data_dict['max'] - data_dict['average']])

def create_pk_npk_dict(pk_npk_data_file):
    header_lines = 3
    with codecs.open(pk_npk_data_file) as f:
        pk_npk_raw_data = np.loadtxt(f, skiprows=header_lines)
    pk_npk_dict = {}
    pk_npk_dict['raw_data'] = pk_npk_raw_data
    pk_npk_dict['pk-spole'] = pk_npk_raw_data[:,1]
    pk_npk_dict['pk-scyl'] = pk_npk_raw_data[:,3]
    pk_npk_dict['npk-spole'] = pk_npk_raw_data[:,5]
    pk_npk_dict['npk-scyl'] = pk_npk_raw_data[:,7]
    return pk_npk_dict

def create_data_dicts(filepath, coil_permutation=None):
    print "Creating data dictionaries."
    if coil_permutation == None: coil_permutation = [1,2,3,4]
    print "Permuting coils with", coil_permutation
    coil_permutation = [inx-1 for inx in coil_permutation]
    header_lines = 3
    with codecs.open(filepath) as f:
        head = [next(f).decode('utf-8', 'ignore') for x in xrange(header_lines)]
    col_names = head[2].strip('\r\n').split('\t')

    with codecs.open(filepath) as f:
        raw_data = np.loadtxt(f, skiprows=header_lines)

    xdict = {}
    xdict['raw_data'] = raw_data[:,1:5]
    xdict['col_names'] = col_names[1:5]
    ydict = {}
    ypicker = map([5,6,7,8].__getitem__,coil_permutation)
    ydict['raw_data'] = raw_data[:,ypicker]
    ydict['col_names'] = map(col_names[5:9].__getitem__,coil_permutation)
    compute_data_dict_avg_min_max(xdict)
    compute_data_dict_avg_min_max(ydict)
    return xdict, ydict

def plot_tf(filepath, args):

    tr_type = args.type
    pk_npk_file = args.pk_npk_file
    single_coils = args.single_coils
    no_plot_average = args.no_average
    no_plot_average_error = args.no_average_error
    legend_location = args.legend_location
    neighbour_shell_averages = not args.no_neighbour_shell_averages
    coil_permutation = args.coil_permutation

    filebase = filepath.split('.txt')[0]
    save_fig = filebase + '.png'

    xdict, ydict = create_data_dicts(filepath, coil_permutation)
    xdata = xdict['average']
    ydata = ydict['average']

    pk_npk_dict = create_pk_npk_dict(pk_npk_file)

    colors = ['b','g','r','c','m','y','k']
    markers = ['o', '^', 'v', '<', '>', '1', '2', '3', '4', 's', 'p', '*', 'h', '+', 'x']

    if no_plot_average:
        print "Plot average of all shell vs pole gauges."
        ax.plot(xdata,ydata,'--ro',label='Meas. Avg.')
        if no_plot_average_error:
            _ = make_error_boxes(ax, xdata, ydata, xdict['error'], ydict['error'], facecolor='g', edgecolor='None', alpha=0.5)

    if single_coils:
        print "Plot single pole gauges",
        for i in range(4):
            if neighbour_shell_averages:
                if i==0:print "vs single shell gauges."
                xdata = xdict['raw_data'][:,i]
                ydata = ydict['raw_data'][:,i]
                label = xdict['col_names'][i]+'-'+ydict['col_names'][i]
            else:
                if i==0:print "vs single shell gauges."
                nof_cols = np.shape(xdict['raw_data'])[1]
                xdata = (xdict['raw_data'][:,i] + xdict['raw_data'][:,(i+1)%nof_cols])/2.
                ydata = ydict['raw_data'][:,i]
                label = xdict['col_names'][i]+xdict['col_names'][(i+1)%nof_cols]+'-'+ydict['col_names'][i]

            ax.plot(xdata, ydata,'--'+colors[i]+markers[i],label=label)

    ax.plot(pk_npk_dict['pk-scyl'],pk_npk_dict['pk-spole'],'-bo',label='FEM PK')
    ax.plot(pk_npk_dict['npk-scyl'],pk_npk_dict['npk-spole'],'-bd',label='FEM NPK')
    ax.set_xlabel('Shell Azimuthal Stress (MPa)')
    ax.set_ylabel('Pole Azimuthal Stress (MPa)')
    ax.grid()

    ax.legend(loc=legend_location)
    if args.show_plot:
        plt.show()
    else:
        plt.savefig(save_fig)

    #lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5,-.1), fancybox=True, shadow=True, ncol=2)
    #plt.savefig(save_fig, bbox_extra_artists=(lgd,), bbox_inches='tight')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot transfer function')
    parser.add_argument('paths', nargs='+', type=str)
    parser.add_argument('-t', '--type', type=int, default=1) 
    parser.add_argument('-pk', '--pk-npk-file', type=str, default='TRANSFER1_PK_NPK_simple.txt') 
    parser.add_argument('-s', '--single-coils', action='store_true', default=False) 
    parser.add_argument('-sp', '--show-plot', action='store_true', default=False) 
    parser.add_argument('-nav', '--no-average', action='store_false', default=True) 
    parser.add_argument('-nerr', '--no-average-error', action='store_false', default=True) 
    parser.add_argument('-nshav', '--no-neighbour-shell-averages', action='store_false', default=True) 
    parser.add_argument('-ll', '--legend-location', type=str, default='best')
    parser.add_argument('-perm', '--coil-permutation', nargs=4, type=int)

    args = parser.parse_args()
    paths = args.paths
    for filepath in paths:
        plt.cla()
        plot_tf(filepath, args)
    


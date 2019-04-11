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
    data_dict['average'] = np.average(data_dict['raw_data'], axis=1)
    data_dict['min'] = np.min(data_dict['raw_data'], axis=1)
    data_dict['max'] = np.max(data_dict['raw_data'], axis=1)
    print data_dict['average']
    print data_dict['min']
    print data_dict['max']
    data_dict['error'] = np.array([data_dict['average'] - data_dict['min'], data_dict['max'] - data_dict['average']])
    print data_dict['error']

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

def create_data_dicts(filepath):
    header_lines = 3
    with codecs.open(filepath) as f:
        raw_data = np.loadtxt(f, skiprows=header_lines)

    xdict = {}
    xdict['raw_data'] = raw_data[:,1:5]
    ydict = {}
    ydict['raw_data'] = raw_data[:,5:9]
    compute_data_dict_avg_min_max(xdict)
    compute_data_dict_avg_min_max(ydict)
    return xdict, ydict

def plot_tf(filepath, pk_npk_data_file = 'TRANSFER1_PK_NPK.txt'):
    filebase = filepath.split('.txt')[0]
    save_fig = filebase + '.png'

    xdict, ydict = create_data_dicts(filepath)
    xdata = xdict['average']
    ydata = ydict['average']

    pk_npk_dict = create_pk_npk_dict(pk_npk_data_file)

    ax.plot(xdata,ydata,'--ro',label='Meas. Avg.')
    _ = make_error_boxes(ax, xdata, ydata, xdict['error'], ydict['error'], facecolor='g', edgecolor='None', alpha=0.5)

    ax.plot(pk_npk_dict['pk-scyl'],pk_npk_dict['pk-spole'],'-bo',label='FEM PK')
    ax.plot(pk_npk_dict['npk-scyl'],pk_npk_dict['npk-spole'],'-bd',label='FEM NPK')
    ax.set_xlabel('Shell Azimuthal Stress (MPa)')
    ax.set_ylabel('Pole Azimuthal Stress (MPa)')
    ax.grid()

    ax.legend()
    plt.savefig(save_fig)

    #lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5,-.1), fancybox=True, shadow=True, ncol=2)
    #plt.savefig(save_fig, bbox_extra_artists=(lgd,), bbox_inches='tight')


if __name__ == '__main__':
    filepath = 'TRANSFER1_MQXFS6.txt'
    parser = argparse.ArgumentParser(description='Plot transfer function')
    parser.add_argument('paths', nargs='+', type=str)
    parser.add_argument('-t', '--type', type=int, default=1) 
    parser.add_argument('-pk', '--pk-npk-file', type=str, default='TRANSFER1_PK_NPK_simple.txt') 

    args = parser.parse_args()
    paths = args.paths
    tr_type = args.type
    pk_npk_file = args.pk_npk_file

    for filepath in paths:
        plt.cla()
        plot_tf(filepath, pk_npk_file)
    


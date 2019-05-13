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
#ax.autoscale(enable=True, axis='y', tight=True)

colors = ['b','g','r','c','m','y','k']
markers = ['o', '^', 'v', '<', '>', '1', '2', '3', '4', 's', 'p', '*', 'h', '+', 'x']

def read_coil_size_data(filepath):
    coil_size_dict = {}
    header_lines = 10
    with codecs.open(filepath) as f:
        head = [next(f).decode('utf-8', 'ignore') for x in xrange(header_lines)]

    for line in head:
        if 'shim' in line: shim = float(line.split()[1])*args.unit_scaling
    print "_________________________"
    print "filepath", filepath
    print "size shim:", shim

    col_names = ['L+R' if 'L+R' in s else s for s in head[-1].strip('\r\n').split('\t')]
    print "column names:", col_names
    coil_size_dict['shim'] = shim
    coil_size_dict['column_names'] = col_names
    coil_size_dict['filepath'] = filepath
    coil_size_dict['column_indices'] = {}
    col_inds = coil_size_dict['column_indices']
    for i,col_name in enumerate(col_names):
        col_inds[col_name] = i

    with codecs.open(filepath) as f:
        raw_data = np.loadtxt(f, skiprows=header_lines)
    raw_data *= args.unit_scaling

    coil_size_dict['raw_data'] = raw_data
    print "Add column L+R+shim"
    add_to_raw_data('L+R+shim', coil_size_dict, raw_data[:,col_inds['L+R']] + shim)

    print "Add centered position"
    pos = raw_data[:,col_inds['Y']]
    pos_avg = np.average(pos)
    add_to_raw_data('Y centered', coil_size_dict, pos - pos_avg)

    print coil_size_dict['column_names']
    print coil_size_dict['column_indices']

    return coil_size_dict

def add_to_raw_data(name, coil_size_dict, data):
    coil_size_dict['column_names'].append(name)
    coil_size_dict['column_indices'][name] = len(coil_size_dict['column_indices'])
    coil_size_dict['raw_data'] = np.c_[coil_size_dict['raw_data'], data]

def add_interp_to_coil_dict(coil_size_dict, args):
    col_names = coil_size_dict['column_names']
    col_inds = coil_size_dict['column_indices']
    
    if not args.no_centering: xind = col_inds['Y']
    else: xind = col_inds['Y centered']

    coil_size_dict['interp'] = []
    x_vec = coil_size_dict['raw_data'][:,xind]

    for col_name in col_names:
        y_vec = coil_size_dict['raw_data'][:,col_inds[col_name]]
        coil_size_dict['interp'].append(scipy.interpolate.InterpolatedUnivariateSpline(x_vec, y_vec,k=1))

def read_coil_size_dicts(args):
    return [read_coil_size_data(filepath) for filepath in args.paths]
     
def plot_coil_sizes(coil_size_dicts, args):
    for i,coil_size_dict in enumerate(coil_size_dicts):
        add_interp_to_coil_dict(coil_size_dict)
        xdata = coil_size_dict['raw_data'][:,0]
        ydata = coil_size_dict['raw_data'][:,1]
        label = coil_size_dict['filepath'] + coil_size_dict['column_names'][1]
        ax.plot(xdata,ydata,'--'+colors[i]+markers[i],label=label)
    plt.show()

def plot_interpolated_coil_sizes(coil_size_dicts, args, av_interp, shimmed_av_interp):
    for i,coil_size_dict in enumerate(coil_size_dicts):
        col_inds = coil_size_dict['column_indices']
        add_interp_to_coil_dict(coil_size_dict, args)
        if not args.no_centering: xind = col_inds['Y']
        else: xind = col_inds['Y centered']
        if i==0: xdata = coil_size_dict['raw_data'][:,xind]
        ydata = coil_size_dict['interp'][col_inds['L+R']](xdata)
        label = coil_size_dict['filepath'].replace('.size','') + coil_size_dict['column_names'][1]
        ax.plot(xdata,ydata,'--'+colors[i]+markers[i],label=label)
        totnum = i + 1
    if av_interp != None:
        ax.plot(xdata,av_interp(xdata),'-b'+markers[totnum],label='av', linewidth=5)
    if shimmed_av_interp != None:
        ax.plot(xdata,shimmed_av_interp(xdata),'-r'+markers[totnum],label='av+shim', linewidth=5)
    ax.legend(loc=args.legend_location)
    if args.show_plot:
        plt.show()
    else:
        save_fig = 'coil_size.png'
        plt.savefig(save_fig)

def get_average_coil_size_interp(coil_size_dicts, args, col_name='L+R'):
    for i,coil_size_dict in enumerate(coil_size_dicts):
        col_inds = coil_size_dict['column_indices']
        add_interp_to_coil_dict(coil_size_dict, args)
        if not args.no_centering: xind = col_inds['Y']
        else: xind = col_inds['Y centered']
        if i==0: 
            xdata = coil_size_dict['raw_data'][:,xind]
            ydata = coil_size_dict['interp'][col_inds[col_name]](xdata)
            LplusRav = coil_size_dict['interp'][col_inds[col_name]](xdata)
        else:
            ydata = coil_size_dict['interp'][col_inds[col_name]](xdata)
            LplusRav += coil_size_dict['interp'][col_inds[col_name]](xdata)
        totnum = i + 1
    LplusRav /= totnum 
    interp = scipy.interpolate.InterpolatedUnivariateSpline(xdata, LplusRav,k=1)
    return interp, xdata

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot coil sizes')
    parser.add_argument('paths', nargs='+', type=str)
    parser.add_argument('-sp', '--show-plot', action='store_true', default=False) 
    parser.add_argument('-ll', '--legend-location', type=str, default='best')
    parser.add_argument('-us', '--unit-scaling', type=float, default=0.001)
    parser.add_argument('-nc', '--no-centering', action='store_false', default=True) 

    args = parser.parse_args()
    print "read coil sizes"
    coil_size_dicts = read_coil_size_dicts(args)
    #plot_coil_sizes(coil_size_dicts, args)
    #plot_interpolated_coil_sizes(coil_size_dicts, args)
    average_coil_size_interp, av_xdata = get_average_coil_size_interp(coil_size_dicts,args,col_name='L+R',)
    average_coil_size_shim_interp, av_xdata = get_average_coil_size_interp(coil_size_dicts,args,col_name='L+R+shim')
    print "plot coil sizes"
    plot_interpolated_coil_sizes(coil_size_dicts, args, av_interp = average_coil_size_interp, shimmed_av_interp=average_coil_size_shim_interp)
    print "Store average coil size to average_coil.size"
    avdata = np.c_[av_xdata, average_coil_size_interp(av_xdata)]
    np.savetxt('average_coil.size', avdata, header='Y\n Average coil size')
    print "Store shimmed coil size to average_shimmed_coil.size"
    avshimdata = np.c_[av_xdata, average_coil_size_shim_interp(av_xdata)]
    np.savetxt('average_shimmed_coil.size', avshimdata, header='Y\n Average shimmed coil size')


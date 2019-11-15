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

def read_coil_size_data(filepath,args):
    coil_size_dict = {}
    header_lines = 10
    with codecs.open(filepath) as f:
        head = [next(f).decode('utf-8', 'ignore') for x in xrange(header_lines)]

    offset = None
    for line in head:
        if 'shim' in line: shim = float(line.split()[1])*args.unit_scaling
        elif 'offset' in line: offset = float(line.split()[1])*args.unit_scaling

    print "_________________________"
    print "filepath", filepath
    print "size shim:", shim

    col_names = ['L+R' if 'L+R' in s else s for s in head[-1].strip('\r\n').split('\t')]
    print "column names:", col_names
    coil_size_dict['shim'] = shim
    coil_size_dict['offset'] = offset
    coil_size_dict['column_names'] = col_names
    coil_size_dict['filepath'] = filepath
    coil_size_dict['radial_size_reduction'] = args.radial_size_reduction
    coil_size_dict['column_indices'] = {}
    col_inds = coil_size_dict['column_indices']
    for i,col_name in enumerate(col_names):
        col_inds[col_name] = i

    with codecs.open(filepath) as f:
        raw_data = np.loadtxt(f, skiprows=header_lines)
    raw_data *= args.unit_scaling

    coil_size_dict['raw_data'] = raw_data
    if coil_size_dict['offset'] != None:
        print "Using offset", coil_size_dict['offset'], "for correcting the L+R data!"
        coil_size_dict['raw_data'][:, col_inds['L+R']] += coil_size_dict['offset']

    print "Add column L+R+shim"
    add_to_raw_data('L+R+shim', coil_size_dict, raw_data[:,col_inds['L+R']] + shim)

    print "Add column L+R+shim-rsr"
    add_to_raw_data('L+R+shim-rsr', coil_size_dict, raw_data[:,col_inds['L+R']] + shim - args.radial_size_reduction)

    print "Add centered position"
    pos = raw_data[:,col_inds['Y']]
    pos_avg = np.average(pos)
    add_to_raw_data('Y centered', coil_size_dict, pos - pos_avg)

    coil_size_dict['L+R average'] = 1e6*np.mean(raw_data, axis=0)[col_inds['L+R']]
    print "L+R average: {:2.0f} um".format(coil_size_dict['L+R average'] )

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
    return [read_coil_size_data(filepath,args) for filepath in args.paths]
     
def plot_coil_sizes(coil_size_dicts, args):
    for i,coil_size_dict in enumerate(coil_size_dicts):
        add_interp_to_coil_dict(coil_size_dict)
        xdata = coil_size_dict['raw_data'][:,0] * args.xunit_plot_scaling
        ydata = coil_size_dict['raw_data'][:,1] * args.yunit_plot_scaling
        label = coil_size_dict['filepath'] + coil_size_dict['column_names'][1]
        ax.plot(xdata,ydata,'--'+colors[i]+markers[i],label=label)
    plt.show()

def plot_interpolated_coil_sizes(coil_size_dicts, args, av_interp, shimmed_av_interp=None, shimmed_rsr_av_interp=None):
    plt.cla()

    if args.set_xlim != None: ax.set_xlim(args.set_xlim)
    if args.set_ylim != None: ax.set_ylim(args.set_ylim)
    ax.grid(args.grid)

    for i,coil_size_dict in enumerate(coil_size_dicts):
        col_inds = coil_size_dict['column_indices']
        add_interp_to_coil_dict(coil_size_dict, args)
        if not args.no_centering: xind = col_inds['Y']
        else: xind = col_inds['Y centered']
        if i==0: xdata = coil_size_dict['raw_data'][:,xind]
        ydata = coil_size_dict['interp'][col_inds['L+R']](xdata)
        label = coil_size_dict['filepath'].replace('.size','')# + coil_size_dict['column_names'][1]jbbb
        label += " ({:2.0f} $\mu$m)".format(coil_size_dict['L+R average'])
        ax.plot(xdata* args.xunit_plot_scaling,ydata* args.yunit_plot_scaling,'--'+colors[i]+markers[i],label=label)
        totnum = i + 1
    ax.set_xlabel("Longitudinal location (m)")
    ax.set_ylabel("L+R ($\mu$m)")
    if av_interp != None:
        ax.plot(xdata* args.xunit_plot_scaling,av_interp(xdata)* args.yunit_plot_scaling,color='black', marker=markers[totnum],label='av', linewidth=3)
    #if shimmed_av_interp != None:
        #ax.plot(xdata* args.xunit_plot_scaling,shimmed_av_interp(xdata)* args.yunit_plot_scaling,'-r'+markers[totnum],label='av+shim', linewidth=5)
    #if shimmed_rsr_av_interp != None and args.radial_size_reduction != 0:
        #ax.plot(xdata* args.xunit_plot_scaling,shimmed_rsr_av_interp(xdata)* args.yunit_plot_scaling,'-y'+markers[totnum],label='av+shim-rsr', linewidth=5)
    ax.legend(loc=args.legend_location)
    if args.show_plot:
        plt.show()
    else:
        save_fig = 'coil_size.png'
        plt.savefig(save_fig)

def plot_interpolated_shimmed_coil_sizes(coil_size_dicts, args, av_interp, shimmed_av_interp=None, shimmed_rsr_av_interp=None):
    plt.cla()

    if args.set_xlim != None: ax.set_xlim(args.set_xlim)
    if args.set_ylim != None: ax.set_ylim(args.set_ylim)
    ax.grid(args.grid)

    for i,coil_size_dict in enumerate(coil_size_dicts):
        col_inds = coil_size_dict['column_indices']
        shim = coil_size_dict['shim']
        add_interp_to_coil_dict(coil_size_dict, args)
        if not args.no_centering: xind = col_inds['Y']
        else: xind = col_inds['Y centered']
        if i==0: xdata = coil_size_dict['raw_data'][:,xind]
        ydata = coil_size_dict['interp'][col_inds['L+R']](xdata) + shim
        label = coil_size_dict['filepath'].replace('.size','')# + coil_size_dict['column_names'][1]
        if not np.isclose(shim,0):
            label += '+shim ({:2.0f} $\mu$m)'.format(1e6*shim)
        ax.plot(xdata* args.xunit_plot_scaling,ydata* args.yunit_plot_scaling,'--'+colors[i]+markers[i],label=label)
        totnum = i + 1
    ax.set_xlabel("Longitudinal location (m)")
    ax.set_ylabel("L+R ($\mu$m)")
    if shimmed_av_interp != None:
        ax.plot(xdata* args.xunit_plot_scaling,shimmed_av_interp(xdata)* args.yunit_plot_scaling,color='black', marker=markers[totnum],label='av+shim', linewidth=3)
    ax.legend(loc=args.legend_location)
    if args.show_plot:
        plt.show()
    else:
        save_fig = 'coil_size_shimmed.png'
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
    parser.add_argument('-xus', '--xunit-plot-scaling', type=float, default=1)
    parser.add_argument('-yus', '--yunit-plot-scaling', type=float, default=1e6)
    parser.add_argument('-nc', '--no-centering', action='store_false', default=True) 
    parser.add_argument('-rsr', '--radial-size-reduction', type=float, default=0.)
    parser.add_argument('--fit', action='store_true', default=False)
    parser.add_argument('--set-xlim', nargs=2, type=float)
    parser.add_argument('--set-ylim', nargs=2, type=float)
    parser.add_argument('--grid', action='store_true', default=False)

    args = parser.parse_args()

    print "read coil sizes"
    coil_size_dicts = read_coil_size_dicts(args)
    #plot_coil_sizes(coil_size_dicts, args)
    average_coil_size_interp, av_xdata = get_average_coil_size_interp(coil_size_dicts,args,col_name='L+R',)
    average_coil_size_shim_interp, av_xdata = get_average_coil_size_interp(coil_size_dicts,args,col_name='L+R+shim')
    average_coil_size_shim_rsr_interp, av_xdata = get_average_coil_size_interp(coil_size_dicts,args,col_name='L+R+shim-rsr')
    print "plot coil sizes"
    plot_interpolated_coil_sizes(coil_size_dicts, args, av_interp = average_coil_size_interp, shimmed_av_interp=average_coil_size_shim_interp, shimmed_rsr_av_interp=average_coil_size_shim_rsr_interp)
    plot_interpolated_shimmed_coil_sizes(coil_size_dicts, args, av_interp = average_coil_size_interp, shimmed_av_interp=average_coil_size_shim_interp, shimmed_rsr_av_interp=average_coil_size_shim_rsr_interp)
    print "Store average coil size to average_coil.size"
    avdata = np.c_[av_xdata, average_coil_size_interp(av_xdata)]
    np.savetxt('average_coil.size', avdata, header='Y\n Average coil size')
    print "Store shimmed coil size to average_shimmed_coil.size"
    avshimdata = np.c_[av_xdata, average_coil_size_shim_interp(av_xdata)]
    np.savetxt('average_shimmed_coil.size', avshimdata, header='Y\n Average shimmed coil size')

    if args.radial_size_reduction != 0:
        print "Store shimmed coil size with radial size reduction to average_shimmed_coil_rsr.size"
        np.savetxt('average_shimmed_coil_rsr.size', avshimrsrdata, header='Y\n Average shimmed rsr coil size')
        avshimrsrdata = np.c_[av_xdata, average_coil_size_shim_rsr_interp(av_xdata)]


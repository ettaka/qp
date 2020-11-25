#!env python

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
#import pandas as pd
import argparse
import codecs
import scipy.interpolate
from scipy.stats import linregress

from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import odr_fit


fig = plt.figure()
ax = fig.add_subplot(111)
#ax.autoscale(enable=True, axis='y', tight=True)

colors = ['b','g','r','c','m','y','k']
markers = ['o', '^', 'v', '<', '>', '1', '2', '3', '4', 's', 'p', '*', 'h', '+', 'x']

def read_coil_size_data(filepath,args):
    coil_size_dict = {}
    header_lines = 10
    with open(filepath) as f:
        head = [next(f) for x in range(header_lines)]

    offset = None
    for line in head:
        if 'mshim' in line: mshim = float(line.split()[1])*args.unit_scaling
        elif 'offset' in line: offset = float(line.split()[1])*args.unit_scaling
        elif 'rshim' in line: rshim = float(line.split()[1])*args.unit_scaling

    if args.rshim != None: rshim = args.rshim

    print("_________________________")
    print("filepath", filepath)
    print("mid-plane shim:", mshim)
    print("radial shim:", rshim)
    print("offset", offset)

    col_names = ['L+R' if ('L+R' in s and 'spread' not in s) else s for s in head[-1].strip('\r\n').split('\t')]
    print("column names:", col_names)
    coil_size_dict['mshim'] = mshim
    coil_size_dict['rshim'] = rshim
    coil_size_dict['offset'] = offset
    coil_size_dict['column_names'] = col_names
    coil_size_dict['filepath'] = filepath
    coil_size_dict['radial_size_reduction'] = args.radial_size_reduction
    coil_size_dict['column_indices'] = {}
    col_inds = coil_size_dict['column_indices']
    for i,col_name in enumerate(col_names):
        print (i, col_name)
        col_inds[col_name] = i

    with open(filepath) as f:
        raw_data = np.loadtxt(f, skiprows=header_lines)
    raw_data *= args.unit_scaling

    coil_size_dict['raw_data'] = raw_data
    if coil_size_dict['offset'] != None:
        print("Using offset", coil_size_dict['offset'], "for correcting the L+R data!")
        coil_size_dict['raw_data'][:, col_inds['L+R']] += coil_size_dict['offset']

    print("Add column L+R+mshim")
    add_to_raw_data('L+R+mshim', coil_size_dict, raw_data[:,col_inds['L+R']] + mshim)

    print("Add column L+R+mshim-rsr")
    add_to_raw_data('L+R+mshim-rsr', coil_size_dict, raw_data[:,col_inds['L+R']] + mshim - args.radial_size_reduction)

    print("Add centered position")
    pos = raw_data[:,col_inds['Y']]
    pos_avg = np.average(pos)
    add_to_raw_data('Y centered', coil_size_dict, pos - pos_avg)

    coil_size_dict['L+R average'] = 1e6*np.mean(raw_data, axis=0)[col_inds['L+R']]
    print("L+R average: {:2.0f} um".format(coil_size_dict['L+R average'] ))

    print(coil_size_dict['column_names'])
    print(coil_size_dict['column_indices'])

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
    plot_station_vertical_lines()
    plt.show()

def plot_station_vertical_lines():
    if args.short_magnet:
        CE = 1.55/2.
        ax.axvline(x=CE, color='black', linestyle='dashed')
    else:
        LE = .607
        CE = 3.407
        RE = 7.007
        ax.axvline(x=LE, color='black', linestyle='dashed')
        ax.axvline(x=CE, color='black', linestyle='dashed')
        ax.axvline(x=RE, color='black', linestyle='dashed')

def plot_interpolated_coil_sizes(coil_size_dicts, args, av_interp, mshimmed_av_interp=None, shimmed_rsr_av_interp=None):
    plt.cla()
    totnum = 0

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
    plot_station_vertical_lines()

    if args.plot_note is not None:
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)

        # place a text box in upper left in axes coords
        ax.text(0.135, 0.95, args.plot_note, transform=ax.transAxes, fontsize=12,
        horizontalalignment='left', verticalalignment='top', bbox=props)

    ax.set_xlabel("Longitudinal location (m)",fontsize=args.font_size)
    ax.set_ylabel("L+R ($\mu$m)", fontsize=args.font_size)
    if av_interp != None:
        ax.plot(xdata* args.xunit_plot_scaling,av_interp(xdata)* args.yunit_plot_scaling,color='black', marker=markers[totnum],label='av', linewidth=3)
    #if mshimmed_av_interp != None:
        #ax.plot(xdata* args.xunit_plot_scaling,mshimmed_av_interp(xdata)* args.yunit_plot_scaling,'-r'+markers[totnum],label='av+mshim', linewidth=5)
    #if shimmed_rsr_av_interp != None and args.radial_size_reduction != 0:
        #ax.plot(xdata* args.xunit_plot_scaling,shimmed_rsr_av_interp(xdata)* args.yunit_plot_scaling,'-y'+markers[totnum],label='av+mshim-rsr', linewidth=5)
    ax.legend(loc=args.legend_location)
    if args.show_plot:
        plt.show()
    else:
        save_fig = 'coil_size.png'
        plt.savefig(save_fig, bbox_inches='tight', numpoints=1, dpi=200)

def plot_interpolated_shimmed_coil_sizes(coil_size_dicts, args, av_interp, mshimmed_av_interp=None, shimmed_rsr_av_interp=None):
    plt.cla()
    totnum = 0

    if args.set_xlim != None: ax.set_xlim(args.set_xlim)
    if args.set_ylim != None: ax.set_ylim(args.set_ylim)
    ax.grid(args.grid)

    for i,coil_size_dict in enumerate(coil_size_dicts):
        col_inds = coil_size_dict['column_indices']
        mshim = coil_size_dict['mshim']
        add_interp_to_coil_dict(coil_size_dict, args)
        if not args.no_centering: xind = col_inds['Y']
        else: xind = col_inds['Y centered']
        if i==0: xdata = coil_size_dict['raw_data'][:,xind]
        ydata = coil_size_dict['interp'][col_inds['L+R']](xdata) + mshim
        label = coil_size_dict['filepath'].replace('.size','')# + coil_size_dict['column_names'][1]
        label += '+mshim ({:2.0f} $\mu$m)'.format(coil_size_dict['L+R average']+1e6*mshim)
        ax.plot(xdata* args.xunit_plot_scaling,ydata* args.yunit_plot_scaling,'--'+colors[i]+markers[i],label=label)
        totnum = i + 1
    plot_station_vertical_lines()
    ax.set_xlabel("Longitudinal location (m)")
    ax.set_ylabel("L+R ($\mu$m)")
    if mshimmed_av_interp != None:
        ax.plot(xdata* args.xunit_plot_scaling,mshimmed_av_interp(xdata)* args.yunit_plot_scaling,color='black', marker=markers[totnum],label='av+mshim', linewidth=3)
    ax.legend(loc=args.legend_location)
    if args.show_plot:
        plt.show()
    else:
        save_fig = 'coil_size_shimmed.png'
        plt.savefig(save_fig, bbox_inches='tight', numpoints=1, dpi=200)

def plot_interpolated_theoretical_load_key(coil_size_dicts, args, av_interp, mshimmed_av_interp=None, shimmed_rsr_av_interp=None, shell_slope=0.12, pole_slope=-0.20):
    plt.cla()
    totnum = 0

    if args.set_theoretical_key_ylim != None: ax.set_ylim(args.set_theoretical_key_ylim)
    ax.grid(args.grid)

    ax2 = ax.twinx()
    ax3 = ax.twinx()
    #if args.set_xlim != None: ax.set_xlim(args.set_xlim)

    for i,coil_size_dict in enumerate(coil_size_dicts):
        col_inds = coil_size_dict['column_indices']
        mshim = coil_size_dict['mshim']
        add_interp_to_coil_dict(coil_size_dict, args)
        if not args.no_centering: xind = col_inds['Y']
        else: xind = col_inds['Y centered']
        if i==0: xdata = coil_size_dict['raw_data'][:,xind]
        ydata = coil_size_dict['interp'][col_inds['L+R']](xdata) + mshim
        label = coil_size_dict['filepath'].replace('.size','')# + coil_size_dict['column_names'][1]
        label += '+mshim ({:2.0f} $\mu$m)'.format(coil_size_dict['L+R average']+1e6*mshim)
        totnum = i + 1
    ax.set_xlabel("Longitudinal location (m)")
    ax.set_ylabel("Size w.r.t. minimum ($\mu$m)")
    ax2.set_ylabel("Azimuthal shell stress w.r.t. min. (MPa)")
    ax3.set_ylabel("Azimuthal pole compression w.r.t. min. (MPa)")
    ax3.spines["right"].set_position(("axes", 1.2))
    fig.subplots_adjust(right=0.75)
    if mshimmed_av_interp != None:
        ax.plot(xdata* args.xunit_plot_scaling,2./np.pi*(np.max(mshimmed_av_interp(xdata)) - mshimmed_av_interp(xdata))* args.yunit_plot_scaling,color='black', marker=markers[totnum],label='Theoretical load key', linewidth=3)
        ax2.plot(xdata* args.xunit_plot_scaling, shell_slope * 2./np.pi*(np.max(mshimmed_av_interp(xdata)) - mshimmed_av_interp(xdata))* args.yunit_plot_scaling,color='black', marker=markers[totnum],label='Theoretical load key', linewidth=3)
        mn, mx = ax.get_ylim()
        ax2.set_ylim(shell_slope * mn, shell_slope * mx)
        ax3.set_ylim(-pole_slope * mn, -pole_slope * mx)
    plot_station_vertical_lines()
    ax.legend(loc=args.legend_location)
    if args.show_plot:
        plt.show()
    else:
        save_fig = 'theoretical_load_key.png'
        plt.savefig(save_fig, bbox_inches='tight', numpoints=1, dpi=200)

def plot_interpolated_coil_pack_z(coil_size_dicts, args, mshimmed_av_interp, gaps=None):
    plt.cla()
    totnum = 0
    if args.set_xlim != None: ax.set_xlim(args.set_gaps_xlim)
    if args.set_ylim != None: ax.set_ylim(args.set_gaps_ylim)
    ax.grid(args.grid)

    for i,coil_size_dict in enumerate(coil_size_dicts):
        col_inds = coil_size_dict['column_indices']
        mshim = coil_size_dict['mshim']
        rshim = coil_size_dict['rshim']
        add_interp_to_coil_dict(coil_size_dict, args)
        if not args.no_centering: xind = col_inds['Y']
        else: xind = col_inds['Y centered']
        if i==0: xdata = coil_size_dict['raw_data'][:,xind]
        ydata = coil_size_dict['interp'][col_inds['L+R']](xdata) + mshim + rshim*math.pi/2.
        label = coil_size_dict['filepath'].replace('.size','')# + coil_size_dict['column_names'][1]
        #label += '+mshim ({:2.0f} $\mu$m)'.format(coil_size_dict['L+R average']+1e6*mshim)
        #totnum = i + 1
    ax.set_xlabel("Longitudinal location (m)")
    ax.set_ylabel("Size w.r.t. nominal ($\mu$m)")
    ax.plot(xdata * args.xunit_plot_scaling, (mshimmed_av_interp(xdata)+ rshim*math.pi/2.) * args.yunit_plot_scaling,color='red', marker=markers[totnum],label='Shimmed coil average sector length', linewidth=3)
    if gaps is not None:
        for j, gap_data in enumerate(gaps):
            gap_file_name = gap_data['file name']
            data_name = gap_file_name.replace('.size', '').replace('gaps','').replace('_','')
            gap_df = gap_data['DataFrame']
            totnum = totnum + 1
            theoretical_gap = 15.
            coilpack_to_pad = 10.
            gap_data_average = gap_df['Gap AVG'] 
            xdata = (gap_df['Y']-coilpack_to_pad)/1000. * args.xunit_plot_scaling
            ydata = (gap_data_average-theoretical_gap)/1000. * args.yunit_plot_scaling
            ax.plot(xdata, ydata, '--'+colors[totnum]+markers[totnum],label='Average collar gap '+data_name, linewidth=3)

            if 'Boltnr' in gap_df:
                for i in range(len(xdata)):
                    x0 = xdata[i]
                    y0 = ydata[i]+50
                    ax.annotate(gap_df['Boltnr'][i], (x0, y0),  ha='center', va='center')#, fontsize=annotation_font_size)
            elif j == 0:
                for i in range(len(xdata)):
                    x0 = xdata[i]
                    y0 = ydata[i]+50
                    ax.annotate(str(i+1), (x0, y0),  ha='center', va='center')#, fontsize=annotation_font_size)
    plot_station_vertical_lines()
    ax.legend(loc=args.legend_location)
    if args.show_plot:
        plt.show()
    else:
        save_fig = 'coil_average_sector_length.png'
        plt.savefig(save_fig, bbox_inches='tight', numpoints=1, dpi=200)

def plot_interpolated_collar_gap_vs_coil_sec_len(coil_size_dicts, args, mshimmed_av_interp, gaps=None):
    plt.cla()
    totnum = 0
    #if args.set_xlim != None: ax.set_xlim(args.set_gaps_xlim)
    #if args.set_ylim != None: ax.set_ylim(args.set_gaps_ylim)
    if args.set_xlim != None: ax.set_xlim(auto=True)
    if args.set_ylim != None: ax.set_ylim(auto=True)
    ax.grid(args.grid)

    for i,coil_size_dict in enumerate(coil_size_dicts):
        col_inds = coil_size_dict['column_indices']
        mshim = coil_size_dict['mshim']
        rshim = coil_size_dict['rshim']
        add_interp_to_coil_dict(coil_size_dict, args)
        if not args.no_centering: xind = col_inds['Y']
        else: xind = col_inds['Y centered']
        if i==0: xdata = coil_size_dict['raw_data'][:,xind]
        ydata = coil_size_dict['interp'][col_inds['L+R']](xdata) + mshim + rshim*math.pi/2.
        label = coil_size_dict['filepath'].replace('.size','')# + coil_size_dict['column_names'][1]
        #label += '+mshim ({:2.0f} $\mu$m)'.format(coil_size_dict['L+R average']+1e6*mshim)
        #totnum = i + 1
    #ax.set_xlabel("Longitudinal location (m)")
    ax.set_xlabel("Coil sector length w.r.t. nominal ($\mu$m)")
    ax.set_ylabel("Collar gap w.r.t. nominal ($\mu$m)")
    #ax.plot(xdata * args.xunit_plot_scaling, (mshimmed_av_interp(xdata)+ rshim*math.pi/2.) * args.yunit_plot_scaling,color='red', marker=markers[totnum],label='Shimmed coil average sector length', linewidth=3)
    if gaps is not None:
        for j, gap_data in enumerate(gaps):
            gap_file_name = gap_data['file name']
            data_name = gap_file_name.replace('.size', '').replace('gaps','').replace('_','')
            gap_df = gap_data['DataFrame']
            totnum = totnum + 1
            theoretical_gap = 15.
            coilpack_to_pad = 10.
            gap_data_average = gap_df['Gap AVG'] 
            xdata_gap = (gap_df['Y']-coilpack_to_pad)/1000. * args.xunit_plot_scaling
            ydata_gap = (gap_data_average-theoretical_gap)/1000. * args.yunit_plot_scaling 
            ydata_coil_len = (mshimmed_av_interp(xdata_gap) + rshim*math.pi/2.)*1e6

            coef,cov = np.polyfit(ydata_coil_len, ydata_gap, deg=1,cov=True)
            slope, intercept, r_value, p_value, std_err = linregress(ydata_coil_len, ydata_gap)
            fit = np.poly1d(coef)
            sigma = np.sqrt(cov)

            fit_label='fit: ({:0.2f}$\pm${:0.1f})x + ({:2.0f}$\pm${:.0f})\n$R^2$={:.2f}'.format(coef[0],sigma[0,0], coef[1],sigma[1,1], r_value**2.)
            ax.plot(ydata_coil_len, ydata_gap, 'o',color=colors[totnum],label='measured', linewidth=3)
            ax.plot(ydata_coil_len, fit(ydata_coil_len), '-',color='black',label=fit_label, linewidth=3)

            #if 'Boltnr' in gap_df:
            #    for i in range(len(xdata)):
            #        x0 = xdata[i]
            #        y0 = ydata[i]+50
            #        ax.annotate(gap_df['Boltnr'][i], (x0, y0),  ha='center', va='center')#, fontsize=annotation_font_size)
            #elif j == 0:
            #    for i in range(len(xdata)):
            #        x0 = xdata[i]
            #        y0 = ydata[i]+50
            #        ax.annotate(str(i+1), (x0, y0),  ha='center', va='center')#, fontsize=annotation_font_size)
    #plot_station_vertical_lines()
    ax.legend(loc=args.legend_location)
    if args.show_plot:
        plt.show()
    else:
        save_fig = 'collar_gap_vs_coil_len.png'
        plt.savefig(save_fig, bbox_inches='tight', numpoints=1, dpi=200)

def plot_interpolated_collar_radius_vs_coil_pack_size(coil_size_dicts, args, mshimmed_av_interp, gaps=None):
    plt.cla()
    totnum = 0
    #if args.set_xlim != None: ax.set_xlim(args.set_gaps_xlim)
    #if args.set_ylim != None: ax.set_ylim(args.set_gaps_ylim)
    if args.set_xlim != None: ax.set_xlim(auto=True)
    if args.set_ylim != None: ax.set_ylim(auto=True)
    ax.grid(args.grid)

    for i,coil_size_dict in enumerate(coil_size_dicts):
        col_inds = coil_size_dict['column_indices']
        mshim = coil_size_dict['mshim']
        rshim = coil_size_dict['rshim']
        add_interp_to_coil_dict(coil_size_dict, args)
        if not args.no_centering: xind = col_inds['Y']
        else: xind = col_inds['Y centered']
        if i==0: xdata = coil_size_dict['raw_data'][:,xind]
        ydata = coil_size_dict['interp'][col_inds['L+R']](xdata) + mshim + rshim*math.pi/2.
        label = coil_size_dict['filepath'].replace('.size','')# + coil_size_dict['column_names'][1]
        #label += '+mshim ({:2.0f} $\mu$m)'.format(coil_size_dict['L+R average']+1e6*mshim)
        #totnum = i + 1
    #ax.set_xlabel("Longitudinal location (m)")
    ax.set_xlabel("Coil pack size w.r.t. nominal ($\mu$m)")
    ax.set_ylabel("Collar radius w.r.t. nominal ($\mu$m)")
    #ax.plot(xdata * args.xunit_plot_scaling, (mshimmed_av_interp(xdata)+ rshim*math.pi/2.) * args.yunit_plot_scaling,color='red', marker=markers[totnum],label='Shimmed coil average sector length', linewidth=3)
    if gaps is not None:
        for j, gap_data in enumerate(gaps):
            gap_file_name = gap_data['file name']
            data_name = gap_file_name.replace('.size', '').replace('gaps','').replace('_','')
            gap_df = gap_data['DataFrame']
            totnum = totnum + 1
            theoretical_gap = 15.
            coilpack_to_pad = 10.
            gap_data_average = gap_df['Gap AVG'] 
            xdata_gap = (gap_df['Y']-coilpack_to_pad)/1000. * args.xunit_plot_scaling
            ydata_gap = (gap_data_average-theoretical_gap)/1000. * args.yunit_plot_scaling 
            ydata_coil_len = (mshimmed_av_interp(xdata_gap) + rshim*math.pi/2.)*1e6

            ydata_coil_pack_size = 2. * ydata_coil_len/math.pi
            collar_radius = 2. * ydata_gap/math.pi

            coef, coef_sd, reduced_chi_square = odr_fit.fit(ydata_coil_pack_size,collar_radius)
            fit = np.poly1d(coef)

            fit_label='fit: ({:0.2f}$\pm${:0.1f})x + ({:2.0f}$\pm${:.0f})\n$\chi_\\nu^2$={:.2f}'.format(coef[0],coef_sd[0], coef[1],coef_sd[1], reduced_chi_square)
            ax.plot(ydata_coil_pack_size, collar_radius, 'o',color=colors[totnum],label='measured', linewidth=3)
            ax.plot(ydata_coil_pack_size, fit(ydata_coil_pack_size), '-',color='black',label=fit_label, linewidth=3)

            #if 'Boltnr' in gap_df:
            #    for i in range(len(xdata)):
            #        x0 = xdata[i]
            #        y0 = ydata[i]+50
            #        ax.annotate(gap_df['Boltnr'][i], (x0, y0),  ha='center', va='center')#, fontsize=annotation_font_size)
            #elif j == 0:
            #    for i in range(len(xdata)):
            #        x0 = xdata[i]
            #        y0 = ydata[i]+50
            #        ax.annotate(str(i+1), (x0, y0),  ha='center', va='center')#, fontsize=annotation_font_size)
    #plot_station_vertical_lines()
    ax.legend(loc=args.legend_location)
    if args.show_plot:
        plt.show()
    else:
        save_fig = 'collar_radius_vs_coil_pack_size.png'
        plt.savefig(save_fig, bbox_inches='tight', numpoints=1, dpi=200)

def plot_interpolated_collar_gaps(coil_size_dicts, args, mshimmed_av_interp, gaps=None):
    plt.cla()
    totnum = 0
    #if args.set_xlim != None: ax.set_xlim(args.set_gaps_xlim)
    #if args.set_ylim != None: ax.set_ylim(args.set_gaps_ylim)
    if args.set_xlim != None: ax.set_xlim(auto=True)
    if args.set_ylim != None: ax.set_ylim(auto=True)
    ax.grid(args.grid)

    for i,coil_size_dict in enumerate(coil_size_dicts):
        col_inds = coil_size_dict['column_indices']
        mshim = coil_size_dict['mshim']
        rshim = coil_size_dict['rshim']
        add_interp_to_coil_dict(coil_size_dict, args)
        if not args.no_centering: xind = col_inds['Y']
        else: xind = col_inds['Y centered']
        if i==0: xdata = coil_size_dict['raw_data'][:,xind]
        ydata = coil_size_dict['interp'][col_inds['L+R']](xdata) + mshim + rshim*math.pi/2.
        label = coil_size_dict['filepath'].replace('.size','')# + coil_size_dict['column_names'][1]
        #label += '+mshim ({:2.0f} $\mu$m)'.format(coil_size_dict['L+R average']+1e6*mshim)
        #totnum = i + 1
    ax.set_xlabel("Longitudinal location (m)")
    ax.set_ylabel("Gap size (mm)")
    #ax.plot(xdata * args.xunit_plot_scaling, (mshimmed_av_interp(xdata)+ rshim*math.pi/2.) * args.yunit_plot_scaling,color='red', marker=markers[totnum],label='Shimmed coil average sector length', linewidth=3)
    if gaps is not None:
        for j, gap_data in enumerate(gaps):
            gap_file_name = gap_data['file name']
            data_name = gap_file_name.replace('.size', '').replace('gaps','').replace('_','')
            gap_df = gap_data['DataFrame']
            totnum = totnum + 1
            theoretical_gap = 15.
            coilpack_to_pad = 10.
            gap_data_a = gap_df['Gap A'] 
            gap_data_b = gap_df['Gap B'] 
            gap_data_c = gap_df['Gap C'] 
            gap_data_d = gap_df['Gap D'] 
            gap_data_avg = gap_df['Gap AVG'] 
            xdata_gap = (gap_df['Y']-coilpack_to_pad)/1000. * args.xunit_plot_scaling

            ax.plot(xdata_gap, gap_data_a, '-'+markers[0],color=colors[0],label='A', linewidth=1)
            ax.plot(xdata_gap, gap_data_b, '-'+markers[1],color=colors[1],label='B', linewidth=1)
            ax.plot(xdata_gap, gap_data_c, '-'+markers[2],color=colors[2],label='C', linewidth=1)
            ax.plot(xdata_gap, gap_data_d, '-'+markers[3],color=colors[3],label='D', linewidth=1)
            ax.plot(xdata_gap, gap_data_avg, '-',color='black',label='Avg.', linewidth=3)

            #if 'Boltnr' in gap_df:
            #    for i in range(len(xdata)):
            #        x0 = xdata[i]
            #        y0 = ydata[i]+50
            #        ax.annotate(gap_df['Boltnr'][i], (x0, y0),  ha='center', va='center')#, fontsize=annotation_font_size)
            #elif j == 0:
            #    for i in range(len(xdata)):
            #        x0 = xdata[i]
            #        y0 = ydata[i]+50
            #        ax.annotate(str(i+1), (x0, y0),  ha='center', va='center')#, fontsize=annotation_font_size)
    #plot_station_vertical_lines()
    ax.legend(loc=args.legend_location)
    if args.show_plot:
        plt.show()
    else:
        save_fig = 'collar_gaps.png'
        plt.savefig(save_fig, bbox_inches='tight', numpoints=1, dpi=200)


def plot_interpolated_pole_key_gaps(coil_size_dicts, args, mshimmed_av_interp, gaps=None):
    plt.cla()
    totnum = 0
    #if args.set_xlim != None: ax.set_xlim(args.set_gaps_xlim)
    #if args.set_ylim != None: ax.set_ylim(args.set_gaps_ylim)
    if args.set_xlim != None: ax.set_xlim(auto=True)
    if args.set_ylim != None: ax.set_ylim(auto=True)
    ax.set_ylim((0,700))
    ax.grid(args.grid)

    ax.set_xlabel("Longitudinal location (m)")
    ax.set_ylabel("Gap size (mm)")
    if gaps is not None:
        for j, gap_data in enumerate(gaps):
            gap_file_name = gap_data['file name']
            data_name = gap_file_name.replace('.size', '').replace('gaps','').replace('_','')
            gap_df = gap_data['DataFrame']
            totnum = totnum + 1
            theoretical_gap = 15.
            pole_key_size = 13.9
            coilpack_to_pad = 10.

            nominal_polekey = 13.9
            nominal_GI = 0.125
            nominal_polekey_gap = (theoretical_gap - nominal_polekey - 4*nominal_GI)/2.

            gap_data_a = gap_df['Gap A'] 
            gap_data_b = gap_df['Gap B'] 
            gap_data_c = gap_df['Gap C'] 
            gap_data_d = gap_df['Gap D'] 
            gap_data_avg = gap_df['Gap AVG'] 
            gap_data_all = gap_df[['Gap A', 'Gap B', 'Gap C', 'Gap D']]
            gap_data_max = gap_data_all.max(axis=1)
            gap_data_min = gap_data_all.min(axis=1)
            xdata_gap = (gap_df['Y']-coilpack_to_pad)/1000. * args.xunit_plot_scaling
            pole_key_gap_a = (gap_data_a-nominal_polekey-4*nominal_GI)/2.*1e3
            pole_key_gap_b = (gap_data_b-nominal_polekey-4*nominal_GI)/2.*1e3
            pole_key_gap_c = (gap_data_c-nominal_polekey-4*nominal_GI)/2.*1e3
            pole_key_gap_d = (gap_data_d-nominal_polekey-4*nominal_GI)/2.*1e3
            pole_key_gap_avg = (gap_data_avg-nominal_polekey-4*nominal_GI)/2.*1e3
            pole_key_gap_max = (gap_data_max-nominal_polekey-4*nominal_GI)/2.*1e3
            pole_key_gap_min = (gap_data_min-nominal_polekey-4*nominal_GI)/2.*1e3
            theoretical_pole_key_gap = nominal_polekey_gap*1e3

            #ax.plot(xdata_gap, pole_key_gap_a, '-'+markers[0],color=colors[0],label='A', linewidth=1)
            #ax.plot(xdata_gap, pole_key_gap_b, '-'+markers[1],color=colors[1],label='B', linewidth=1)
            #ax.plot(xdata_gap, pole_key_gap_c, '-'+markers[2],color=colors[2],label='C', linewidth=1)
            #ax.plot(xdata_gap, pole_key_gap_d, '-'+markers[3],color=colors[3],label='D', linewidth=1)

            ax.plot(xdata_gap, pole_key_gap_max, '-',color='black',linestyle = 'dashed', label='Max.', linewidth=3)
            ax.plot(xdata_gap, pole_key_gap_avg, '-',color='black',label='Avg.', linewidth=3)
            ax.plot(xdata_gap, pole_key_gap_min, '-',color='black',linestyle = 'dashed', label='Min.', linewidth=3)
            
            ax.plot(xdata_gap, np.zeros_like(xdata_gap)+theoretical_pole_key_gap, '-',color='green',linestyle = 'dashed', label='Nom.', linewidth=2)

            #if 'Boltnr' in gap_df:
            #    for i in range(len(xdata)):
            #        x0 = xdata[i]
            #        y0 = ydata[i]+50
            #        ax.annotate(gap_df['Boltnr'][i], (x0, y0),  ha='center', va='center')#, fontsize=annotation_font_size)
            #elif j == 0:
            #    for i in range(len(xdata)):
            #        x0 = xdata[i]
            #        y0 = ydata[i]+50
            #        ax.annotate(str(i+1), (x0, y0),  ha='center', va='center')#, fontsize=annotation_font_size)
    #plot_station_vertical_lines()
    ax.legend(loc=args.legend_location)
    if args.show_plot:
        plt.show()
    else:
        save_fig = 'polekey_gaps.png'
        plt.savefig(save_fig, bbox_inches='tight', numpoints=1, dpi=200)

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

def get_shimming_info(coil_size_dicts, info_locations=None, location_length=.05, location_length_points=100, coilpackdiff = -100., gaps=None):
    if info_locations == None:
        short_shell = 0.3415
        long_shell = 2. * short_shell
        #LE = short_shell + long_shell
        #CE = l / 2.
        #RE = l - LE
        if args.short_magnet:
            l = 1.510
            CE = l/2.
            info_locations = {'AVG':None, 'CE':CE}
        else:
            l = 2*short_shell + 10*long_shell
            LE = .607
            CE = 3.407
            RE = 7.007
            info_locations = {'AVG':None, 'LE':LE, 'CE':CE, 'RE':RE}

    shimming_info = {"Location name":[],
                    "Location":[],
                    "Loc len":[],
                    "Coil":[],
                    "L+R":[],
                    "mshim":[], 
                    "L+R+mshim":[],
                    "dR":[],
                    "rshim": [],
                    "Coilpack": [],
                    "Coil average sector length": [],
                    "Coilpack e": [],
                    "Pole key gap": []
                    }
    
    if gaps is not None:
        gaps_interp = gaps[-1]['Gap AVGinterp']

    nominal_collar_gap = 15
    nominal_polekey = 13.9
    nominal_GI = 0.125

    for loc_name in info_locations:
        for i,coil_size_dict in enumerate(coil_size_dicts):
            shimming_info['Location name'].append(loc_name)
            location = info_locations[loc_name]
            col_inds = coil_size_dict['column_indices']
            mshim = 1e6*coil_size_dict['mshim']
            rshim = 1e6*coil_size_dict['rshim']
            offset = 1e6*coil_size_dict['rshim']
            if loc_name == 'AVG':
                coil_size = coil_size_dict['L+R average']
                shimming_info['Location'].append(1e3*l/2.)
                shimming_info['Loc len'].append(1e3*l)
                location_start = coil_size_dict['raw_data'][0,0]
                location_end = coil_size_dict['raw_data'][-1,0]
                location_area = np.linspace(location_start, location_end, location_length_points)
                if gaps is not None:
                    collar_gap = np.average(gaps_interp(1e3*location_area))
                    pole_key_gap = (collar_gap - 4*nominal_GI - nominal_polekey)/2
            else:
                location_start = location-location_length/2.
                location_end = location+location_length/2.
                location_area = np.linspace(location_start, location_end, location_length_points)
                coil_size = 1e6*(np.average(coil_size_dict['interp'][col_inds['L+R']](location_area)))
                if gaps is not None:
                    collar_gap = np.average(gaps_interp(1e3*location_area))
                    pole_key_gap = (collar_gap - 4*nominal_GI - nominal_polekey)/2
                shimming_info['Location'].append(1e3*location)
                shimming_info['Loc len'].append(1e3*location_length)

            shimming_info['Coil'].append(coil_size_dict['filepath'].replace('.size',''))
            shimming_info['L+R'].append(coil_size)
            shimming_info['mshim'].append(mshim)
            shimming_info['rshim'].append(rshim)
            shimming_info['L+R+mshim'].append(coil_size+mshim)
            shimming_info['dR'].append(2./math.pi *(coil_size+mshim))
            shimming_info['Coilpack'].append(2./math.pi *(coil_size+mshim)+rshim)
            coil_pack = shimming_info['Coilpack'][-1]
            shimming_info['Coil average sector length'].append(coil_size+mshim+rshim*math.pi/2.)
            shimming_info['Coilpack e'].append(2./math.pi *(coil_size+mshim)+rshim+coilpackdiff)
            if gaps is not None:
                shimming_info['Pole key gap'].append(1e3*pole_key_gap)
            else:
                shimming_info['Pole key gap'].append(None)


    sort_values = ['Loc len', 'Location']
    df = pd.DataFrame(shimming_info)
    cols = ["Coil", "Location name", "Location", "Loc len", "L+R", "mshim", "L+R+mshim", "dR", "rshim", "Coilpack", "Coil average sector length", "Coilpack e", "Pole key gap"]

    df_out = df[cols].sort_values(sort_values)
    #print(df_out.to_string(float_format="{:0.0f}".format, index=False))
    df_out.to_csv('shimming_info.csv')

    avg_table = df.groupby(['Location name','Location'], as_index=False, group_keys=True).mean()
    avg_table['Coil'] = ['AVG' for dR in avg_table['dR']]
    #print(avg_table[cols].sort_values(sort_values).to_string(float_format="{:0.0f}".format, index=False))
    avg_table.to_csv('shimming_info_avg_rshim'+str(rshim)+'.csv')

    joined_table = df_out.append(avg_table)

    #print(joined_table[cols].sort_values(sort_values).to_string(float_format="{:0.0f}".format, index=False))
    joined_table[cols].sort_values(sort_values).to_csv('shimming_info_joined'+str(rshim)+'.csv')

    LpR = {}
    locations = ['AVG','LE','CE','RE']
    LpR['Coil'] = list(joined_table[joined_table['Location name'].str.contains('LE')]['Coil'])
    for loc in locations:
        LpR[loc] = list(joined_table[joined_table['Location name'].str.contains(loc)]['L+R'])

    
    df_LpR = pd.DataFrame(LpR)

    df_LpR.to_csv('LplusR_out'+str(rshim)+'.csv')


    return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot coil sizes')
    parser.add_argument('paths', nargs='+', type=str)
    parser.add_argument('-sp', '--show-plot', action='store_true', default=False) 
    parser.add_argument('-ll', '--legend-location', type=str, default='best')
    parser.add_argument('--plot-note', type=str, default=None)
    parser.add_argument('-us', '--unit-scaling', type=float, default=0.001)
    parser.add_argument('-xus', '--xunit-plot-scaling', type=float, default=1)
    parser.add_argument('-yus', '--yunit-plot-scaling', type=float, default=1e6)
    parser.add_argument('-rshim', type=float, default=None)
    parser.add_argument('-nc', '--no-centering', action='store_false', default=True) 
    parser.add_argument('-rsr', '--radial-size-reduction', type=float, default=0.)
    parser.add_argument('--fit', action='store_true', default=False)
    parser.add_argument('--set-xlim', nargs=2, type=float)
    parser.add_argument('--set-ylim', nargs=2, type=float)
    parser.add_argument('--set-gaps-xlim', nargs=2, type=float)
    parser.add_argument('--set-gaps-ylim', nargs=2, type=float)
    parser.add_argument('--set-theoretical-key-ylim', nargs=2, type=float)
    parser.add_argument('--shell-slope', type=float, default = 0.12)
    parser.add_argument('--pole-slope', type=float, default = -0.20)
    parser.add_argument('--grid', action='store_true', default=False)
    parser.add_argument('--gaps', nargs='+', type=str, default=['gaps_FUJI.size','gaps_final.size'])
    parser.add_argument('--font-size', type=float, default=12)
    parser.add_argument('--short-magnet', action='store_true', default=False) 
    parser.add_argument('--fix-gaps-kapton-size', type=float, default =None)


    args = parser.parse_args()

    try: 
        gaps = [{'file name':gaps_file_name, 'DataFrame':pd.read_table(gaps_file_name)} for gaps_file_name in args.gaps]

        if args.fix_gaps_kapton_size is not None: 
            print("Take kapton layer into account, total kapton size", args.fix_gaps_kapton_size)
            for gap_dict in gaps:
                for key in gap_dict['DataFrame']:
                    if 'A' in key or 'B' in key or 'C' in key or 'D' in key:
                        gap_dict['DataFrame'][key]+=args.fix_gaps_kapton_size
        for gap_dict in gaps:
            gap_df = gap_dict['DataFrame']
            gap_df['Gap AVG'] = gap_df.filter(regex='Gap.*').mean(axis=1)
            for key in gap_df:
                print('computing interpolation function for:', key)
                gap_dict[key+'interp'] = scipy.interpolate.InterpolatedUnivariateSpline(gap_dict['DataFrame']['Y'], gap_dict['DataFrame'][key],k=1)
    except:
        print("Gaps data not found.")
        gaps = None
        pass

    plt.rcParams.update({'font.size':args.font_size})

    print("read coil sizes")
    coil_size_dicts = read_coil_size_dicts(args)
    #plot_coil_sizes(coil_size_dicts, args)
    average_coil_size_interp, av_xdata = get_average_coil_size_interp(coil_size_dicts,args,col_name='L+R',)
    mshimmed_av_interp, av_xdata = get_average_coil_size_interp(coil_size_dicts,args,col_name='L+R+mshim')
    average_coil_size_shim_rsr_interp, av_xdata = get_average_coil_size_interp(coil_size_dicts,args,col_name='L+R+mshim-rsr')
    print("plot coil sizes")
    plot_interpolated_coil_sizes(coil_size_dicts, args, av_interp = average_coil_size_interp, mshimmed_av_interp=mshimmed_av_interp, shimmed_rsr_av_interp=average_coil_size_shim_rsr_interp)
    plot_interpolated_shimmed_coil_sizes(coil_size_dicts, args, av_interp = average_coil_size_interp, mshimmed_av_interp=mshimmed_av_interp, shimmed_rsr_av_interp=average_coil_size_shim_rsr_interp)
    plot_interpolated_coil_pack_z(coil_size_dicts, args, mshimmed_av_interp, gaps)
    plot_interpolated_collar_gap_vs_coil_sec_len(coil_size_dicts, args, mshimmed_av_interp, gaps)
    plot_interpolated_collar_radius_vs_coil_pack_size(coil_size_dicts, args, mshimmed_av_interp, gaps)
    plot_interpolated_collar_gaps(coil_size_dicts, args, mshimmed_av_interp, gaps)
    plot_interpolated_pole_key_gaps(coil_size_dicts, args, mshimmed_av_interp, gaps)
    plot_interpolated_theoretical_load_key(coil_size_dicts, args, av_interp = average_coil_size_interp, mshimmed_av_interp=mshimmed_av_interp, shimmed_rsr_av_interp=average_coil_size_shim_rsr_interp, shell_slope = args.shell_slope, pole_slope = args.pole_slope)

    shimming_info = get_shimming_info(coil_size_dicts, gaps=gaps)
    print(shimming_info)
    print("Store average coil size to average_coil.size")
    avdata = np.c_[av_xdata, average_coil_size_interp(av_xdata)]
    np.savetxt('average_coil.size', avdata, header='Y\n Average coil size')
    print("Store shimmed coil size to average_shimmed_coil.size")
    avshimdata = np.c_[av_xdata, mshimmed_av_interp(av_xdata)]
    np.savetxt('average_shimmed_coil.size', avshimdata, header='Y\n Average shimmed coil size')

    if args.radial_size_reduction != 0:
        print("Store shimmed coil size with radial size reduction to average_shimmed_coil_rsr.size")
        np.savetxt('average_shimmed_coil_rsr.size', avshimrsrdata, header='Y\n Average shimmed rsr coil size')
        avshimrsrdata = np.c_[av_xdata, average_coil_size_shim_rsr_interp(av_xdata)]


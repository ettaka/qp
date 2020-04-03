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
        if 'mshim' in line: mshim = float(line.split()[1])*args.unit_scaling
        elif 'offset' in line: offset = float(line.split()[1])*args.unit_scaling
        elif 'rshim' in line: rshim = float(line.split()[1])*args.unit_scaling

    if args.rshim != None: rshim = args.rshim

    print "_________________________"
    print "filepath", filepath
    print "mid-plane shim:", mshim
    print "radial shim:", rshim
    print "offset", offset

    col_names = ['L+R' if 'L+R' in s else s for s in head[-1].strip('\r\n').split('\t')]
    print "column names:", col_names
    coil_size_dict['mshim'] = mshim
    coil_size_dict['rshim'] = rshim
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

    print "Add column L+R+mshim"
    add_to_raw_data('L+R+mshim', coil_size_dict, raw_data[:,col_inds['L+R']] + mshim)

    print "Add column L+R+mshim-rsr"
    add_to_raw_data('L+R+mshim-rsr', coil_size_dict, raw_data[:,col_inds['L+R']] + mshim - args.radial_size_reduction)

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
    ax.set_xlabel("Longitudinal location (m)")
    ax.set_ylabel("L+R ($\mu$m)")
    if av_interp != None:
        ax.plot(xdata* args.xunit_plot_scaling,av_interp(xdata)* args.yunit_plot_scaling,color='black', marker=markers[totnum],label='av', linewidth=3)
    #if shimmed_av_interp != None:
        #ax.plot(xdata* args.xunit_plot_scaling,shimmed_av_interp(xdata)* args.yunit_plot_scaling,'-r'+markers[totnum],label='av+mshim', linewidth=5)
    #if shimmed_rsr_av_interp != None and args.radial_size_reduction != 0:
        #ax.plot(xdata* args.xunit_plot_scaling,shimmed_rsr_av_interp(xdata)* args.yunit_plot_scaling,'-y'+markers[totnum],label='av+mshim-rsr', linewidth=5)
    ax.legend(loc=args.legend_location)
    if args.show_plot:
        plt.show()
    else:
        save_fig = 'coil_size.png'
        plt.savefig(save_fig)

def plot_interpolated_shimmed_coil_sizes(coil_size_dicts, args, av_interp, shimmed_av_interp=None, shimmed_rsr_av_interp=None):
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
    ax.set_xlabel("Longitudinal location (m)")
    ax.set_ylabel("L+R ($\mu$m)")
    if shimmed_av_interp != None:
        ax.plot(xdata* args.xunit_plot_scaling,shimmed_av_interp(xdata)* args.yunit_plot_scaling,color='black', marker=markers[totnum],label='av+mshim', linewidth=3)
    ax.legend(loc=args.legend_location)
    if args.show_plot:
        plt.show()
    else:
        save_fig = 'coil_size_shimmed.png'
        plt.savefig(save_fig)

def plot_interpolated_theoretical_load_key(coil_size_dicts, args, av_interp, shimmed_av_interp=None, shimmed_rsr_av_interp=None, shell_slope=0.12, pole_slope=-0.20):
    plt.cla()
    totnum = 0

    ax2 = ax.twinx()
    ax3 = ax.twinx()
    #if args.set_xlim != None: ax.set_xlim(args.set_xlim)
    #if args.set_ylim != None: ax.set_ylim(args.set_ylim)
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
        totnum = i + 1
    ax.set_xlabel("Longitudinal location (m)")
    ax.set_ylabel("Size w.r.t. minimum ($\mu$m)")
    ax2.set_ylabel("Azimuthal shell stress w.r.t. minimum (MPa)")
    ax3.set_ylabel("Azimuthal pole stress w.r.t. minimum (MPa)")
    ax3.spines["right"].set_position(("axes", 1.2))
    fig.subplots_adjust(right=0.75)
    if shimmed_av_interp != None:
        ax.plot(xdata* args.xunit_plot_scaling,2./np.pi*(np.max(shimmed_av_interp(xdata)) - shimmed_av_interp(xdata))* args.yunit_plot_scaling,color='black', marker=markers[totnum],label='Theoretical load key', linewidth=3)
        ax2.plot(xdata* args.xunit_plot_scaling, shell_slope * 2./np.pi*(np.max(shimmed_av_interp(xdata)) - shimmed_av_interp(xdata))* args.yunit_plot_scaling,color='black', marker=markers[totnum],label='Theoretical load key', linewidth=3)
        mn, mx = ax.get_ylim()
        ax2.set_ylim(shell_slope * mn, shell_slope * mx)
        ax3.set_ylim(pole_slope * mn, pole_slope * mx)
    ax.legend(loc=args.legend_location)
    if args.show_plot:
        plt.show()
    else:
        save_fig = 'theoretical_load_key.png'
        plt.savefig(save_fig)

def plot_interpolated_coil_pack_z(coil_size_dicts, args, av_interp, gaps=None):
    plt.cla()
    totnum = 0
    #if args.set_xlim != None: ax.set_xlim(args.set_xlim)
    #if args.set_ylim != None: ax.set_ylim(args.set_ylim)
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
    ax.plot(xdata * args.xunit_plot_scaling, (av_interp(xdata)+ mshim + rshim*math.pi/2.) * args.yunit_plot_scaling,color='red', marker=markers[totnum],label='Shimmed coil average sector length', linewidth=3)
    if gaps is not None:
        for j, gap_data in enumerate(gaps):
            gap_file_name = gap_data['file name']
            data_name = gap_file_name.replace('.size', '').replace('gaps','').replace('_','')
            gap_df = gap_data['DataFrame']
            totnum = totnum + 1
            theoretical_gap = 15.
            coilpack_to_pad = 10.
            gap_data_average = gap_df.filter(regex='Gap.*').mean(axis=1)
            xdata = (gap_df['Y']-coilpack_to_pad)/1000. * args.xunit_plot_scaling
            ydata = (gap_data_average-theoretical_gap)/1000. * args.yunit_plot_scaling
            ax.plot(xdata, ydata, '--'+colors[totnum]+markers[totnum],label='Average collar gap '+data_name, linewidth=3)

            if j == 0:
                for i in range(len(xdata)):
                    x0 = xdata[i]
                    y0 = ydata[i]+50
                    ax.annotate(str(i+1), (x0, y0),  ha='center', va='center')#, fontsize=annotation_font_size)
    ax.legend(loc=args.legend_location)
    if args.show_plot:
        plt.show()
    else:
        save_fig = 'coil_average_sector_length.png'
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

def get_shimming_info(coil_size_dicts, info_locations=None, location_length=.05, location_length_points=100., coilpackdiff = -100.):
    if info_locations == None:
        short_shell = 0.3415
        long_shell = 2. * short_shell
        l = 2*short_shell + 10*long_shell
        #LE = short_shell + long_shell
        #CE = l / 2.
        #RE = l - LE
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
                    "Coilpack e": []
                    }

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
            else:
                location_start = location-location_length/2.
                location_end = location+location_length/2.
                location_area = np.linspace(location_start, location_end, location_length_points)
                coil_size = 1e6*(np.average(coil_size_dict['interp'][col_inds['L+R']](location_area)))
                shimming_info['Location'].append(1e3*location)
                shimming_info['Loc len'].append(1e3*location_length)

            shimming_info['Coil'].append(coil_size_dict['filepath'].replace('.size',''))
            shimming_info['L+R'].append(coil_size)
            shimming_info['mshim'].append(mshim)
            shimming_info['rshim'].append(rshim)
            shimming_info['L+R+mshim'].append(coil_size+mshim)
            shimming_info['dR'].append(2./math.pi *(coil_size+mshim))
            shimming_info['Coilpack'].append(2./math.pi *(coil_size+mshim)+rshim)
            shimming_info['Coil average sector length'].append(coil_size+mshim+rshim*math.pi/2.)
            shimming_info['Coilpack e'].append(2./math.pi *(coil_size+mshim)+rshim+coilpackdiff)


    sort_values = ['Loc len', 'Location']
    df = pd.DataFrame(shimming_info)
    cols = ["Coil", "Location name", "Location", "Loc len", "L+R", "mshim", "L+R+mshim", "dR", "rshim", "Coilpack", "Coil average sector length", "Coilpack e"]

    df_out = df[cols].sort_values(sort_values)
    print df_out.to_string(float_format="{:0.0f}".format, index=False)

    avg_table = df.groupby(['Location name','Location'], as_index=False, group_keys=True).mean()
    avg_table['Coil'] = ['AVG' for dR in avg_table['dR']]
    print avg_table[cols].sort_values(sort_values).to_string(float_format="{:0.0f}".format, index=False)

    joined_table = df_out.append(avg_table)

    print joined_table[cols].sort_values(sort_values).to_string(float_format="{:0.0f}".format, index=False)

    return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot coil sizes')
    parser.add_argument('paths', nargs='+', type=str)
    parser.add_argument('-sp', '--show-plot', action='store_true', default=False) 
    parser.add_argument('-ll', '--legend-location', type=str, default='best')
    parser.add_argument('-us', '--unit-scaling', type=float, default=0.001)
    parser.add_argument('-xus', '--xunit-plot-scaling', type=float, default=1)
    parser.add_argument('-yus', '--yunit-plot-scaling', type=float, default=1e6)
    parser.add_argument('-rshim', type=float, default=None)
    parser.add_argument('-nc', '--no-centering', action='store_false', default=True) 
    parser.add_argument('-rsr', '--radial-size-reduction', type=float, default=0.)
    parser.add_argument('--fit', action='store_true', default=False)
    parser.add_argument('--set-xlim', nargs=2, type=float)
    parser.add_argument('--set-ylim', nargs=2, type=float)
    parser.add_argument('--shell-slope', type=float, default = 0.12)
    parser.add_argument('--pole-slope', type=float, default = -0.20)
    parser.add_argument('--grid', action='store_true', default=False)
    parser.add_argument('--gaps', nargs='+', type=str, default=['gaps_FUJI.size','gaps_final.size'])
    parser.add_argument('--font-size', type=float, default=30)


    args = parser.parse_args()

    try: 
        gaps = [{'file name':gaps_file_name, 'DataFrame':pd.read_table(gaps_file_name)} for gaps_file_name in args.gaps]
    except:
        print "Gaps data not found."
        gaps = None
        pass

    plt.rcParams.update({'font.size':args.font_size})

    print "read coil sizes"
    coil_size_dicts = read_coil_size_dicts(args)
    #plot_coil_sizes(coil_size_dicts, args)
    average_coil_size_interp, av_xdata = get_average_coil_size_interp(coil_size_dicts,args,col_name='L+R',)
    average_coil_size_shim_interp, av_xdata = get_average_coil_size_interp(coil_size_dicts,args,col_name='L+R+mshim')
    average_coil_size_shim_rsr_interp, av_xdata = get_average_coil_size_interp(coil_size_dicts,args,col_name='L+R+mshim-rsr')
    print "plot coil sizes"
    plot_interpolated_coil_sizes(coil_size_dicts, args, av_interp = average_coil_size_interp, shimmed_av_interp=average_coil_size_shim_interp, shimmed_rsr_av_interp=average_coil_size_shim_rsr_interp)
    plot_interpolated_shimmed_coil_sizes(coil_size_dicts, args, av_interp = average_coil_size_interp, shimmed_av_interp=average_coil_size_shim_interp, shimmed_rsr_av_interp=average_coil_size_shim_rsr_interp)
    plot_interpolated_coil_pack_z(coil_size_dicts, args, average_coil_size_interp, gaps)
    plot_interpolated_theoretical_load_key(coil_size_dicts, args, av_interp = average_coil_size_interp, shimmed_av_interp=average_coil_size_shim_interp, shimmed_rsr_av_interp=average_coil_size_shim_rsr_interp, shell_slope = args.shell_slope, pole_slope = args.pole_slope)

    shimming_info = get_shimming_info(coil_size_dicts)
    print shimming_info
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


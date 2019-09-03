#!env python

import numpy as np
#import matplotlib
import matplotlib.pyplot as plt
import time
#import pandas as pd
import argparse
import codecs

from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

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

    data_shape = np.shape(data_dict['raw_data'])
    
    if len(data_shape) == 1:
        data_dict['average'] = np.average(data_dict['raw_data'])
        data_dict['min'] = np.average(data_dict['raw_data'])
        data_dict['max'] = np.average(data_dict['raw_data'])
        data_dict['error'] = np.array([data_dict['average'] - data_dict['min'], data_dict['max'] - data_dict['average']])
    else:
        for i in range(data_shape[1]):
            if np.any(data_dict['raw_data'][:,i]!=0): use_cols.append(i)

        data_dict['average'] = np.average(data_dict['raw_data'][:,use_cols], axis=1)
        data_dict['min'] = np.min(data_dict['raw_data'][:,use_cols], axis=1)
        data_dict['max'] = np.max(data_dict['raw_data'][:,use_cols], axis=1)
        data_dict['error'] = np.array([data_dict['average'] - data_dict['min'], data_dict['max'] - data_dict['average']])

def create_pk_npk_dict(pk_npk_data_file):
    pk_npk_dict = {}
    try:
        fh = open(pk_npk_data_file)
        header_lines = 3
        with codecs.open(pk_npk_data_file) as f:
            pk_npk_raw_data = np.loadtxt(f, skiprows=header_lines)
        pk_npk_dict['raw_data'] = pk_npk_raw_data
        pk_npk_dict['npk-spole'] = pk_npk_raw_data[:,1]
        pk_npk_dict['npk-scyl'] = pk_npk_raw_data[:,3]
        pk_npk_dict['pk-spole'] = pk_npk_raw_data[:,5]
        pk_npk_dict['pk-scyl'] = pk_npk_raw_data[:,7]
    except:
        pk_npk_dict = None
        pass


    return pk_npk_dict

def fix_zero_cols_with_average(cols):
    zero_col_indx=[]
    average = np.zeros_like(cols[:,0])
    for i, col in enumerate(cols.T):
        if np.count_nonzero(col) == 0:
            zero_col_indx.append(i)
            print 'Column', i, 'detected as a zero column'
        else:
            average+=col
    average /= i
    for indx in zero_col_indx:
        cols[:,indx]=average

def create_data_dicts(filepath, args, coil_permutation=None):
    print "Creating data dictionaries."
    if coil_permutation == None: coil_permutation = [1,2,3,4]
    print "Permuting coils with", coil_permutation
    coil_permutation = [inx-1 for inx in coil_permutation]
    header_lines = 3
    with codecs.open(filepath) as f:
        head = [next(f).decode('utf-8', 'ignore') for x in xrange(header_lines)]

    col_names = head[2].strip('\r\n').split('\t')
    print len(col_names)

    with codecs.open(filepath) as f:
        raw_data = np.loadtxt(f, skiprows=header_lines)

    key_dict = {}
    shell_dict = {}
    coil_dict = {}
    key_dict['axis label'] = 'Loading Key Thickness (mm)'
    shell_dict['axis label'] = 'Shell Azimuthal Stress (MPa)'
    coil_dict['axis label'] = 'Pole Azimuthal Stress (MPa)'
    if len(col_names) == 9 or len(col_names) == 10:
        key_dict['raw_data'] = raw_data[:,0]
        key_dict['col_names'] = col_names[0]
        shell_dict['raw_data'] = raw_data[:,1:5]
        shell_dict['col_names'] = col_names[1:5]
        ypicker = map([5,6,7,8].__getitem__,coil_permutation)
        coil_dict['raw_data'] = raw_data[:,ypicker]
        coil_dict['col_names'] = map(col_names[5:9].__getitem__,coil_permutation)
    elif len(col_names) == 12:
        key_dict['raw_data'] = raw_data[:,0:4]
        key_dict['col_names'] = col_names[0:4]
        shell_dict['raw_data'] = raw_data[:,1+3:5+3]
        shell_dict['col_names'] = col_names[1+3:5+3]
        ypicker = map([5+3,6+3,7+3,8+3].__getitem__,coil_permutation)
        coil_dict['raw_data'] = raw_data[:,ypicker]
        coil_dict['col_names'] = map(col_names[5+3:9+3].__getitem__,coil_permutation)
    else:
        print "You have", str(col_names), "columns in your datafile! (it should be 9, 10 or 12)"
        exit()

    if args.fix_gauges_with_average:
        print "Fixing shell gauges with average value"
        fix_zero_cols_with_average(shell_dict['raw_data'])
        print "Fixing coil gauges with average value"
        fix_zero_cols_with_average(coil_dict['raw_data'])

    if args.remove_coil_deltas != None:
        print "removing coil deltas of points with indices:", args.remove_coil_deltas 
        for ind in args.remove_coil_deltas.split():
            index = int(ind) - 1
            if index < np.shape(coil_dict['raw_data'])[0]:
                for i in range(4):
                    dy = coil_dict['raw_data'][index,i] - coil_dict['raw_data'][index-1,i]
                    coil_dict['raw_data'][index:, i] -= dy

    if args.pick_only_last_points:
        coil_dict['raw_data']=np.array([coil_dict['raw_data'][-1, :]])
        key_dict['raw_data']=np.array([key_dict['raw_data'][-1]])
        shell_dict['raw_data']=np.array([shell_dict['raw_data'][-1, :]])

    compute_data_dict_avg_min_max(key_dict)
    compute_data_dict_avg_min_max(shell_dict)
    compute_data_dict_avg_min_max(coil_dict)

    return key_dict, shell_dict, coil_dict

def plot_tf(ax, times_called, filepath, args):

    tr_type = args.type
    pk_npk_file = args.pk_npk_file
    single_coils = args.single_coils
    no_plot_average = args.no_average
    no_plot_average_error = args.no_average_error
    neighbour_shell_averages = not args.no_neighbour_shell_averages
    coil_permutation = args.coil_permutation

    filebase = filepath.split('.txt')[0]

    key_dict, shell_dict, coil_dict = create_data_dicts(filepath, args, coil_permutation)

    if args.key_shell:
        xdict = key_dict
        ydict = shell_dict
        plotname = filebase + '_key_shell'
    elif args.key_pole:
        xdict = key_dict
        ydict = coil_dict
        plotname = filebase + '_key_pole'
    else:
        xdict = shell_dict
        ydict = coil_dict
        plotname = filebase

    if args.print_final_stresses:
        print "Final stresses:"
        print "---------------"
        print "\tmin\tmax\taverage"
        print "shell\t" + str(xdict['min'][-1]) + "\t" + str(xdict['max'][-1]) + "\t" + str(xdict['average'][-1])
        print "pole\t" + str(ydict['min'][-1])  + "\t" + str(ydict['max'][-1]) + "\t" + str(ydict['average'][-1])
        print "---------------"
        print "::filepath\tshell min\tshell max\tshell average\tpole min\tpole max\tpole average"
        print ":" + filepath + "\t" + str(xdict['min'][-1]) + "\t" + str(xdict['max'][-1]) + "\t" + str(xdict['average'][-1]) + "\t" + str(ydict['min'][-1])  + "\t" + str(ydict['max'][-1]) + "\t" + str(ydict['average'][-1])
        print "---------------"
 
    else:

        if args.remove_coil_deltas != None:
            plotname += "_remove_coil_deltas" + args.remove_coil_deltas.replace(' ', '_')

        xdata = xdict['average']
        ydata = ydict['average']
     
        pk_npk_dict = create_pk_npk_dict(pk_npk_file)

        colors = ['r','g','b','c','m','y','k']
        markers = ['o', '^', 'v', '<', '>', 's', 'p', '*', 'h', 'd', '1', '2', '3', '4']

        if times_called < 1 and pk_npk_dict != None:
            ax.plot(pk_npk_dict['pk-scyl'],pk_npk_dict['pk-spole'],'-bo',label='FEM3D PK', markersize=args.marker_size, linewidth=args.line_width)
            ax.plot(pk_npk_dict['npk-scyl'],pk_npk_dict['npk-spole'],'-bd',label='FEM3D NPK', markersize=args.marker_size, linewidth=args.line_width)

        if args.label_type == 'filename':
            data_label=filepath.replace('.txt','').replace('TRANSFER1_','')
        else:
            data_label='Meas. Av.'

        data_color=colors[times_called%len(colors)]
        data_marker=markers[times_called%len(markers)]

        linestyle = '--'
        if args.pick_only_last_points:
            linestyle = ''

        if no_plot_average:
            print "Plot average of all shell vs pole gauges."
            ax.plot(xdata,ydata,linestyle+data_marker, color=data_color,label=data_label, markersize=args.marker_size, linewidth=args.line_width)
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

                ax.plot(xdata, ydata,'--'+colors[i]+markers[i],label=label, markersize=args.marker_size, linewidth=args.line_width)

        if args.fit:
            lower_limit, upper_limit = tuple(float(s) for s in args.fit_range.split())
            lower_index = np.argwhere(xdata>lower_limit)[0][0]
            upper_index = np.argwhere(xdata<=upper_limit)[-1][0]
            fit_plot_xdata = np.insert(xdata, 0, [args.fit_plot_lower])
            #print xdata, ydata
            #print lower_limit, upper_limit
            #print lower_index, upper_index
            #print xdata, lower_limit, lower_index, upper_limit, upper_index   

            fit = np.poly1d(np.polyfit(xdata[lower_index:upper_index], ydata[lower_index:upper_index], deg=1))
            ax.plot(fit_plot_xdata,fit(fit_plot_xdata), color='black', linestyle='--', label=data_label+' fit', linewidth=args.line_width, markersize=args.marker_size)

            fit_label = 'Fitted initial thickness = {:2.1f}'.format(fit.r[0])+' mm\n'
            fit_label += 'Fitted slope = {:2.0f}'.format(fit[1]) + ' MPa/mm\n'
            #ax.text(0.1, 0.1, fit_label, transform = ax.transAxes)
            fitfilename = plotname + '.fit'
            fitfile = open(fitfilename, 'w')
            print "writing fit parameters to file: ", fitfilename
            fitfile.write(fit_label)


        ax.set_xlabel(xdict['axis label'])
        ax.set_ylabel(ydict['axis label'])


        #lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5,-.1), fancybox=True, shadow=True, ncol=2)
        #plt.savefig(plotname, bbox_extra_artists=(lgd,), bbox_inches='tight')
    return plotname

def set_ax_parameters(ax, args):
    name_suffix = ''
    ax.legend(loc=args.legend_location, numpoints=1)

    if args.no_xaxis:
        ax.get_xaxis().set_visible()
        name_suffix += '_no-xaxis'
    if args.no_yaxis:
        ax.get_yaxis().set_visible()
        name_suffix += '_no-yaxis'
    if args.no_xticklabels: 
        ax.set_xticklabels([])
        name_suffix += '_no-xticklabels'
    if args.no_yticklabels: 
        ax.set_yticklabels([])
        name_suffix += '_no-yticklabels'
    if args.no_xlabel: 
        ax.set_xlabel('')
        name_suffix += '_no-xlabel'
    if args.no_ylabel: 
        ax.set_ylabel('')
        name_suffix += '_no-ylabel'
    ax.grid()

    if args.set_xticks != None:
        xticks = [float(tic) for tic in args.set_xticks.split()]
        ax.set_xticks(xticks)

    if args.set_yticks != None:
        yticks = [float(tic) for tic in args.set_yticks.split()]
        ax.set_yticks(yticks)

    if args.set_xlim!= None:
        xlim= [float(rang) for rang in args.set_xlim.split()]
        ax.set_xlim(xlim)

    if args.set_ylim!= None:
        ylim= [float(rang) for rang in args.set_ylim.split()]
        ax.set_ylim(ylim)

    if args.legend_outside:
        lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5,-.1), fancybox=True, shadow=True, ncol=2, numpoints=1)
    else:
        lgd = ax.legend()

    return lgd, name_suffix


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
    parser.add_argument('--set-xticks', type=str)
    parser.add_argument('--set-yticks', type=str)
    parser.add_argument('--set-xlim', type=str)
    parser.add_argument('--set-ylim', type=str)
    parser.add_argument('--key-shell', action='store_true', default=False) 
    parser.add_argument('--key-pole', action='store_true', default=False) 
    parser.add_argument('--remove-coil-deltas', type=str)
    parser.add_argument('--label-type', type=str)
    parser.add_argument('--print-final-stresses', action='store_true', default=False) 
    parser.add_argument('--fit', action='store_true', default=False)
    parser.add_argument('--fit-range', type=str)
    parser.add_argument('--fit-plot-lower', type=float, default=13.2)
    parser.add_argument('--pick-only-last-points', action='store_true', default=False) 
    parser.add_argument('--legend-outside', action='store_true', default=False) 
    parser.add_argument('--no-xaxis', action='store_true', default=False) 
    parser.add_argument('--no-yaxis', action='store_true', default=False) 
    parser.add_argument('--no-xticklabels', action='store_true', default=False) 
    parser.add_argument('--no-yticklabels', action='store_true', default=False) 
    parser.add_argument('--no-xlabel', action='store_true', default=False) 
    parser.add_argument('--no-ylabel', action='store_true', default=False) 
    parser.add_argument('--font-size', type=float, default=16)
    parser.add_argument('--marker-size', type=float, default=8)
    parser.add_argument('--line-width', type=float, default=1.5)
    parser.add_argument('--fig-height', type=float, default=8)
    parser.add_argument('--fig-width', type=float, default=8)
    parser.add_argument('--title', type=str, default='')
    parser.add_argument('--image-name', type=str, default='')
    parser.add_argument('-fgwa', '--fix-gauges-with-average', action='store_true', default=False) 


    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.autoscale(enable=True, axis='y', tight=True)

    args = parser.parse_args()
    fig.set_figheight(args.fig_height)
    fig.set_figwidth(args.fig_width)
    plt.rcParams.update({'font.size':args.font_size})
    plt.title(args.title)

    paths = args.paths
    plotnames = []
    for i, filepath in enumerate(paths):
        #plt.cla()
        print "filepath", filepath
        plotname = plot_tf(ax, i, filepath, args)
        plotnames.append(plotname)
    
    lgd, name_suffix = set_ax_parameters(ax, args)

    if not args.print_final_stresses:
        if args.show_plot:
            plt.show()
        else:
            imagename = '_'.join(plotnames)
            if 'TRANSFER1' in imagename:
                imagename = '_'.join(plotnames).replace('TRANSFER1_MQXF','')
                imagename = 'TRANSFER1_MQXF' + imagename
            imagename += name_suffix
            if args.image_name != '': imagename = args.image_name
            print "creating image file", imagename+'.png'
            if args.legend_outside:
                plt.savefig(imagename + '.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
            else:
                plt.savefig(imagename + '.png', bbox_inches='tight')



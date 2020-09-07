#!env python

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import time
#import pandas as pd
import argparse
import codecs
import pandas as pd
import re
import os

from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

from scipy.stats import linregress

from channel_utils import create_channel_dict_list
from ansys_parse_utils import parse_ansys_2d_files
from ansys_parse_utils import parse_ansys_3d_files
from ansys_parse_utils import parse_ansys_2d_master_to_master
from ansys_parse_utils import octname_to_Qname

markers = ['o', '^', 'v', '<', '>', 's', 'p', '*', 'h', 'd', '1', '2', '3', '4']

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
        try:
            data_dict['min'] = np.min(data_dict['raw_data'][:,use_cols], axis=1)
        except:
            data_dict['min'] = data_dict['raw_data'][:,0]
            
        try:
            data_dict['max'] = np.max(data_dict['raw_data'][:,use_cols], axis=1)
        except:
            data_dict['max'] = data_dict['raw_data'][:,0]
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
            print('Column', i, 'detected as a zero column')
        else:
            average+=col
    average /= i
    for indx in zero_col_indx:
        cols[:,indx]=average

def fill_header(header_line):
    for i, string in enumerate(header_line):
        if pd.isnull(string): header_line[i] = ''

        header_line[i] = re.sub(r"\(.*\)", "", header_line[i])

        if i>0 and header_line[i]=='':
            header_line[i] = header_line[i-1]

def replace_header_strings(header, replace_dict):
    for replace_word in replace_dict:
        for i, element in enumerate(header):
            header[i] = element.replace(replace_word, replace_dict[replace_word])

def load_csv(filepath, times_called, args):
    header_lines = 3
    with codecs.open(filepath) as f:
        head = [next(f) for x in range(header_lines)]

    for i,line in enumerate(head):
        if 'date' in line.lower() and 'time' in line.lower():
            header = i

    df = pd.read_csv(filepath)
    fill_header(df.iloc[0])
    fill_header(df.iloc[1])
    new_header = df.iloc[0] + df.iloc[1]
    replace_dict = {
            'Shell Stress ': '',
            'Coil Stress ': '',
            'Rod Stress ': ''
            }
    replace_header_strings(new_header, replace_dict)

    df = df[2:]
    df.columns = new_header

    df['average key'] = df.filter(regex='Keys size.*').astype(float).mean(axis=1)
    
    df_out = df
    if args.all_points:
        if not args.bladders:
            df_out = df[~df['Operation'].str.contains('Bladder')]
    elif args.bladders:
        df_out = df[df['Operation'] == 'Bladder']
    else:

        if args.operations is not None:
            not_ignore_regex = '|'.join(args.operations)
            df_out = df_out[(df['Operation'].str.contains(not_ignore_regex))] 
        if args.no_mixed_keys:
            df_out = df_out[~df_out['Operation'].str.contains('Mixed')]
        if args.no_idle_points:
            df_out = df_out[~df_out['Operation'].str.contains('Idling')]
        if args.operations_ignore is not None:
            for label in args.operations_ignore:
                df_out = df_out[~df_out['Operation'].str.contains(label)]

    if args.use_fibre != None and args.use_fibre[times_called]==1:
        keys_data = df_out.filter(regex='Keys.*$|C...T_FBG$|SH.T')
    else:
        keys_data = df_out.filter(regex='Keys.*$|C...T$|SH.T')


    raw_data = keys_data.astype(float).values
    col_names = keys_data.columns
    channel_dict_list = create_channel_dict_list(col_names, col_names, old_format=False)

    print(len(col_names))

    return raw_data, col_names, df_out, channel_dict_list
 
def load_txt(filepath, args):
    header_lines = 3
    with codecs.open(filepath) as f:
        head = [next(f) for x in range(header_lines)]

    col_names = head[2].strip('\r\n').split('\t')
    print(len(col_names))

    with codecs.open(filepath) as f:
        raw_data = np.loadtxt(f, skiprows=header_lines)

    return raw_data, col_names


def create_data_dicts(filepath, times_called, file_extension, args, coil_permutation=None):

    print("Creating data dictionaries.")
    if coil_permutation == None: coil_permutation = [1,2,3,4]
    print("Permuting coils with", coil_permutation)
    coil_permutation = [inx-1 for inx in coil_permutation]

    if file_extension == '.txt':
        raw_data, col_names = load_txt(filepath, args)
        if not args.TF2: raw_data = raw_data[:-1,:]
        df = None
        channel_dict_list = None
    elif file_extension == '.csv':
        raw_data, col_names, df, channel_dict_list = load_csv(filepath, times_called, args)

    key_dict = {}
    shell_dict = {}
    coil_dict = {}
    key_dict['axis label'] = 'Loading Key Thickness (mm)'
    shell_dict['axis label'] = 'Shell Azimuthal Stress (MPa)'
    coil_dict['axis label'] = 'Pole Azimuthal Stress (MPa)'
    if len(col_names) == 9 or len(col_names) == 10:
        key_dict['raw_data'] = raw_data[:,0]
        key_dict['col_names'] = col_names[0]
        key_dict['channel_dict_list'] = None
        shell_dict['raw_data'] = raw_data[:,1:5]
        shell_dict['col_names'] = col_names[1:5]
        ypicker = list(map([5,6,7,8].__getitem__,coil_permutation))
        coil_dict['raw_data'] = raw_data[:,ypicker]
        coil_dict['col_names'] = list(map(col_names[5:9].__getitem__,coil_permutation))
        if channel_dict_list is not None:
            shell_dict['channel_dict_list'] = channel_dict_list[1:5]
            coil_dict['channel_dict_list'] = list(map(channel_dict_list[5:9].__getitem__,coil_permutation))
        else:
            shell_dict['channel_dict_list'] = None
            coil_dict['channel_dict_list'] = None
    elif len(col_names) == 12:
        key_dict['raw_data'] = raw_data[:,0:4]
        key_dict['col_names'] = col_names[0:4]
        key_dict['channel_dict_list'] = None
        shell_dict['raw_data'] = raw_data[:,1+3:5+3]
        shell_dict['col_names'] = col_names[1+3:5+3]
        ypicker = list(map([5+3,6+3,7+3,8+3].__getitem__,coil_permutation))
        coil_dict['raw_data'] = raw_data[:,ypicker]
        coil_dict['col_names'] = list(map(col_names[5+3:9+3].__getitem__,coil_permutation))
        if channel_dict_list is not None:
            shell_dict['channel_dict_list'] = channel_dict_list[1+3:5+3]
            coil_dict['channel_dict_list'] = list(map(channel_dict_list[5+3:9+3].__getitem__,coil_permutation))
        else:
            shell_dict['channel_dict_list'] = None
            coil_dict['channel_dict_list'] = None
    else:
        print("You have", str(col_names), "columns in your datafile! (it should be 9, 10 or 12)")
        exit()

    if args.fix_gauges_with_average:
        print("Fixing shell gauges with average value")
        fix_zero_cols_with_average(shell_dict['raw_data'])
        print("Fixing coil gauges with average value")
        fix_zero_cols_with_average(coil_dict['raw_data'])

    if args.remove_coil_deltas != None:
        print("removing coil deltas of points with indices:", args.remove_coil_deltas) 
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

    return key_dict, shell_dict, coil_dict, df

def plot_ansys_data(ax, args):
    ansys_file_data_list = parse_ansys_2d_files(args)
    if ansys_file_data_list is not None:
        data_by_parent={}
        parent_names = []
        for i, data in enumerate(ansys_file_data_list):
            parent_name = data['parent name']
            if not parent_name in parent_names:
                parent_names.append(parent_name)
            if not parent_name in data_by_parent:
                data_by_parent[parent_name] = {}
                data_by_parent[parent_name]['children'] = []

            data_by_parent[parent_name]['children'].append(data)

        for parent_name in parent_names:
            data_by_parent[parent_name]['scyl'] = {}
            data_by_parent[parent_name]['spole'] = {}
            data_by_parent[parent_name]['interf'] = {}
            scyl = []
            spole = []
            interf = []
            names = []
            for child in data_by_parent[parent_name]['children']:
                df = child['DataFrame']
                df = df[~df['name'].str.contains('cur')]
                if not args.TF2:
                    df = df[~df['name'].str.contains('cold')]
                scyl.append(df['scyl'])
                spole.append(df['spole'])
                interf.append(df['interf'])
                names.append(child['name'])

            data_by_parent[parent_name]['scyl']['raw_data'] = np.array(scyl).transpose()
            data_by_parent[parent_name]['spole']['raw_data'] = np.array(spole).transpose()
            data_by_parent[parent_name]['interf']['raw_data'] = np.array(interf).transpose()
            data_by_parent[parent_name]['names'] = names
            scyl_dict = data_by_parent[parent_name]['scyl']
            spole_dict = data_by_parent[parent_name]['scyl']
            interf_dict = data_by_parent[parent_name]['interf']
            names = data_by_parent[parent_name]['names']

            #def compute_data_dict_avg_min_max(data_dict):
            
        for i,parent_name in enumerate(parent_names):
            scyl = data_by_parent[parent_name]['scyl']
            spole = data_by_parent[parent_name]['spole']
            interf = data_by_parent[parent_name]['interf']
            names = data_by_parent[parent_name]['names']

            compute_data_dict_avg_min_max(scyl)
            compute_data_dict_avg_min_max(spole)
            compute_data_dict_avg_min_max(interf)

            data_color=args.curve_colors[i%len(args.curve_colors)]

            #ax.plot(scyl['average'], spole['average'], marker='d', markersize=args.marker_size, label=parent_name)
            if args.key_pole:
                if args.ansys_single_coils:
                    for j in range(4):
                        name = names[j]
                        name = octname_to_Qname(name)
                        data_color=args.curve_colors[j%len(args.curve_colors)]
                        data_marker=markers[j%len(markers)]
                        ax.plot(13.+interf['average']/1000., spole['raw_data'][:,j], marker=data_marker, markersize=3, label=name, color=data_color, linewidth=1, linestyle='-')
                else:
                    ax.plot(13.+interf['average']/1000., spole['average'], marker='d', markersize=1, label=parent_name, color=data_color, linewidth=2.5)
                    _ = make_error_boxes(ax, 13.+interf['average']/1000., spole['average'], interf['error']/1000., spole['error'], facecolor='b', edgecolor='None', alpha=0.3)
            elif args.key_shell:
                if args.ansys_single_coils:
                    for j in range(4):
                        name = names[j]
                        name = octname_to_Qname(name)
                        data_color=args.curve_colors[j%len(args.curve_colors)]
                        data_marker=markers[j%len(markers)]
                        ax.plot(13.+interf['average']/1000., scyl['raw_data'][:,j], marker=data_marker, markersize=3, label=name, color=data_color, linewidth=1, linestyle='-')
                else:
                    ax.plot(13.+interf['average']/1000., scyl['average'], marker='d', markersize=1, label=parent_name, color=data_color, linewidth=2.5)
                    _ = make_error_boxes(ax, 13.+interf['average']/1000., scyl['average'], interf['error']/1000., scyl['error'], facecolor='b', edgecolor='None', alpha=0.3)
            else:
                if args.ansys_single_coils:
                    for j in range(4):
                        name = names[j]
                        name = octname_to_Qname(name)
                        data_color=args.curve_colors[j%len(args.curve_colors)]
                        data_marker=markers[j%len(markers)]
                        ax.plot(scyl['raw_data'][:,j], spole['raw_data'][:,j], marker=data_marker, markersize=3, label=name, color=data_color, linewidth=1, linestyle='-')
                else:
                    ax.plot(scyl['average'], spole['average'], marker='d', markersize=1, label=parent_name, color=data_color, linewidth=2.5)
                    _ = make_error_boxes(ax, scyl['average'], spole['average'], scyl['error'], spole['error'], facecolor='b', edgecolor='None', alpha=0.3)

        for i, data in enumerate(ansys_file_data_list):
            name = data['name']
            df = data['DataFrame']

            #ax.plot(df['scyl'], df['spole'], marker='d', markersize=args.marker_size, label=name)

    ansys3d_file_data_list = parse_ansys_3d_files(args)
    if ansys3d_file_data_list is not None:
        for i, data in enumerate(ansys3d_file_data_list):
            name = data['name']

            df = data['DataFrame']
            ax.plot([0]+list(df['scyl']), [0]+list(df['spole']), marker='d', markersize=args.marker_size, label=name)

def plot_ansys_bladders(ax, args):
    ansys_file_data_list = parse_ansys_2d_master_to_master(args)
    
    if ansys_file_data_list is not None:
        for i, data in enumerate(ansys_file_data_list):
            name = data['name']
            df = data['DataFrame']
            xdata = 13.+df['bdisp1']/1e3
            ydata =  df['bpres']
            ansplot_name = 'interf. ' + name.split('i')[1].split('_')[0] + r' $\mu$m'
            marker=markers[0]
            if 'q1' in name: 
                ansplot_name += ' quadrant'
                marker=markers[1]
            else:
                ansplot_name += ' all'
            ax.plot(xdata, ydata, marker=marker, markersize=args.marker_size, label=ansplot_name)

def ax_plot(*args, **kwargs):
    if 'parseargs' in kwargs:
        parseargs = kwargs['parseargs']
        if not parseargs.no_measurements_plot:
            kwargs.pop('parseargs')
            ax.plot(*args, **kwargs)

def plot_tf(ax, times_called, filepath, args):
    tr_type = args.type
    pk_npk_file = args.pk_npk_file
    single_coils = args.single_coils
    no_plot_average = args.no_average
    no_plot_average_error = args.no_average_error
    neighbour_shell_averages = not args.no_neighbour_shell_averages
    coil_permutation = args.coil_permutation

    if '.csv' in filepath:
        file_extension = '.csv'
    else:
        file_extension = '.txt'

    filebase = filepath.split(file_extension)[0]
    key_dict, shell_dict, coil_dict, df = create_data_dicts(filepath, times_called, file_extension, args, coil_permutation)

    if args.bladders:
        key_cols = [col for col in df.keys() if 'Keys size' in col]
        for col in key_cols:
            df[col] = pd.to_numeric(df[col])

        bladder_cols = [col for col in df.keys() if 'Bladder' in col]
        for col in bladder_cols:
            df[col] = pd.to_numeric(df[col])
        xaxis = df[key_cols].mean(axis=1)
        yaxis = df[bladder_cols].max(axis=1)
        label = filebase
        ax_plot(xaxis, yaxis, marker='d', markersize=args.marker_size, label=label, linestyle = "None", parseargs=args)
        ax.set_xlabel('Key size (mm)')
        ax.set_ylabel('Bladder Pressure (Bars)')
        plotname = filebase + '_bladders'
        return plotname
    if args.sg_vs_fbg and df is not None:
        plt.close()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.autoscale(enable=True, axis='y', tight=True)

        FBG_list = [key for key in df if "FBG" in key]
        for key in df:
            df[key] = pd.to_numeric(df[key], errors='ignore')
        for fbg_key in FBG_list:
            ax.cla()
            key = fbg_key.replace('_FBG','')
            df_not_nan = df.dropna(subset=[key,fbg_key])
            if len(df_not_nan) == 0: df_not_nan=df.fillna(0)
            #if args.fit_range is not None: 
            #    lower_limit, upper_limit = tuple(float(s) for s in args.fit_range.split())
            #    df_not_nan = df_not_nan[df_not_nan[key] > lower_limit]
            #    df_not_nan = df_not_nan[df_not_nan[key] < upper_limit]
            if args.fit_point_range is not None:
                lower_limit, upper_limit = tuple(int(s) for s in args.fit_point_range.split())
                print("upper/lower", lower_limit, upper_limit)
                df_fit = df_not_nan[lower_limit:upper_limit]
            else:
                df_fit = df_not_nan
            ax.plot(df_not_nan[key], df_not_nan[fbg_key], marker='d', markersize=args.marker_size, label='__nolegend__', linestyle = "None")

            try:
                coef, cov = np.polyfit(df_fit[key], df_fit[fbg_key], 1, cov=True)
                slope, intercept, r_value, p_value, std_err = linregress(df_fit[key], df_fit[fbg_key])
            except:
                coef = (0,0)
                cov = ((0,0),(0,0))
                r_value = 0
                pass
            sigma = np.sqrt(cov)
            fit = np.poly1d(coef)
            ax.plot(df_not_nan[key], fit(df_not_nan[key]), label='y=({:0.1f}$\pm${:0.1f})x + ({:0.0f}$\pm${:.0f})\n$R^2$={:.2f}'.format(coef[0],sigma[0,0], coef[1],sigma[1,1], r_value**2.), marker="None", color='black', linestyle='-')  
            ax.plot(df_fit[key], df_fit[fbg_key], color='red', linestyle="None", marker='d', markersize=args.marker_size, label='__nolegend__')  
            ax.set_xlabel(key)
            ax.set_ylabel(fbg_key)
            plt.title(args.title)
            imagename = filepath+'_sg-vs-fbg_' + key
            imagename = imagename.replace('.','_').replace(' ','_')
            print("creating image file", imagename+'.png')
            if args.legend_outside:
                lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5,-.13), fancybox=True, shadow=True, ncol=2)
                plt.savefig(imagename + '.png', bbox_extra_artists=(lgd,), bbox_inches='tight', numpoints=1, dpi=200)
            else:
                lgd = ax.legend(loc='best')
                plt.savefig(imagename + '.png', bbox_inches='tight', numpoints=1, dpi=200)
        exit()

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
        print("Final stresses:")
        print("---------------")
        print("\tmin\tmax\taverage")
        print("shell\t" + str(xdict['min'][-1]) + "\t" + str(xdict['max'][-1]) + "\t" + str(xdict['average'][-1]))
        print("pole\t" + str(ydict['min'][-1])  + "\t" + str(ydict['max'][-1]) + "\t" + str(ydict['average'][-1]))
        print("---------------")
        print("::filepath\tshell min\tshell max\tshell average\tpole min\tpole max\tpole average")
        print(":" + filepath + "\t" + str(xdict['min'][-1]) + "\t" + str(xdict['max'][-1]) + "\t" + str(xdict['average'][-1]) + "\t" + str(ydict['min'][-1])  + "\t" + str(ydict['max'][-1]) + "\t" + str(ydict['average'][-1]))
        print("---------------")
 
    else:

        if args.remove_coil_deltas != None:
            plotname += "_remove_coil_deltas" + args.remove_coil_deltas.replace(' ', '_')

        xdata = xdict['average']
        ydata = ydict['average']
     
        pk_npk_dict = create_pk_npk_dict(pk_npk_file)

        colors = args.curve_colors

        if args.curve_markers==None:
            markers = ['o', '^', 'v', '<', '>', 's', 'p', '*', 'h', 'd', '1', '2', '3', '4']
        else:
            markers = args.curve_markers

        if times_called < 1 and pk_npk_dict != None:
            if args.no_pk_legend:
                pklabel = ''
                npklabel = ''
            else:
                pklabel = 'PK'
                npklabel = 'NPK'
            ax.plot(pk_npk_dict['pk-scyl'],pk_npk_dict['pk-spole'],'-bo',label=pklabel, markersize=args.marker_size, linewidth=args.line_width)
            ax.plot(pk_npk_dict['npk-scyl'],pk_npk_dict['npk-spole'],'-bd',label=npklabel, markersize=args.marker_size, linewidth=args.line_width)

        if args.label_type == 'filename':
            data_label=filepath.replace('.txt','').replace('TRANSFER1_','').replace('.csv','')
        elif args.label_type == 'imagename':
            data_label=args.image_name
        elif args.meas_legend_label != '':
            data_label=args.meas_legend_label
        else:
            data_label='Meas. Av.'

        data_color=colors[times_called%len(colors)]
        data_marker=markers[times_called%len(markers)]

        linestyle = '--'
        if args.pick_only_last_points:
            linestyle = ''

        if no_plot_average:
            print("Plot average of all shell vs pole gauges.")
            ax_plot(xdata,ydata,linestyle+data_marker, color=data_color,label=data_label, markersize=args.marker_size, linewidth=args.line_width, parseargs=args)
            if no_plot_average_error:
                _ = make_error_boxes(ax, xdata, ydata, xdict['error'], ydict['error'], facecolor='g', edgecolor='None', alpha=0.3)

        if single_coils:
            ax.set_title(data_label)
            print("Plot single pole gauges", end=' ')
            for i in range(4):
                if neighbour_shell_averages:
                    if i==0:print("vs single shell gauges.")
                    #xdata = xdict['raw_data'][:,i]
                    xdata = np.mean(xdict['raw_data'][:,1:4],axis=1)
                    ydata = ydict['raw_data'][:,i]
                    if xdict['channel_dict_list'] is not None:
                        labelx = xdict['channel_dict_list'][i]['location']
                    else:
                        labelx = xdict['col_names'][i]
                    if ydict['channel_dict_list'] is not None:
                        labely = ydict['channel_dict_list'][i]['location']
                    else:
                        labely = ydict['col_names'][i]
                    label = labelx+'-'+labely
                else:
                    if i==0:print("vs single shell gauges.")
                    nof_cols = np.shape(xdict['raw_data'])[1]
                    #xdata = (xdict['raw_data'][:,i] + xdict['raw_data'][:,(i+1)%nof_cols])/2.
                    xdata = np.mean(xdict['raw_data'][:,1:4],axis=1)
                    ydata = ydict['raw_data'][:,i]
                    if xdict['channel_dict_list'] is not None:
                        labelx = xdict['channel_dict_list'][i]['location']+xdict['channel_dict_list'][(i+1)%nof_cols]['location']
                    else:
                        labelx = xdict['col_names'][i]+xdict['col_names'][(i+1)%nof_cols]
                    if ydict['channel_dict_list'] is not None:
                        labely = ydict['channel_dict_list'][i]['location']
                    else:
                        labely = ydict['col_names'][i]
                    label = labelx+'-'+labely

                label = label.replace('Shell ','').replace('Keys size ','')

                ax_plot(xdata, ydata,'--'+colors[i]+markers[i],label=label, markersize=args.marker_size, linewidth=args.line_width, parseargs=args)

        try:
            if args.fit: fit = fit_data(ax, xdata, ydata, args.fit_range, data_label+' fit', args)
            if args.fit2: fit = fit_data(ax, xdata, ydata, args.fit2_range, data_label+' fit 2', args)

            if args.fit or args.fit2:
                if args.key_shell or args.key_pole:
                    fit_text = 'Fitted initial thickness = {:2.2f}'.format(fit.r[0])+' mm\n'
                    fit_text += 'Fitted slope = {:2.0f}'.format(fit[1]) + ' MPa/mm\n'
                else: 
                    fit_text = 'Fitted initial stress = {:2.2f}'.format(fit.r[0])+' MPa\n'
                    fit_text += 'Fitted slope = {:2.2f}'.format(fit[1]) + ' MPa/MPa\n'
                fitfilename = plotname + '.fit'
                fitfile = open(fitfilename, 'w')
                print("writing fit parameters to file: ", fitfilename)
                fitfile.write(fit_text)
        except:
            print("fitting failed!")


        ax.set_xlabel(xdict['axis label'])
        ax.set_ylabel(ydict['axis label'])

        if args.annotate_points:
            if df is not None:
                values = df['average key'].values
                value_counts = df['average key'].value_counts()
                annotations = df.index
            else:
                annotations = list(range(len(xdata)))

            if args.annotation_list:
                textlist = '\n'.join(annotations)
                x0 = min(ax.get_xlim())
                y0 = max(ax.get_ylim())
                ax.text(x0, y0, textlist)
            else:
                if args.bladders:
                    annotation_shift_percent = 0.02
                    step_percent = 0.02
                    annotation_font_size = 25
                else:
                    annotation_shift_percent = 0.05
                    step_percent = 0.05
                    annotation_font_size = 16
                value_count_dict = {}
                ylim = ax.get_ylim()
                yrange = np.max(ylim) - np.min(ylim)
                step = step_percent*yrange
                annotation_shift = annotation_shift_percent*yrange
                for i, txt in enumerate(annotations):
                    if values[i] not in value_count_dict:
                        value_count_dict[values[i]] = 0
                    value_count_dict[values[i]] += 1

                    highest_y = np.max(ydata[np.argwhere(xdata==xdata[i]).transpose()[0]])
                    x0 = xdata[i]
                    y0 = highest_y +annotation_shift+ step * (value_counts[values[i]]-value_count_dict[values[i]])
                    #y0 = annotation_shift + min(ax.get_ylim()) + step * (value_counts[values[i]]-value_count_dict[values[i]])
                    ax.text(x0, y0, str(txt), ha='center', va='center', fontsize=annotation_font_size)



        #lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5,-.1), fancybox=True, shadow=True, ncol=2)
        #plt.savefig(plotname, bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=200)
    return plotname

def fit_data(ax, xdata, ydata, fit_range, fit_label, args):
    lower_limit, upper_limit = tuple(float(s) for s in fit_range.split())
    lower_index = np.argwhere(xdata>lower_limit)[0][0]
    upper_index = np.argwhere(xdata<=upper_limit)[-1][0]
    fit_plot_xdata = np.insert(xdata, 0, [args.fit_plot_lower])
    fit = np.poly1d(np.polyfit(xdata[lower_index:upper_index], ydata[lower_index:upper_index], deg=1))
    ax.plot(fit_plot_xdata,fit(fit_plot_xdata), color='black', linestyle='--', label=fit_label, linewidth=args.line_width, markersize=args.marker_size)
    return fit
 
def set_ax_parameters(ax, args):
    name_suffix = ''
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

#    if args.legend_outside:
#        lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5,-.1), fancybox=True, shadow=True, ncol=2, numpoints=1)
#    else:
#        lgd = ax.legend(loc='best')

    if args.legend_location == 'right outside':
        lgd = ax.legend(loc='upper left', bbox_to_anchor=(1.04,1), fancybox=True, shadow=True, numpoints=1)
    elif args.legend_location == 'bottom outside':
        lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5,-.18), fancybox=True, shadow=True, ncol=args.legend_ncol, numpoints=1)
    elif args.legend_location == 'nolegend':
        lgd = None
    else:
        lgd = ax.legend(loc=args.legend_location, numpoints=1)

    return lgd, name_suffix


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot transfer function')
    parser.add_argument('paths', nargs='+', type=str)
    parser.add_argument('-t', '--type', type=int, default=1) 
    parser.add_argument('-pk', '--pk-npk-file', type=str, default=None) 
    parser.add_argument('--ansys-2d-files', type=str, nargs='+') 
    parser.add_argument('--ansys-3d-files', type=str, nargs='+') 
    parser.add_argument('--ansys-2d-bladder-files', type=str, nargs='+') 
    parser.add_argument('--ansys-show-x', action='store_true', default=False)
    parser.add_argument('--ansys-single-coils', action='store_true', default=False)
    parser.add_argument('-s', '--single-coils', action='store_true', default=False) 
    parser.add_argument('-sp', '--show-plot', action='store_true', default=False) 
    parser.add_argument('-nmp', '--no-measurements-plot', action='store_true', default=False) 
    parser.add_argument('--TF2', action='store_true', default=False) 
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
    parser.add_argument('--fit-point-range', type=str)
    parser.add_argument('--fit2', action='store_true', default=False)
    parser.add_argument('--fit2-range', type=str)
    parser.add_argument('--fit-plot-lower', type=float, default=13.2)
    parser.add_argument('--pick-only-last-points', action='store_true', default=False) 
    parser.add_argument('--legend-outside', action='store_true', default=False) 
    parser.add_argument('--no-xaxis', action='store_true', default=False) 
    parser.add_argument('--no-yaxis', action='store_true', default=False) 
    parser.add_argument('--no-xticklabels', action='store_true', default=False) 
    parser.add_argument('--no-yticklabels', action='store_true', default=False) 
    parser.add_argument('--no-xlabel', action='store_true', default=False) 
    parser.add_argument('--no-ylabel', action='store_true', default=False) 
    parser.add_argument('--font-size', type=float, default=10)
    parser.add_argument('--legend-font-size', type=float, default=10)
    parser.add_argument('--marker-size', type=float, default=4)
    parser.add_argument('--line-width', type=float, default=1.5)
    parser.add_argument('--fig-height', type=float, default=3)
    parser.add_argument('--fig-width', type=float, default=4)
    parser.add_argument('--title', type=str, default='')
    parser.add_argument('--image-name', type=str, default='')
    parser.add_argument('-fgwa', '--fix-gauges-with-average', action='store_true', default=False) 
    parser.add_argument('--no-pk-legend', action='store_true', default=False) 
    parser.add_argument('--no-legend', action='store_true', default=False) 
    parser.add_argument('--legend-ncol', type=int, default=2)
    parser.add_argument('--meas-legend-label', type=str, default='')
    parser.add_argument('--plot-style', type=str, default='')
    parser.add_argument('--curve-colors', nargs='+', type=str, default=None)
    parser.add_argument('--curve-markers', nargs='+', type=str, default=None)
    parser.add_argument('--operations-ignore', nargs='+', type=str, default=None)
    parser.add_argument('--operations', nargs='+', type=str, default=['Idling','Key','Initial'])
    parser.add_argument('-uf', '--use-fibre', nargs='+', type=int) 
    parser.add_argument('--no-mixed-keys', action='store_true', default=False) 
    parser.add_argument('--no-idle-points', action='store_true', default=False) 
    parser.add_argument('--all-points', action='store_true', default=False) 
    parser.add_argument('--annotate-points', action='store_true', default=False) 
    parser.add_argument('--annotation-list', action='store_true', default=False) 
    parser.add_argument('--annotation-text', action='store_true', default=False) 
    parser.add_argument('--bladders', action='store_true', default=False) 
    parser.add_argument('--sg-vs-fbg', action='store_true', default=False) 
    parser.add_argument('--delete-dev-higher', type=float, default=1.)
    parser.add_argument('--errorbar-capsize', type=float, default=2.)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.autoscale(enable=True, axis='y', tight=True)

    args = parser.parse_args()

    if args.curve_colors is None:
        args.curve_colors = ['r','g','k','c','m','y','b']

    if args.plot_style == 'TF paper':
        args.set_xlim = "0 120"
        args.set_ylim = "-140 0" 
        args.font_size = 12
        args.legend_font_size = 12
        args.no_pk_legend = True
        args.fig_width = 8
        args.marker_size = 10
    elif args.plot_style == 'TF2 paper':
        args.set_xlim = "0 200"
        args.set_ylim = "-200 0" 
        args.font_size = 30 
        args.no_pk_legend = True
        args.fig_width = 12
        args.marker_size = 10

    if args.annotate_points and args.bladders:
        args.fig_width = 30
        args.fig_height = 20
        
    fig.set_figheight(args.fig_height)
    fig.set_figwidth(args.fig_width)
    plt.rcParams.update({'errorbar.capsize':args.errorbar_capsize})
    #plt.rcParams.update({'font.size':30})
    #plt.title(args.title)

    #params = {'legend.fontsize': '20',
    #      'figure.figsize': (1, 1),
    #     'axes.labelsize': '200',
    #     'axes.titlesize':'x-large',
    #     'xtick.labelsize':'x-large',
    #     'ytick.labelsize':'x-large'}
    #matplotlib.rcParams.update(params)
    
    #style = {axes.titlesize : 24,
    #axes.labelsize : 20,
    #xtick.labelsize : 16,
    #ytick.labelsize : 16}

    #matplotlib.rcParams.update(style)
    

    if args.ansys_2d_files is not None or args.ansys_3d_files is not None:
        plot_ansys_data(ax, args)

    if args.ansys_2d_bladder_files is not None:
        plot_ansys_bladders(ax, args)

    paths = args.paths
    plotnames = []

    for i, filepath in enumerate(paths):
        #plt.cla()
        print("filepath", filepath)
        plotname = plot_tf(ax, i, filepath, args)
        plotnames.append(plotname)
    
    lgd, name_suffix = set_ax_parameters(ax, args)

    if args.no_legend: lgd.remove()

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
            print("creating image file", imagename+'.png')
            if args.legend_outside:
                plt.savefig(imagename + '.png', bbox_extra_artists=(lgd,), bbox_inches='tight', numpoints=1, dpi=200)
            else:
                plt.savefig(imagename + '.png', bbox_inches='tight', numpoints=1, dpi=200)



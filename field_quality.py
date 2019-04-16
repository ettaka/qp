#!env python

import numpy as np
import matplotlib.pyplot as plt
import time
#import pandas as pd
import argparse
import codecs
import scipy.interpolate
import scipy.optimize

from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

fig = plt.figure()
ax = fig.add_subplot(111)
ax.autoscale(enable=True, axis='y', tight=True)

def read_value_table(table_name):
    with codecs.open(casedir + '/' + table_name) as f:
        rows = [line.split('\t') for line in f.read().decode('utf-8', 'ignore').replace('\r','').split('\n')]

    row_dicts = {}
    for row in rows:
        if len(row)>3:
            row_dicts[row[0]]={}
            row_dict = row_dicts[row[0]]
            row_dict['name'] = row[0]
            row_dict['unit'] = row[1]
            try:
                row_dict['values'] = np.array([float(value) for value in row[2:]])
            except:
                row_dict['values'] = np.zeros(len(row[2:]))
            row_dict['argmax'] = np.argmax(row_dict['values'])
            row_dict['argmin'] = np.argmin(row_dict['values'])
            row_dict['argmin'] = np.argmin(row_dict['values'])

    return row_dicts
 
def fetch_row_dict(name, row_dicts):
    for row_dict in row_dicts:
        if row_dict['name'] == name:
            return row_dict 
    return None

def add_interp_to_row_dicts(data_dict):
    x_vec = data_dict['data rows']['pos']['values']
    for key in data_dict['data rows']:
        row_dict = data_dict['data rows'][key]
        try:
            row_dict['interp'] = scipy.interpolate.InterpolatedUnivariateSpline(x_vec, row_dict['values'],k=1)
        except:
            row_dict['interp'] = None

def add_integrated_values(data_dicts):
    for key in data_dicts['data rows']:
        row_dict = data_dicts['data rows'][key]
        try:
            row_dict['integrated values'] = np.array([row_dict['interp'].integral(1,i+1)/np.size(row_dict['values']) for i,tmp in enumerate(row_dict['values'])])
        except:
            row_dict['integrated values'] = None

def get_center_integral(pos0, pos_last, magnTFdict):
    def center_integral(pos):
        magnTFinterp = magnTFdict['interp']
        return magnTFinterp.integral(1,pos)-magnTFinterp.integral(1,pos_last)/2.
    return center_integral
    
def find_center_position(data_dict):
    pos_dict = data_dict['data rows']['pos']
    TF_dict = data_dict['data rows']['TF']
    p0 = pos_dict['values'][0]
    p1 = pos_dict['values'][-1]
    center_integral_func = get_center_integral(p0,p1,TF_dict)
    center_position_solution = scipy.optimize.root(center_integral_func, 1)
    center_position = center_position_solution['x'][0]
    tot_integral = TF_dict['interp'].integral(p0,p1) 
    cent_pos_integral = TF_dict['interp'].integral(p0,center_position) 
    print "center position", center_position
    pos_dict['center'] = center_position

    z_dict = data_dict['data rows']['z']
    z_dict['center'] = z_dict['interp'](center_position)
    print "center in meters", z_dict['center']
    print "total integral", tot_integral 
    print "integral until center position", cent_pos_integral
    print "center integral ratio", cent_pos_integral/tot_integral

def add_centered_z(data_dict):
    #z_dict = fetch_row_dict('z', data_dict)
    z_dict = data_dict['data rows']['z']
    z_centered_dict = {}
    data_dict['data rows']['z centered'] = z_centered_dict
    z_centered_dict['name'] = 'z centered'
    z_centered_dict['values'] = z_dict['values'] - z_dict['center']
    z_centered_dict['interp'] = scipy.interpolate.InterpolatedUnivariateSpline(data_dict['data rows']['pos']['values'], z_centered_dict['values'],k=1)

class RowPlotter:

    colors = ['b','g','r','c','m','y','k']
    markers = ['o', '^', 'v', '<', '>', '1', '2', '3', '4', 's', 'p', '*', 'h', '+', 'x']

    def __init__(self, data_dict, args):
        self.i=0
        self.args=args
        self.data_dict = data_dict
        self.xdata = self.data_dict['data rows']['z centered']['values']

    def plot(self, row, plot_integrated_values=False):

        self.i+=1
        if plot_integrated_values: 
            ydata = row['integrated values'] 
            name = row['name'] + ' integ.'
        else:
            ydata = [row['interp'](pos+1) for pos,tmp in enumerate(row['values'])]
            name = row['name']

        save_fig = self.data_dict['casedir'] + '_' + name + '.png'
        ax.plot(self.xdata, ydata,'--'+self.colors[self.i]+self.markers[self.i],label=name)

def plot_data(data_dict, args):
    row_plotter = RowPlotter(data_dict, args)
    row_plotter.plot(data_dict['data rows']['TF'])
    row_plotter.plot(data_dict['data rows']['TF'], True)

    ax.legend(loc=args.legend_location)
    if args.show_plot:
        plt.show()
    else:
        plt.savefig(save_fig)

    return

def read_case(casedir, args):
    #meta not read
    data_dict = {}
    data_dict['casedir'] = casedir
    data_dict['data rows'] = {}
    main_field = read_value_table('main_field')
    normal_multipoles = read_value_table('normal_multipoles')
    skew_multipoles = read_value_table('skew_multipoles')
    data_dict['data rows'].update(main_field)
    data_dict['data rows'].update(normal_multipoles)
    data_dict['data rows'].update(skew_multipoles)

    add_interp_to_row_dicts(data_dict)
    find_center_position(data_dict)
    add_centered_z(data_dict)
    add_integrated_values(data_dict)
    
    plot_data(data_dict, args)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot transfer function')
    parser.add_argument('casedirs', nargs='+', type=str)
    parser.add_argument('-sp', '--show-plot', action='store_true', default=False) 
    parser.add_argument('-ll', '--legend-location', type=str, default='best')

    args = parser.parse_args()
    for casedir in args.casedirs:
        read_case(casedir, args)
    

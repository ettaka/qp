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

def read_value_table(table_name):
    with codecs.open(casedir + '/' + table_name) as f:
        rows = [line.split('\t') for line in f.read().decode('utf-8', 'ignore').replace('\r','').split('\n')]

    row_dicts = []
    for row in rows:
        if len(row)>3:
            row_dicts.append({})
            row_dict = row_dicts[-1]
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

def add_interp_to_row_dicts(row_dicts, x_vec):
    for row_dict in row_dicts:
        try:
            row_dict['interp'] = scipy.interpolate.InterpolatedUnivariateSpline(x_vec, row_dict['values'],k=1)
        except:
            row_dict['interp'] = None

def get_center_integral(pos0, pos_last, magnTFdict):
    def center_integral(pos):
        magnTFinterp = magnTFdict['interp']
        return magnTFinterp.integral(1,pos)-magnTFinterp.integral(1,pos_last)/2.
    return center_integral
    
def find_center_position(main_field):
    pos_dict = fetch_row_dict('pos', main_field)
    TF_dict = fetch_row_dict('TF', main_field)
    p0 = pos_dict['values'][0]
    p1 = pos_dict['values'][-1]
    center_integral_func = get_center_integral(p0,p1,TF_dict)
    center_position_solution = scipy.optimize.root(center_integral_func, 1)
    center_position = center_position_solution['x'][0]
    tot_integral = TF_dict['interp'].integral(p0,p1) 
    cent_pos_integral = TF_dict['interp'].integral(p0,center_position) 
    print "center position", center_position
    pos_dict['center'] = center_position

    z_dict = fetch_row_dict('z', main_field)
    z_dict['center'] = z_dict['interp'](center_position)
    print "center in meters", z_dict['center']
    print "total integral", tot_integral 
    print "integral until center position", cent_pos_integral
    print "center integral ratio", cent_pos_integral/tot_integral


def read_case(casedir, args):
       
    #meta not read
    main_field = read_value_table('main_field')
    normal_multipoles = read_value_table('normal_multipoles')
    skew_multipoles = read_value_table('skew_multipoles')
    add_interp_to_row_dicts(main_field, main_field[0]['values'])
    add_interp_to_row_dicts(normal_multipoles, main_field[0]['values'])
    add_interp_to_row_dicts(skew_multipoles, main_field[0]['values'])
    find_center_position(main_field)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot transfer function')
    parser.add_argument('casedirs', nargs='+', type=str)

    args = parser.parse_args()
    for casedir in args.casedirs:
        read_case(casedir, args)
    

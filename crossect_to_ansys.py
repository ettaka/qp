#!env python

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse
from scipy import interpolate
from scipy import signal

interp1d = interpolate.interp1d
savgol_filter = signal.savgol_filter
#window = 41
polyorder = 1
#npoints = 50
interpkind=1

plt.figure(figsize=(15,10))
mpl.rcParams['lines.markersize'] = 2

def filter_cs(cs, location='mid-planeR', eps=0.001, quadrant=0):
    theta0=quadrant*np.pi/2.
    pole_filter = (45.-41.2)/180.*np.pi
    mp_rlim=113.
    if location=='mid-planeR':
        cs_filtered = cs[(cs['theta']-theta0<eps) & (cs['theta']-theta0>-eps)]
        cs_filtered = cs_filtered[(cs_filtered['r']<mp_rlim)]
    elif location=='mid-planeL':
        cs_filtered = cs[(cs['theta']-theta0-np.pi/2.<eps) & (cs['theta']-theta0-np.pi/2.>-eps)]
        cs_filtered = cs_filtered[(cs_filtered['r']<mp_rlim)]
    elif location=='outerR':
        cs_filtered = cs[(cs['theta']-theta0<np.pi/4.-eps-pole_filter) & (cs['theta']-theta0>eps)]
        cs_filtered = cs_filtered.sort_values('theta')
    elif location=='outerL':
        cs_filtered = cs[(cs['theta']-theta0<np.pi/2.-eps) & (cs['theta']-theta0>np.pi/4.+eps+pole_filter)]
        cs_filtered = cs_filtered.sort_values('theta')
    return cs_filtered

def rotate_vectors(x, z, theta):
    xrot = x * np.cos(theta) - z * np.sin(theta)
    zrot = x * np.sin(theta) + z * np.cos(theta)
    return xrot, zrot

def read_cs(filename, quadrant, args):
    "quadrant is 0, 1, 2 or 3"
    names=['x','y','z','dx','dy','dz','+']
    if args.convert_long:
        cs = pd.read_csv(filename,header=None,names=names,delim_whitespace=True)
    else:
        cs = pd.read_csv(filename,header=None,names=names,delimiter=',')


    if args.convert_long:
        # DPulikowski:
        # Each file consists of 6 columns: x+dx | y+dy | z+dz | -dx | -dy | -dz
        # so we need to correct these
        cs['x'] = cs['x'] + cs['dx']
        cs['y'] = cs['y'] + cs['dy']
        cs['z'] = cs['z'] + cs['dz']
        cs['dx'] = -cs['dx']
        cs['dy'] = -cs['dy']
        cs['dz'] = -cs['dz']
        LE = 607
        CE = 3407
        RE = 7007
        if args.longitudinal_position == 'LE': 
            args.ymin = LE-20
            args.ymax = LE+20
        elif args.longitudinal_position == 'CE': 
            args.ymin = CE-20
            args.ymax = CE+20
        elif args.longitudinal_position == 'RE': 
            args.ymin = RE-20
            args.ymax = RE+20
        elif args.longitudinal_position == 'Station1': 
            args.ymin = 50.
            args.ymax = 1800.
        elif args.longitudinal_position == 'Station2': 
            args.ymin = 2000.
            args.ymax = 4000.
        elif args.longitudinal_position == 'Station3': 
            args.ymin = 4200.
            args.ymax = 6000.
        elif args.longitudinal_position == 'Station4': 
            args.ymin = 6200.
            args.ymax = 7450.
        else:
            args.ymin = float(args.longitudinal_position) - 20.
            args.ymax = float(args.longitudinal_position) + 20.

    cs = cs[(cs['y']>=args.ymin) & (cs['y']<=args.ymax)]

    if not args.convert_long:
        # turn x-axis around
        cs['x']=-cs['x']
        cs['dx']=-cs['dx']

    # on midplane we have some space for shimming
    cs['x'] -= cs['x'].min() - 0.125 
    cs['z'] -= cs['z'].min() - 0.125 

    cs['xmeas'] = cs['x'] + cs['dx']
    cs['zmeas'] = cs['z'] + cs['dz']

    #rotate values to their correct quadrant
    theta = float(quadrant)*np.pi/2.
    cs['x'], cs['z'] = rotate_vectors(cs['x'], cs['z'], theta)

    cs['xmeas'], cs['zmeas'] = rotate_vectors(cs['xmeas'], cs['zmeas'], theta)

    cs['dx'] = cs['xmeas'] - cs['x']
    cs['dz'] = cs['zmeas'] - cs['z']

    cs['r']=np.sqrt(cs['x']**2.+cs['z']**2.)
    cs['theta']=np.arctan2(cs['z'],cs['x'])
    if quadrant>1: cs['theta'] += 2.*np.pi
    cs['theta_deg']=cs['theta']/np.pi*180.

    cs['rmeas']=np.sqrt((cs['x']+cs['dx'])**2.+(cs['z']+cs['dz'])**2.)
    cs['dr']=cs['rmeas']-cs['r']
    cs = cs[cs['r']>80.]

    return cs

def produce_figure(show=False, figname='test.png'):
    if not show:
        plt.savefig(figname, bbox_inches='tight', numpoints=1, dpi=200)
    else:
        plt.show()

class PlotCounter:
    def __init__(self, count=0):
        self.reset(count)

    def next(self):
        self.count += 1
        return self.count
    
    def reset(self, count=0):
        self.count = count

def read_quadrants(filenames, rshim, args):
    coil_dict_list = []
    for quadrant, filename in enumerate(filenames):
        filebase = filename.replace('.txt','')
        with open('../'+filebase+'.size') as f:
            line = f.readline()
            if 'mshim' in line:
                mshim = float(line.split()[1])
                print('Found mshim for', filename, ':', mshim)
        
        coil_dict_list.append({})
        coil_dict = coil_dict_list[-1]
        coil_dict['mshim'] = mshim
        coil_dict['rshim'] = rshim
        coil_dict['cs'] = read_cs(filename, quadrant, args=args)
        cs = coil_dict['cs']
        coil_dict['quadrant'] = quadrant
        coil_dict['filename'] = filename
        coil_dict['filebase'] = filebase
        coil_dict['csf_mpr'] = filter_cs(cs, location='mid-planeR', eps=eps_mp, quadrant=quadrant)
        coil_dict['csf_mpr']['dxshimmed'] = coil_dict['csf_mpr']['dx']
        coil_dict['csf_mpr']['dzshimmed'] = coil_dict['csf_mpr']['dz'] + mshim/2.

        coil_dict['csf_mpl'] = filter_cs(cs, location='mid-planeL', eps=eps_mp, quadrant=quadrant)
        coil_dict['csf_mpl']['dxshimmed'] = coil_dict['csf_mpl']['dx'] + mshim/2.
        coil_dict['csf_mpl']['dzshimmed'] = coil_dict['csf_mpl']['dz']

        # add the midplane shim
        coil_dict['csf_or'] = filter_cs(cs, location='outerR', eps=eps_o, quadrant=quadrant)
        coil_dict['csf_or']['drshimmed'] = coil_dict['csf_or']['dr'] + rshim

        # add the midplane shim
        coil_dict['csf_ol'] = filter_cs(cs, location='outerL', eps=eps_o, quadrant=quadrant)
        coil_dict['csf_ol']['drshimmed'] = coil_dict['csf_ol']['dr'] + rshim

    return coil_dict_list

def set_shimming_and_interpolation(coil_dict_list):
    for quadrant, coil_dict in enumerate(coil_dict_list):
        cs = coil_dict['cs']
        csf_mpr = coil_dict['csf_mpr']
        csf_mpl = coil_dict['csf_mpl']
        csf_or = coil_dict['csf_or'] 
        csf_ol = coil_dict['csf_ol'] 

        # set delta to the normal direction
        if quadrant == 0:
            csf_mpr['dev'] = -csf_mpr['dz']
            csf_mpl['dev'] = -csf_mpl['dx']
            csf_mpr['devshimmed'] = csf_mpr['dzshimmed']
            csf_mpl['devshimmed'] = csf_mpl['dxshimmed']
        elif quadrant == 1:
            csf_mpr['dev'] = csf_mpr['dx']
            csf_mpl['dev'] = -csf_mpl['dz']
            csf_mpr['devshimmed'] = csf_mpr['dzshimmed']
            csf_mpl['devshimmed'] = csf_mpl['dxshimmed']
        elif quadrant == 2:
            csf_mpr['dev'] = csf_mpr['dz']
            csf_mpl['dev'] = csf_mpl['dx']
            csf_mpr['devshimmed'] = csf_mpr['dxshimmed']
            csf_mpl['devshimmed'] = csf_mpl['dzshimmed']
        elif quadrant == 3:
            csf_mpr['dev'] = -csf_mpr['dx']
            csf_mpl['dev'] = csf_mpl['dz']
            csf_mpr['devshimmed'] = csf_mpr['dxshimmed']
            csf_mpl['devshimmed'] = csf_mpl['dzshimmed']

        csf_or['dev'] = csf_or['dr']
        csf_ol['dev'] = csf_ol['dr']

        csf_or['devshimmed'] = csf_or['drshimmed']
        csf_ol['devshimmed'] = csf_ol['drshimmed']

        coil_dict['mpr'] = {}
        coil_dict['mpl'] = {}
        coil_dict['or'] = {}
        coil_dict['ol'] = {}

        coil_dict['mpr']['dev_interp'] = interp1d(csf_mpr['r'], csf_mpr['dev'],kind=interpkind)
        coil_dict['mpl']['dev_interp'] = interp1d(csf_mpl['r'], csf_mpl['dev'],kind=interpkind)
        coil_dict['or']['dev_interp'] = interp1d(csf_or['theta_deg'], csf_or['dev'],kind=interpkind)
        coil_dict['ol']['dev_interp'] = interp1d(csf_ol['theta_deg'], csf_ol['dev'],kind=interpkind)

        coil_dict['mpr']['devshimmed_interp'] = interp1d(csf_mpr['r'], csf_mpr['devshimmed'],kind=interpkind)
        coil_dict['mpl']['devshimmed_interp'] = interp1d(csf_mpl['r'], csf_mpl['devshimmed'],kind=interpkind)
        coil_dict['or']['devshimmed_interp'] = interp1d(csf_or['r'], csf_or['devshimmed'],kind=interpkind)
        coil_dict['ol']['devshimmed_interp'] = interp1d(csf_ol['r'], csf_ol['devshimmed'],kind=interpkind)

def get_mid_plane_dicts(coil_dict_list, window, args):
    mid_plane_dicts = []
    for coil_dict in coil_dict_list:
        mid_plane_dicts.append({})
        mid_plane_dict = mid_plane_dicts[-1]
        cs = coil_dict['cs']
        quadrant = coil_dict['quadrant']
        mid_plane_dict['quadrant'] = quadrant
        quadrant_prev = (quadrant-1)%4

        coil_dict_prev = coil_dict_list[quadrant_prev]
        mid_plane_dict['coil_dict'] = coil_dict
        mid_plane_dict['coil_dict_prev'] = coil_dict_prev
        mid_plane_dict['shim'] = (coil_dict['mshim'] + coil_dict_prev['mshim'])/2.

        csf_mpr = coil_dict['csf_mpr']
        csf_mpl = coil_dict['csf_mpl']
        csf_or = coil_dict['csf_or'] 
        csf_ol = coil_dict['csf_ol'] 

        csf_mpr_prev = coil_dict_prev['csf_mpr']
        csf_mpl_prev = coil_dict_prev['csf_mpl']
        csf_or_prev = coil_dict_prev['csf_or'] 
        csf_ol_prev = coil_dict_prev['csf_ol'] 

        rmax = csf_mpl_prev['r'].max()
        rmin = csf_mpl_prev['r'].min()

        rmax2 = csf_mpr['r'].max()
        rmin2 = csf_mpr['r'].min()

        rmax = min([rmax,rmax2])
        rmin = max([rmin,rmin2])

        mid_plane_dict['rmax'] = rmax
        mid_plane_dict['rmin'] = rmin
        mid_plane_dict['npoints'] = args.npoints

        rspace = np.linspace(rmin, rmax,num=args.npoints)
        mid_plane_dict['rspace'] = rspace
        mpl_mpr_sum = coil_dict_prev['mpl']['dev_interp'](rspace) + coil_dict['mpr']['dev_interp'](rspace)
        mpl_mpr_average = (mpl_mpr_sum)/2.
        mid_plane_dict['mpl_mpr_average'] = mpl_mpr_average
        mid_plane_dict['mpl_mpr_sum'] = mpl_mpr_sum

        mid_plane_dict['mpl_mpr_average_savgol'] = savgol_filter(mpl_mpr_average, window, polyorder)
        mid_plane_dict['mpl_mpr_sum_savgol'] = savgol_filter(mpl_mpr_sum, window, polyorder)

        mid_plane_dict['mpl_interp'] = savgol_filter(coil_dict_prev['mpl']['dev_interp'](rspace), window, polyorder)
        mid_plane_dict['mpr_interp'] = savgol_filter(coil_dict['mpr']['dev_interp'](rspace), window, polyorder)

        
        coil_dict['or']['theta_deg_space'] = np.linspace(csf_or['theta_deg'].min(), csf_or['theta_deg'].max(),num=args.npoints)
        coil_dict['ol']['theta_deg_space'] = np.linspace(csf_ol['theta_deg'].min(), csf_ol['theta_deg'].max(),num=args.npoints)

    return mid_plane_dicts

def plot_mid_planes(mid_plane_dicts):
    for mid_plane_dict in mid_plane_dicts:
        quadrant = mid_plane_dict['quadrant']
        quadrant_prev = (quadrant-1)%4

        coil_dict_prev = mid_plane_dict['coil_dict_prev']
        coil_dict = mid_plane_dict['coil_dict']

        csf_mpr = coil_dict['csf_mpr']
        csf_mpl = coil_dict['csf_mpl']
        csf_or = coil_dict['csf_or'] 
        csf_ol = coil_dict['csf_ol'] 

        csf_mpr_prev = coil_dict_prev['csf_mpr']
        csf_mpl_prev = coil_dict_prev['csf_mpl']
        csf_or_prev = coil_dict_prev['csf_or'] 
        csf_ol_prev = coil_dict_prev['csf_ol'] 

        rmax = mid_plane_dict['rmax']
        rmin = mid_plane_dict['rmin']
        rspace = mid_plane_dict['rspace']
        mpl_mpr_average = mid_plane_dict['mpl_mpr_average']
        mpl_mpr_sum = mid_plane_dict['mpl_mpr_sum']
        mpl_mpr_average_savgol = mid_plane_dict['mpl_mpr_average_savgol']
        mpl_mpr_sum_savgol = mid_plane_dict['mpl_mpr_sum_savgol']
        shim = mid_plane_dict['shim']

        plt.subplot(3,4,plot_counter.next())

        plt.plot(rspace, mid_plane_dict['mpl_interp'], label='__no_legend__', color='C0')
        plt.plot(rspace, mid_plane_dict['mpr_interp'], label='__no_legend__', color='C1')

        plt.plot(csf_mpl_prev['r'], csf_mpl_prev['dev'], label='mpl {}'.format(coil_dict_prev['filebase']), linestyle='None', marker='o')
        plt.plot(csf_mpr['r'], csf_mpr['dev'], label='mpr {}'.format(coil_dict['filebase']), linestyle='None', marker='d')
        plt.plot(rspace, mpl_mpr_average_savgol, label='av filter')
        plt.plot(rspace, mpl_mpr_sum_savgol, label='sum filter')
        plt.plot(rspace, mpl_mpr_average, label='av')
        plt.plot(rspace, mpl_mpr_sum, label='sum')

        plt.ylim(-0.35,0.5)

        if quadrant == 0:
            plt.ylabel('dev (mm)')
        plt.xlabel('r (mm)')
        plt.legend()

def plot_coil_cross_sections(coil_dict_list, args, row_plot=True):

    if args.magnet_name is not None:
        plt.text(0,50, args.magnet_name, horizontalalignment='center', verticalalignment='center')

    for quadrant, coil_dict in enumerate(coil_dict_list):
        if row_plot:
            plt.subplot(3,4,plot_counter.next())

        csf_mpr = coil_dict['csf_mpr']
        csf_mpl = coil_dict['csf_mpl']
        csf_or = coil_dict['csf_or'] 
        csf_ol = coil_dict['csf_ol'] 

        mear_colors = ['red','blue']
        err_scaling = args.error_scaling
        ref_mpl_space_x, ref_mpl_space_z  = rotate_vectors(np.sqrt(2)*0.125, np.sqrt(2)*0.125, quadrant*np.pi/2.)
        datas = [csf_mpr, csf_mpl, csf_or, csf_ol]
        for i,plot_data in enumerate(datas):
            if row_plot:
                xdata, zdata = rotate_vectors(plot_data['x'], plot_data['z'], -quadrant*np.pi/2.)
            else:
                xdata, zdata = plot_data['x']+ref_mpl_space_x*err_scaling, plot_data['z']+ref_mpl_space_z*err_scaling

            name = '__no_legend__'
            if row_plot and i == 0: name = "ref"
            if (not row_plot) and (i == 0 and quadrant == 0): name = "ref"
            plt.plot(xdata, zdata, label=name, color = 'green', zorder=2)

            if not row_plot: plt.text(ref_mpl_space_x*450., ref_mpl_space_z*450., coil_dict['filebase'],horizontalalignment='center', verticalalignment='center')

        for i, plot_data in enumerate(datas):
            if row_plot:
                xdata, zdata = rotate_vectors(plot_data['x'] + plot_data['dx']*err_scaling, plot_data['z'] + plot_data['dz']*err_scaling, -quadrant*np.pi/2.)
            else:
                xdata, zdata = plot_data['x']+ref_mpl_space_x*err_scaling + plot_data['dx']*err_scaling, plot_data['z'] +ref_mpl_space_z*err_scaling+ plot_data['dz']*err_scaling
            name = '__no_legend__'
            if row_plot and i == 0: name = coil_dict['filebase']
            if (not row_plot) and (i == 0 and (quadrant == 0 or quadrant == 1)): name = "meas"
            meas_color = 'red'
            if (not row_plot) and (quadrant%2==1): 
                meas_color = 'blue'
            plt.plot(xdata, zdata, label=name, color=meas_color, zorder=1)

        if quadrant == 0:
            plt.ylabel('z (mm)')
            plt.xlabel('x (mm)')
        plt.legend()

def plot_shimmed_mid_planes(mid_plane_dicts):
    midpd = mid_plane_dicts
    for quadrant,midpd in enumerate(mid_plane_dicts):
        plt.subplot(3,4,plot_counter.next())

        rmax = midpd['rmax']
        rmin = midpd['rmin']
        rspace = midpd['rspace']
        mpl_mpr_average_savgol = midpd['mpl_mpr_average_savgol']
        mpl_mpr_sum_savgol = midpd['mpl_mpr_sum_savgol']
        mpl_mpr_average = midpd['mpl_mpr_average']
        mpl_mpr_sum = midpd['mpl_mpr_sum']
        shim = midpd['shim']

        plt.plot(rspace, mpl_mpr_average_savgol, label='av', color='green')
        plt.plot(rspace, mpl_mpr_sum_savgol+shim, label='sum + shim ({:0.0f} $\mu$m)'.format(1000.*shim,0), color='blue')
        plt.plot(rspace, mpl_mpr_average, label='av filter', color='green')
        plt.plot(rspace, mpl_mpr_sum+shim, label='sum filter + shim ({:0.0f} $\mu$m)'.format(1000.*shim,0), color='blue')

        plt.ylim(-0.15,0.5)

        if quadrant == 0:
            plt.ylabel('dev (mm)')
        plt.xlabel('r (mm)')
        plt.legend()

def plot_outer_radia(coil_dict_list, window):
    csf_o_list = []
    csf_mp_list = []
    for quadrant,coil_dict in enumerate(coil_dict_list):
        cs = coil_dict['cs']

        csf_mpr = coil_dict['csf_mpr']
        csf_mpl = coil_dict['csf_mpl']
        csf_or = coil_dict['csf_or'] 
        csf_ol = coil_dict['csf_ol'] 

        csf_o_list.append(csf_or)
        csf_o_list.append(csf_ol)

        plt.subplot(3,4,plot_counter.next())

        or_interp = coil_dict['or']['dev_interp'] 
        ol_interp = coil_dict['ol']['dev_interp'] 

        theta_deg_space_ol = coil_dict['ol']['theta_deg_space'] 
        theta_deg_space_or = coil_dict['or']['theta_deg_space'] 

        plt.plot(theta_deg_space_or, savgol_filter(or_interp(theta_deg_space_or),window,polyorder), label='__no_legend__ ', color='C0')
        plt.plot(theta_deg_space_ol, savgol_filter(ol_interp(theta_deg_space_ol),window,polyorder), label='__no_legend__ ', color='C1')
        plt.plot(csf_or['theta_deg'], csf_or['dev'], label='or '+coil_dict['filebase'], linestyle='None', marker='o')
        plt.plot(csf_ol['theta_deg'], csf_ol['dev'], label='ol '+coil_dict['filebase'], linestyle='None', marker='d')

        plt.ylim(-0.15,0.15)

        if quadrant == 0:
            plt.ylabel('$\Delta$r (mm)')
        plt.xlabel('theta ($^\circ$)')
        plt.legend()

def plot_shimmed_outer_radia(coil_dict_list, window):
    for quadrant,coil_dict in enumerate(coil_dict_list):
        cs = coil_dict['cs']
        rshim = coil_dict['rshim']

        csf_or = coil_dict['csf_or'] 
        csf_ol = coil_dict['csf_ol'] 

        or_interp = coil_dict['or']['dev_interp'] 
        ol_interp = coil_dict['ol']['dev_interp'] 

        theta_deg_space_ol = coil_dict['ol']['theta_deg_space'] 
        theta_deg_space_or = coil_dict['or']['theta_deg_space'] 

        plt.subplot(3,4,plot_counter.next())

        plt.plot(theta_deg_space_or, savgol_filter(or_interp(theta_deg_space_or),window,polyorder), label='__no_legend__ ', color='green')
        plt.plot(theta_deg_space_ol, savgol_filter(ol_interp(theta_deg_space_ol),window,polyorder), label='__no_legend__ ', color='green')
        plt.plot(theta_deg_space_or, savgol_filter(or_interp(theta_deg_space_or) + rshim,window,polyorder), label='__no_legend__', color='blue')
        plt.plot(theta_deg_space_ol, savgol_filter(ol_interp(theta_deg_space_ol) + rshim,window,polyorder), label='shim {:0.0f} $\mu$m'.format(1000.*rshim), color='blue')


        plt.ylim(-0.25,0.15)

        if quadrant == 0:
            plt.ylabel('$\Delta$r (mm)')
        plt.xlabel('theta ($^\circ$)')
        plt.legend()


def get_savgol_outer_radia(coil_dict_list, window):
    csf_savgol_o_list = []
    for quadrant,coil_dict in enumerate(coil_dict_list):
        cs = coil_dict['cs']
        coil_dict['rshim'] = rshim

        or_interp = coil_dict['or']['dev_interp'] 
        ol_interp = coil_dict['ol']['dev_interp'] 

        theta_deg_space_or = coil_dict['or']['theta_deg_space'] 
        theta_deg_space_ol = coil_dict['ol']['theta_deg_space'] 

        csf_savgol_or_list = {}
        csf_savgol_or_list['theta_deg_space'] = theta_deg_space_or
        csf_savgol_or_list['savgol_interp'] = savgol_filter(or_interp(theta_deg_space_or), window, polyorder) + rshim
        csf_savgol_or = pd.DataFrame(csf_savgol_or_list)

        csf_savgol_ol_list = {}
        csf_savgol_ol_list['theta_deg_space'] = theta_deg_space_ol
        csf_savgol_ol_list['savgol_interp'] = savgol_filter(ol_interp(theta_deg_space_ol), window, polyorder) + rshim
        csf_savgol_ol = pd.DataFrame(csf_savgol_ol_list)

        csf_savgol_o_list.append(csf_savgol_or)
        csf_savgol_o_list.append(csf_savgol_ol)

    csf_savgol_o = pd.concat(csf_savgol_o_list)
    return csf_savgol_o

def write_shimmed_mid_planes(mid_plane_dicts, longposstr):
    for quadrant,midpd in enumerate(mid_plane_dicts):
        mshim = midpd['shim']
        midpd_df_build = {}
        midpd_df_build['rspace_meters'] = midpd['rspace']/1000.
        midpd_df_build['mpl_mpr_average_shimmed_meters'] = (midpd['mpl_mpr_average_savgol'] + mshim)/1000.
        midpd_df_build['mpl_mpr_sum_shimmed_meters'] = (midpd['mpl_mpr_sum_savgol'] + mshim)/1000.
        midpd_df = pd.DataFrame(midpd_df_build)

        profile_name = longposstr + 'mpl_profile'+str(quadrant+1)
        np.savetxt(profile_name + '.dat',midpd_df[['rspace_meters','mpl_mpr_sum_shimmed_meters']], fmt='%a',header='r(m)\tgap(m)')

        f = open('nof_points_' + profile_name + '.dat','w')
        f.write(str(midpd_df.shape[0]))
        f.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Coil crossections to ansys')
    parser.add_argument('paths', nargs='+', type=str)
    parser.add_argument('--rshim', type=float, default=None)
    parser.add_argument('--ymin', type=float, default=-700.)
    parser.add_argument('--ymax', type=float, default=-600.)
    parser.add_argument('--convert-long', action='store_true', default=False)
    parser.add_argument('--error-scaling', type=float, default=100.)
    parser.add_argument('--npoints', type=int, default=50)
    parser.add_argument('-sp', '--show-plot', action='store_true', default=False)
    parser.add_argument('--longitudinal-position', type=str, default="")
    parser.add_argument('--savgol-window', type=int, default=11)
    parser.add_argument('--savgol-window-mpl', type=int, default=11)
    parser.add_argument('--magnet-name', type=str, default=None)

    args = parser.parse_args()
    window = args.savgol_window
    if args.convert_long and args.longitudinal_position == '':
        args.longitudinal_position = 'CE'

    longposstr = args.longitudinal_position
    if args.convert_long:
        longposstr += '_'

    if "Station" in args.longitudinal_position:
        longposstr += 'Station'

    plot_counter = PlotCounter()

    #coil_dict_order = ['204','210','203','212']
    rshim = args.rshim/1000.
    filebase = longposstr + 'l2p2_f_prof_2d_'+str(1000*rshim)
    #rshim -= 0.1 # experiments show that the coils reduce 100um of radial size
    print ("radial shimming:", rshim)
    #filenames = ['204.txt','210.txt','203.txt','212.txt']
    filenames = args.paths
    eps_mp = 0.005
    eps_o = 0.01

    coil_dict_list = read_quadrants(filenames, rshim, args)
    set_shimming_and_interpolation(coil_dict_list)
    mid_plane_dicts = get_mid_plane_dicts(coil_dict_list, args.savgol_window_mpl, args)

    plot_coil_cross_sections(coil_dict_list,args)
    plot_mid_planes(mid_plane_dicts)
    plot_outer_radia(coil_dict_list, window)

    produce_figure(args.show_plot, longposstr + 'coil_sizes_'+str(rshim)+'rshim.png')

    plt.clf()
    plot_counter.reset()

    plot_coil_cross_sections(coil_dict_list, args)
    plot_shimmed_mid_planes(mid_plane_dicts)
    plot_shimmed_outer_radia(coil_dict_list, window)

    produce_figure(args.show_plot, longposstr + 'coil_sizes_shimming_'+str(rshim)+'rshim.png')

    csf_savgol_o = get_savgol_outer_radia(coil_dict_list, window)
    csf_savgol_o['savgol_interp_meters'] = csf_savgol_o['savgol_interp']/1000.
    csf_savgol_o['savgol_savgol_interp'] = savgol_filter(csf_savgol_o['savgol_interp'], window, polyorder*2)
    csf_savgol_o['savgol_savgol_interp_meters'] = csf_savgol_o['savgol_savgol_interp']/1000.
    np.savetxt(filebase+'.dat', csf_savgol_o[['theta_deg_space','savgol_savgol_interp_meters']], fmt='%a',header='theta(deg)\tgap(m)')

    f = open(filebase + 'nof_points.dat','w')
    f.write(str(csf_savgol_o.shape[0]))
    f.close()

    write_shimmed_mid_planes(mid_plane_dicts, longposstr)

    plt.cla()
    plt.subplot(1,1,1)
    plt.plot(csf_savgol_o['theta_deg_space'], csf_savgol_o['savgol_savgol_interp'], label='ol shimmed')
    plt.xlabel('theta ($^\circ$)')
    plt.ylabel('$\Delta$r (mm)')
    plt.legend()
    produce_figure(args.show_plot, longposstr + 'shimmed_theta_deg_'+str(rshim)+'rshim.png')

    plt.clf()
    plt.figure(figsize=(3.5,3.5))
    plt.subplot(1,1,1)
    plt.xlim((-150,150))
    plt.xticks((-150,-100,-50,0,50,100,150))
    plt.ylim((-150,150))
    plt.yticks((-150,-100,-50,0,50,100,150))
    plot_coil_cross_sections(coil_dict_list, args, row_plot=False)
    plt.xlabel('')
    plt.ylabel('')
    magnet_name = ''
    if args.magnet_name is not None:
        magnet_name = args.magnet_name.replace('/','')
    produce_figure(args.show_plot, magnet_name + longposstr + '_crossection.png')

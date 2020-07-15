#!/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import argparse

def read_path(path):
    with open(path,'r') as f:
        data = f.readlines()

    values = {}
    time=0.
    values['fconv']={}
    values['cpval']={}
    values['fconvcpval']={}
    values['fconvval']={}
    values['fconvcrit']={}
    fconv_changed = False
    for line in data:
        if 'Time=' in line:
            try:
                time_new = float(line.split()[1])
            except:
                pass
            if time_new != time: 
                time = time_new
                print("time", time, "detected")
                values['fconv'][time]=[]
                values['fconvval'][time]=[]
                values['fconvcrit'][time]=[]
                values['fconvcpval'][time]=[]
                values['cpval'][time]=[]

        if fconv_changed:
            try:
                values['fconvcpval'][time].append(float(line.split('CP=')[1]))
                fconv_changed=False
            except:
                pass

        if 'F Convergence Norm' in line:
            values['fconv'][time].append(float(line.split()[3]))
        if 'FORCE CONVERGENCE VALUE' in line:
            values['fconvval'][time].append(float(line.split()[4]))
            values['fconvcrit'][time].append(float(line.split()[6]))
            fconv_changed = True
        if 'Iter.' in line and 'CP=' in line and time!=0:
            try:
                values['cpval'][time].append(float(line.split('CP=')[1]))
            except:
                values['cpval'][time].append(0.)

    for time in values['fconv']:
        values['fconv'][time] = np.array(values['fconv'][time])
        values['fconvval'][time] = np.array(values['fconvval'][time])
        values['fconvcrit'][time] = np.array(values['fconvcrit'][time])
        values['fconvcpval'][time] = np.array(values['fconvcpval'][time])
    return values


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='plot ansys convergence')
    parser.add_argument('paths', nargs='+', type=str)
    parser.add_argument('--times', type=float, nargs='+', default=[1.0])
    parser.add_argument('--show-criteria', action='store_true', default=False)
    parser.add_argument('-np','--nof-processes', type=float, nargs='+',default=None)

    plt.yscale('log')

    args = parser.parse_args()
    if args.nof_processes is None:
        args.nof_processes = [1. for i in args.paths]
    if len(args.nof_processes) != len(args.nof_processes):
        print ("The number of processes needs to be given for all paths individually.")
        exit()
    for i,path in enumerate(args.paths):
        print ("path:",path)
        values = read_path(path)
        print ("reading done")
        nofp = args.nof_processes[i]

        for time in args.times:
            if time in values['fconv'].keys():
                name = path.split()[0]
        #        xdata = (values['fconvcpval'][time] - values['fconvcpval'][time][0])/3600./nofp
                ydata = values['fconvval'][time]
                #plt.plot(xdata,ydata, label='{} conv time {}'.format(name,time))
                plt.plot(ydata, label='{} conv time {}'.format(name,time))
                if args.show_criteria:
                    #plt.plot(xdata,values['fconvcrit'][time], label='{} crit time {}'.format(name, time))
                    plt.plot(values['fconvcrit'][time], label='{} crit time {}'.format(name, time))
    plt.legend()
    plt.show()


import os
import pandas as pd
import argparse
from os import listdir
from os.path import isfile, join

def parse_ansys_2d_master_to_master(args):
    parsed_file_data_list = None
    if args.ansys_2d_bladder_files is not None:
        parsed_file_data_list = []
        for path in args.ansys_2d_bladder_files:
            if os.path.isdir(path):
                filename = 'master_to_master_dist.txt'
                name = '_'.join(path.split("QXF_2D")[1].split("-")[1].split("/")[0].split("_")[1:])
                df = pd.read_csv(path+'/'+filename,delim_whitespace=True)
                parsed_file_data_list.append({
                    'DataFrame':df,
                    'name':name,
                    'path':path
                    })
    return parsed_file_data_list

def parse_ansys_2d_files(args):
    """
    possible octants: ['','_2','_3','_4','X','_2X','_3X','_4X']
    """
    parsed_file_data_list = None
    #possible_octants = ['','_2','_3','_4','X','_2X','_3X','_4X']
    possible_octants = ['','_2','_3','_4']
    if args.ansys_show_x:
        possible_octants = ['','_2','_3','_4','X','_2X','_3X','_4X']
    if args.ansys_2d_files is not None:
        parsed_file_data_list = []
        for path in args.ansys_2d_files:
            if os.path.isdir(path):
                filebase = 'MQXF_2d_mech_results_output'
                name = '_'.join(path.split("QXF_2D")[1].split("-")[1].split("/")[0].split("_")[1:])
                testfiles = [f for f in listdir(path) if isfile(join(path, f)) and filebase in f]
                octantfiles = []
                for octant in possible_octants:
                    for testfile in testfiles:
                        if str.lower(filebase+octant+'.txt') in testfile.lower():
                            octantfiles.append((octant,testfile))
                for octant,octantfile in octantfiles:
                    octpath = path+'/' + octantfile
                    print (octpath)
                    octname = '_'.join(octpath.split("QXF_2D")[1].split("-")[1].split("/")[0].split("_")[1:])+'_oct{}'.format(octant)
                    f = open(octpath, 'r')
                    file_dict=dict()
                    for line in f:
                        if 'Step' in line:
                            colnames = line.split()
                            for i, col_name in enumerate(colnames):
                                file_dict[col_name]=[]
                        else:
                            for i,value in enumerate(line.split()):
                                colname = colnames[i]
                                file_dict[colname].append(value)
                    parsed_file_data_list.append({
                        'DataFrame':pd.DataFrame(file_dict),
                        'name':octname,
                        'parent name':name,
                        'path':octpath
                        })
                    df = parsed_file_data_list[-1]['DataFrame']
                    for col in list(df.keys()):
                        try:
                            df[col] = pd.to_numeric(df[col])
                        except:
                            pass

    return parsed_file_data_list

def get_var_in_inp(inpfile, name):

    with open(inpfile) as fd:
        lines = fd.readlines()

    var = None
    for line in lines:
        line = line.split('!')[0]
        line = line.replace(' ','')
        if name+'=' in line:
            try:
                var = line.strip().split(name+'=')[1]
            except:
                pass

    return var

def parse_ansys_3d_files(args):
    parsed_file_data_list = []
    if args.ansys_3d_files is not None:
        for path in args.ansys_3d_files:
            paths = []
            pathroot = str(path)
            path += '/path/'
            file_datas=dict()
            file_datas['filenames']=[]
            file_datas['position']=[]
            file_datas['Step']=[]
            file_datas['datas']=[]
            if os.path.isdir(path): 
                for root, dirs, filenames in os.walk(path+'.'):
                    for fn in filenames:
                        print('filename:',fn)
                        file_datas['filenames'].append(fn)
                        pos = fn.replace('MQXF_3D_mech_path_','').replace('.txt','')
                        Step = pos.split('_')[-1]
                        pos = pos.replace('_'+Step,'')
                        file_datas['position'].append(pos)
                        file_datas['Step'].append(Step)
                        file_datas['datas'].append(pd.read_csv(path+fn,delim_whitespace=True))
                for root, dirs, filenames in os.walk(pathroot+'/rod/.'):
                    for fn in filenames:
                        if not 'force' in fn:
                            print('filename:',fn)
                            file_datas['filenames'].append(fn)
                            pos = 'rod' + fn.replace('MQXF_3D_mech_rod_','').replace('.txt','')
                            Step = pos.split('_')[-1]
                            pos = pos.replace('_'+Step,'')
                            file_datas['position'].append(pos)
                            file_datas['Step'].append(Step)
                            file_datas['datas'].append(pd.read_csv(pathroot+'/rod/'+fn,delim_whitespace=True))
                soluinppath = pathroot+'/inp_bkp/MQXF_solution_3d.inp'
                nomstr = 'if_nom'
                ultstr = 'if_ult'

                with open(soluinppath) as fd:
                    solution_inp = fd.readlines()

                for line in solution_inp:
                    coilcur_found = False
                    if 'Iq(' in line and '=' in line:
                        coilcur_found = True
                        coilcur = [float(cur) for cur in line.split('=')[1].split(',')]
                if not coilcur_found:
                    coilcur = []
                    if_nom = get_var_in_inp(soluinppath, nomstr)
                    if if_nom is None:
                        soluinppath = pathroot+'/inp_bkp/MQXF_parameters.inp'
                        nomstr = 'if_grad'
                        ultstr = 'if_grad_ult'

                    if_nom = get_var_in_inp(soluinppath, nomstr)
                    if if_nom is not None and int(if_nom):
                        coilcur.append(16.47)

                    if_ult = get_var_in_inp(soluinppath, ultstr)
                    if if_ult is not None and int(if_ult):
                        coilcur.append(17.89)

            name = '_'.join(path.split("QXF_3D")[1].split("-")[1].split("/")[0].split("_")[1:])

            file_dict = dict()
            file_dict['Step'] = []
            file_dict['coilcur'] = []
            file_dict['position'] = []
            rod_dict = dict()
            rod_dict['Step'] = []
            rod_dict['coilcur'] = []
            rod_dict['position'] = []
            cols = ['S', 'ZG', 'STH', 'SZ', 'ETH', 'EZ']
            rod_cols = ['S', 'XG', 'YG', 'ZG', 'EZ', 'SZ']
            for col in cols:
                file_dict[col+'_Z0'] = []
            for col in rod_cols:
                rod_dict[col+'_Z0'] = []
            for i, data in enumerate(file_datas['datas']):
                accept = True
                accept_rod = False
                for col in cols:
                    if not col in data: accept = False
                if not accept:
                    accept_rod = True
                    for col in rod_cols:
                        if not col in data: accept_rod = False

                if accept:
                    pos = file_datas['position'][i]
                    Step = file_datas['Step'][i]
                    file_dict['position'].append(pos)
                    file_dict['Step'].append(Step)
                    for col in cols:
                        file_dict[col+'_Z0'].append(data[col].values[0])

                    if 'field' in Step:
                        curstep = int(Step.split('f')[0])-4
                        file_dict['coilcur'].append(float(coilcur[curstep]))
                    else:
                        file_dict['coilcur'].append(0.)
                elif accept_rod:
                    pos = file_datas['position'][i]
                    Step = file_datas['Step'][i]
                    rod_dict['position'].append(pos)
                    rod_dict['Step'].append(Step)
                    for col in rod_cols:
                        rod_dict[col+'_Z0'].append(data[col].values[0])

                    if 'field' in Step:
                        curstep = int(Step.split('f')[0])-4
                        rod_dict['coilcur'].append(coilcur[curstep])
                    else:
                        rod_dict['coilcur'].append(0.)

            df_raw = pd.DataFrame(file_dict)
            df_raw_rod = pd.DataFrame(rod_dict)

            df_shell = df_raw[df_raw['position'].str.contains('shell_15')]
            df_pole = df_raw[df_raw['position'].str.contains('pole')]
            df_rod = df_raw_rod[df_raw_rod['position'].str.contains('rod')]
            df = dict()
            df['Step'] = list(df_pole.sort_values('Step')['Step'].values)
            df['coilcur'] = list(df_pole.sort_values('Step')['coilcur'].values)
            df['scyl'] = list(df_shell.sort_values('Step')['STH_Z0'].values/1e6)
            df['spole'] = list(df_pole.sort_values('Step')['STH_Z0'].values/1e6)
            df['srod'] = list(df_rod.sort_values('Step')['SZ_Z0'].values/1e6)
            df['ecyl'] = list(df_shell.sort_values('Step')['ETH_Z0'].values*1e6)
            df['epole'] = list(df_pole.sort_values('Step')['ETH_Z0'].values*1e6)
            df['erod'] = list(df_rod.sort_values('Step')['EZ_Z0'].values*1e6)
            df = pd.DataFrame(df)
                #'DataFrame':pd.DataFrame(file_dict),

            parsed_file_data_list.append({
                'RawDataFrame':df_raw,
                'DataFrame':df,
                'name':name,
                'path':path
                })

    return parsed_file_data_list

def octname_to_Qname(octname):
    num = octname.split('oct')[1]
    if num == '':
        num = '_1'
    return 'Q'+num[1]

def test3d_parse():
    args = argparse.Namespace(ansys_3d_files=['QXF_3D_20200408-113003_MQXFB_i450_npk_alrod'])
    parsed_file_data_list = parse_ansys_3d_files(args)
    df = parsed_file_data_list[0]['DataFrame']
    print(df)

if __name__ == '__main__':
    test3d_parse()

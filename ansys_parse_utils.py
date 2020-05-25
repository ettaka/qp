import os
import pandas as pd
import argparse
from os import listdir
from os.path import isfile, join

def parse_ansys_2d_files(args):
    """
    possible octants: ['','_2','_3','_4','X','_2X','_3X','_4X']
    """
    parsed_file_data_list = None
    possible_octants = ['','_2','_3','_4','X','_2X','_3X','_4X']
    if args.ansys_2d_files is not None:
        parsed_file_data_list = []
        for path in args.ansys_2d_files:
            if os.path.isdir(path):
                filebase = 'MQXF_2d_mech_results_output'
                testfiles = [f for f in listdir(path) if isfile(join(path, f)) and filebase in f]
                octantfiles = []
                for octant in possible_octants:
                    for testfile in testfiles:
                        if str.lower(filebase+octant+'.txt') in testfile.lower():
                            octantfiles.append((octant,testfile))
                for octant,octantfile in octantfiles:
                    octpath = path+'/' + octantfile
                    print (octpath)
                    name = '_'.join(octpath.split("QXF_2D")[1].split("-")[1].split("/")[0].split("_")[1:])+'_oct{}'.format(octant)
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
                        'name':name,
                        'path':octpath
                        })
                    df = parsed_file_data_list[-1]['DataFrame']
                    for col in list(df.keys()):
                        try:
                            df[col] = pd.to_numeric(df[col])
                        except:
                            pass

    return parsed_file_data_list

def parse_ansys_3d_files(args):
    parsed_file_data_list = []
    if args.ansys_3d_files is not None:
        for path in args.ansys_3d_files:
            paths = []
            path += '/path/'
            file_datas=dict()
            file_datas['filenames']=[]
            file_datas['position']=[]
            file_datas['Step']=[]
            file_datas['datas']=[]
            if os.path.isdir(path): 
                for root, dirs, filenames in os.walk(path+'.'):
                    for fn in filenames:
                        file_datas['filenames'].append(fn)
                        pos = fn.replace('MQXF_3D_mech_path_','').replace('.txt','')
                        Step = pos.split('_')[-1]
                        pos = pos.replace('_'+Step,'')
                        file_datas['position'].append(pos)
                        file_datas['Step'].append(Step)
                        file_datas['datas'].append(pd.read_csv(path+fn,delim_whitespace=True))

            name = '_'.join(path.split("QXF_3D")[1].split("-")[1].split("/")[0].split("_")[1:])

            file_dict = dict()
            file_dict['Step'] = []
            file_dict['position'] = []
            cols = ['S', 'ZG', 'STH', 'SZ', 'ETH', 'EZ']
            for col in cols:
                file_dict[col+'_Z0'] = []
            for i, data in enumerate(file_datas['datas']):
                accept = True
                for col in cols:
                    if not col in data: accept = False

                if accept:
                    pos = file_datas['position'][i]
                    Step = file_datas['Step'][i]
                    file_dict['position'].append(pos)
                    file_dict['Step'].append(Step)
                    for col in cols:
                        file_dict[col+'_Z0'].append(data[col].values[0])

            df_raw = pd.DataFrame(file_dict)

            df_shell = df_raw[df_raw['position'].str.contains('shell_15')]
            df_pole = df_raw[df_raw['position'].str.contains('pole')]
            df = dict()
            df['scyl'] = df_shell.sort_values('Step')['STH_Z0'].values/1e6
            df['spole'] = df_pole.sort_values('Step')['STH_Z0'].values/1e6
            df['Step'] = df_pole.sort_values('Step')['Step'].values
                #'DataFrame':pd.DataFrame(file_dict),

            parsed_file_data_list.append({
                'RawDataFrame':df_raw,
                'DataFrame':df,
                'name':name,
                'path':path
                })

    return parsed_file_data_list

def test3d_parse():
    args = argparse.Namespace(ansys_3d_files=['QXF_3D_20200408-113003_MQXFB_i450_npk_alrod'])
    parsed_file_data_list = parse_ansys_3d_files(args)
    df = parsed_file_data_list[0]['DataFrame']
    print(df)

if __name__ == '__main__':
    test3d_parse()

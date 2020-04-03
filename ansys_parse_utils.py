import os
import pandas as pd

def parse_ansys_2d_files(args):
    parsed_file_data_list = None
    if args.ansys_2d_files is not None:
        parsed_file_data_list = []
        for path in args.ansys_2d_files:
            if os.path.isdir(path): path += '/MQXF_2d_mech_results_output.txt'
            name = '_'.join(path.split("QXF_2D")[1].split("-")[1].split("/")[0].split("_")[1:])
            f = open(path, 'r')
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
                'path':path
                })
            df = parsed_file_data_list[-1]['DataFrame']
            for col in df.keys():
                try:
                    df[col] = pd.to_numeric(df[col])
                except:
                    pass

    return parsed_file_data_list

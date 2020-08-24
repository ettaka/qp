import tdm_loader
import pandas as pd
import numpy as np

SEC_DIFF_1970_1904 = 2082844800.
SEC_DIFF_1970_0001 = 63713606400.
SEC_DIFF_1970_LV = 63763654966.-1596522166.
SEC_DIFF_UTC_GENEVA = 2*24*3600.

def LV_Geneva_timestamp_to_unix(LV_Geneva_timestamp):
    return LV_Geneva_timestamp - SEC_DIFF_1970_LV

to_unix = LV_Geneva_timestamp_to_unix

def load_sm18_tdm(filename, timekey='NTP Time'):
    print ("Loading sm18 TDM file: ", filename)
    data_file = tdm_loader.OpenFile(filename)
    channel_names = []
    channel_dict = data_file.channel_dict(0)
    t = 0
    for i, key in enumerate(channel_dict):
        if not 'time' in key.lower():
            print ("detecting a channel:", key)
            channel_names.append(key)
        else:
            print ("detecting a time channel:", key)
            if t<1: 
                max_time = np.max(channel_dict[key])
                min_time = np.min(channel_dict[key])
                max_nof_rows = np.size(channel_dict[key])
            else:
                max_time_tmp = np.max(channel_dict[key])
                min_time_tmp = np.min(channel_dict[key])
                max_nof_rows_tmp = np.size(channel_dict[key])

                max_time = np.max([max_time, max_time_tmp])
                min_time = np.min([min_time, min_time_tmp])
                max_nof_rows = np.max([max_nof_rows, max_nof_rows_tmp])
            t=t+1

    # time space is not the same for all, let's change that
    timespace = np.linspace(min_time, max_time, max_nof_rows)

    sm18_dict = {}
    sm18_dict[timekey] = to_unix(timespace)
    for chname in channel_names:
        print("Interpolating a channel:", chname)
        interp = np.interp(timespace, channel_dict[chname+'_time'], channel_dict[chname])
        sm18_dict[chname] = interp

    print ("Loading sm18 TDM file done!")

    sm18_df = pd.DataFrame(sm18_dict)

    return sm18_df

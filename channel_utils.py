
def replace_with_dict(string, repdict):
    for key in repdict:
        string = string.replace(key, repdict[key]).replace(key.lower(), repdict[key])
    return string

def create_ss_channel(channel_dict):
    """
    S1CS4_Comp
    S1CS1_Comp
    S1NCS2_Comp
    S1NCS4_Comp
    S1NCS3_Comp
    S1NCS1_Comp
    S1CS3_Comp
    S1CS6_Comp
    S1M5_Comp
    S1M2_Comp
    S1NCS7_Comp
    S1CS7_Comp
    S1CS2_Comp
    S1M3_Comp
    S1M6_Comp
    S1M1_Comp
    S1CS5_Comp
    """
    naming_convention = {
            'S1CS8':'Top Right SC Inside',
            'S1CS7':'Top Right SC Outside',
            'S1CS6':'Top Left SC Inside',
            'S1CS5':'Top Left SC Outside',
            'S1CS4':'Bottom Left SC Inside',
            'S1CS3':'Bottom Left SC Outside',
            'S1CS2':'Bottom Right SC Inside',
            'S1CS1':'Bottom Right SC Outside',
            'S1M8':'Top Right M Inside',
            'S1M7':'Top Right M Outside',
            'S1M6':'Top Left m Inside',
            'S1M5':'Top Left M Outside',
            'S1M4':'Bottom Left M Inside',
            'S1M3':'Bottom Left M Outside',
            'S1M2':'Bottom Right M Inside',
            'S1M1':'Bottom Right M Outside',
            'S1NCS8':'Top Right NCS Inside',
            'S1NCS7':'Top Right NCS Outside',
            'S1NCS6':'Top Left NCS Inside',
            'S1NCS5':'Top Left NCS Outside',
            'S1NCS4':'Bottom Left NCS Inside',
            'S1NCS3':'Bottom Left NCS Outside',
            'S1NCS2':'Bottom Right NCS Inside',
            'S1NCS1':'Bottom Right NCS Outside',
            }
    long_pos_replace = {
            'NCS':'RE',
            'SC':'LE',
            'M':'CE'
            }
    ss_channel = False
    for key in naming_convention:
        if key in channel_dict['name']:
            ss_channel = True
            naming_list = naming_convention[key].split()
            channel_dict['sensor type'] = ''
            channel_dict['location'] = 'Vessel ' + naming_list[0] + ' ' + naming_list[1] + ' '+ naming_list[3]
            channel_dict['unit'] = 'um/m'
            channel_dict['thermal compensator'] = ''
            channel_dict['direction'] = '' 
            channel_dict['physical_quantity'] = 'Strain'
            long_loc = naming_list[2]
            long_loc = replace_with_dict(long_loc, long_pos_replace)
            channel_dict['longitudinal location'] = long_loc
            channel_dict['longname'] = get_channel_longname(channel_dict)
            channel_dict['name'] = get_channel_longname(channel_dict)

    return ss_channel, channel_dict

def create_channel_dict_list(channel_names, channel_units, old_format=False, long_magnet=True):
        channel_dict_list = []
        for i, channel_name in enumerate(channel_names):
                channel_dict_list.append({})
                channel_dict = channel_dict_list[-1]

                channel_dict['signal'] = 'raw'

                #if '_Stress' in channel_name:
                    #channel_name = channel_name.replace('_Stress','')
                    #channel_dict['signal'] = 'post computation'
                if 'Strain_' in channel_name:
                    channel_name = channel_name.replace('Strain_','')
                if '_cryo' in channel_name:
                    channel_name = channel_name.replace('_cryo','')
                if 'FBG_' in channel_name:
                    channel_name = channel_name.replace('FBG_','')+'_FBG'


                if '_inter' in channel_name:
                    channel_dict['interpolation'] = True
                    channel_name = channel_name.split('_interpolation')[0]
                elif '_Inter' in channel_name:
                    channel_dict['interpolation'] = True
                    channel_name = channel_name.split('_Interpolation')[0]

                channel_dict['name'] = channel_name
                is_ss_channel, ss_channel_dict = create_ss_channel(channel_dict)
                if is_ss_channel:
                    channel_dict = ss_channel_dict
                else:
                    channel_dict['sensor type'] = ''
                    channel_dict['location'] = ''
                    channel_dict['unit'] = channel_units[i]
                    channel_dict['thermal compensator'] = ''
                    channel_dict['direction'] = '' 
                    channel_dict['physical_quantity'] = ''

                    channel_dict['longitudinal location'] = ''
                    if '_LE' in channel_name or 'LE_' in channel_name:
                        channel_dict['longitudinal location'] = 'LE'
                    elif '_CE' in channel_name or '_MI' in channel_name or 'CE_' in channel_name or 'MI_' in channel_name:
                        channel_dict['longitudinal location'] = 'CE'
                    elif '_RE' in channel_name or 'RE_' in channel_name:
                        channel_dict['longitudinal location'] = 'RE'

                    if old_format and '_comp' in channel_name.lower(): 
                        print("Using old format!")
                        if '_compensated' in channel_name.lower(): 
                            channel_dict['thermal compensator'] = 'compensated'
                        else:
                            channel_dict['thermal compensator'] = 'compensator'

                    if 'SG_' == channel_name[0:3]:
                        channel_name = channel_name.replace('SG_','')
                        channel_dict['sensor type'] = 'Resistive'
                        channel_dict['sensor type'] = ''

                    if channel_name.lower() == 'current':
                        channel_dict['name'] = 'Current'
                        channel_dict['physical_quantity'] = 'Current'
                    elif channel_name == 'Time':
                        channel_dict['physical_quantity'] = 'Time'
                    elif 'NTP' in channel_name and 'TIME' in channel_name:
                        channel_dict['physical_quantity'] = 'Time'
                        channel_dict['sensor type'] = 'NTP'
                        channel_dict['longname'] = 'NTP Time'
                        channel_dict['name'] = 'NTP Time'
                    elif 'pressure' in channel_name.lower() or ('bladder' in channel_name.lower()) or 'rods_p' in channel_name.lower() or 'piston' in channel_name.lower():
                        channel_dict['unit'] = 'Bar'
                        channel_dict['physical_quantity'] = 'Pressure'
                        if 'rod' in channel_name.lower() or 'piston' in channel_name.lower():
                            channel_dict['name'] = 'Piston Pressure'
                            channel_dict['longname'] = 'Piston Pressure'
                        elif 'bladder' in channel_name.lower():
                            channel_dict['name'] = 'Bladder Pressure'
                            channel_dict['longname'] = 'Bladder Pressure'
                    else:
                        if channel_name[0:2] == 'SH':
                                shell_location = ''
                                if channel_name[2] == 'R': shell_location = 'Right'
                                elif channel_name[2] == 'L': shell_location = 'Left'
                                elif channel_name[2] == 'T': shell_location = 'Top'
                                elif channel_name[2] == 'B': shell_location = 'Bottom'
                                if channel_name[3] == 'Z': channel_dict['direction'] = 'Longitudinal'
                                elif channel_name[3] == 'T': channel_dict['direction'] = 'Azimuthal'
                                channel_dict['location'] = 'Shell ' + shell_location
                                channel_dict['physical_quantity'] = 'Strain'
                                channel_dict['material'] = 'aluminium'
                        elif channel_name[0] == 'R' and 'n' not in channel_name: 
                                channel_dict['location'] = 'Rod ' + channel_name[1]
                                channel_dict['physical_quantity'] = 'Strain'
                                channel_dict['material'] = 'aluminium'
                                channel_dict['direction'] =  'Longitudinal'
                                channel_dict['rod diameter'] = 0.036
                                if long_magnet: 
                                    print ("Long magnets have stainless steel rods, setting material for:", channel_name)
                                    channel_dict['material'] = 'stainless steel'
                                    channel_dict['rod diameter'] = 0.035
                        elif channel_name[0:2] == 'CO':
                                channel_dict['location'] = 'Coil ' + channel_name[2:5]
                                channel_dict['physical_quantity'] = 'Strain'
                                channel_dict['material'] = 'titanium'
                        elif channel_name == 'Current':
                                channel_dict['physical_quantity'] = 'Current'
                                channel_dict['location'] = 'None'
                        else:
                                channel_dict['location'] = channel_name 
                                channel_dict['physical_quantity'] = 'None'

                        if channel_dict['sensor type'] == '':
                            if '_SG' in channel_name in channel_name:
                                    channel_dict['sensor type'] = 'Resistive'
                                    channel_dict['sensor type'] = ''
                            elif '_FBG' in channel_name:
                                    channel_dict['sensor type'] = 'Optical'
                                    channel_name = channel_name.replace('_FBG', '')
                            else: channel_dict['sensor type'] = ''
                        if '_Z' in channel_name or channel_name[-1] == 'Z' or channel_name[-2:] == 'Z_':
                            channel_dict['direction'] =  'Longitudinal'
                        elif '_T' in channel_name or channel_name[-1] == 'T' or channel_name[-2:] == 'T_':
                            channel_dict['direction'] = 'Azimuthal'

                        if channel_dict['signal'] == 'post computation':
                            channel_dict['physical_quantity'] = 'Stress'
                            channel_dict['unit'] = 'MPa'

                        channel_dict['longname'] = get_channel_longname(channel_dict)
                        channel_dict['name'] = get_channel_longname(channel_dict)

        return channel_dict_list


def get_channel_longname(channel_dict):
    longname = channel_dict['location']
    if not channel_dict['direction'] == '': longname += ' ' + channel_dict['direction']
    if not channel_dict['physical_quantity'] == '': longname += ' ' + channel_dict['physical_quantity']
    if not channel_dict['sensor type'] == '': longname += ' ' + channel_dict['sensor type']
    if not channel_dict['thermal compensator'] == '': longname += ' ' + channel_dict['thermal compensator']
    if not channel_dict['longitudinal location'] == '': longname += ' ' + channel_dict['longitudinal location']
    if not channel_dict['signal'] == '' and channel_dict['signal'] == 'post computation': longname += ' ' + channel_dict['signal']

    return longname

def get_channel_name(channel_dict):
    return channel_dict['location'] + ' ' + channel_dict['direction'] + ' ' + channel_dict['physical_quantity']

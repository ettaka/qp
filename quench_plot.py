import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd
import cPickle as pickle
import os
import codecs
import argparse
import tabulate

"""
 %% Compute Stress Option 
 E_ti = 130e6/10e5/10e2 % GPa,muStrain,MPa - Warm = Cold
 v_ti = 0.30
 k_plane_ti = E_ti/(1-v_ti^2)
 E_al = 79e6/10e5/10e2
 % Exp. Values
 Nq = length(DB.Current)
 dmax = 13
    for iq = 1 : Nq
      DB.Stress.Coil.C105.Theta{iq} = k_plane_ti*(DB.Coil.C105.Theta{iq}+v_ti*DB.Coil.C105.Zeta{iq})
      DB.Stress.Coil.C106.Theta{iq} = k_plane_ti*(DB.Coil.C106.Theta{iq}+v_ti*DB.Coil.C106.Zeta{iq})
      DB.Stress.Coil.C107.Theta{iq} = k_plane_ti*(DB.Coil.C107.Theta{iq}+v_ti*DB.Coil.C107.Zeta{iq})
      DB.Stress.Coil.C007.Theta{iq} = k_plane_ti*(DB.Coil.C007.Theta{iq}+v_ti*DB.Coil.C007.Zeta{iq})
  
"""
MATERIAL_DATA = { 
		'titanium' : {
			'elastic_modulus':'130e9',
			'poisson_ratio':'0.30'
			},
		'aluminium' : {
			'elastic_modulus':'79e6',
			'poisson_ratio':'0.34'
			}
		}

def compute_stress(strain_theta, strain_z, material, axis):
	"""
	stress_z     = E/(1-v**2)*(strain_z     * v * strain_theta)
	stress_theta = E/(1-v**2)*(strain_theta * v * strain_z)
	"""
	E = MATERIAL_DATA[material]['elastic_modulus']
	v = MATERIAL_DATA[material]['poisson_ratio']
	if axis=='theta': 
		strain_1 = strain_theta
		strain_2 = strain_z
	elif axis=='z': 
		strain_1 = strain_z
		strain_2 = strain_theta
	return E/(1.-v**2.)*(strain_1 + v * strain_2)

def timeit(func):
	def wrapper(*args, **kwargs):
		time_start = time.time()
		o=func(*args, **kwargs)
		print 'timeit: ', func.__name__, ':{0:0.2f}s'.format(time.time()-time_start)
		return o
	return wrapper

def tol_cols(array):
	return np.array([np.max(col)-np.min(col) for col in array.transpose()])

@timeit
def detect_max_current(mtbop_data):
	raw_data = mtbop_data['raw_data']
	current_array = raw_data[:, find_channel_index(mtbop_data['channel_names'], 'Current')]
	grad_array = np.gradient(current_array)
	qinx = np.argmin(grad_array)
	qinx = np.argmax(current_array)
	mtbop_data['max_current_index'] = qinx
	mtbop_data['max_current'] = current_array[qinx]
	mtbop_data['max_current_grad'] = grad_array[qinx]
	print "mtbop_data['max_current_index']", mtbop_data['max_current_index']
	print "mtbop_data['max_current_grad']", mtbop_data['max_current_grad']
	print "mtbop_data['max_current']", mtbop_data['max_current']

def detect_first_values(mtbop_data):
	raw_data = mtbop_data['raw_data']
	filtered_array = mtbop_data['filtered_array']
	initial_values = None
	for i, row in enumerate(raw_data):
		if np.all(filtered_array[i]==0): 
			initial_values = raw_data[i]
			break
	mtbop_data['initial_values'] = initial_values
	

@timeit
def dist_filter(array, abs_dist, filtered_array, turn_on=False):
	"""
	Filter Numpy array column wise. Collect first point of the array and then every point for which the 
	column wise value deviates a distance from the last collected point.
	
	Parameters
	__________
	array : numpy.ndarray
		data array to be filtered
	abs_tol : numpy.ndarray
		column wise distances for the filtering

	Returns
	-------
	filtered_array : numpy.ndarray (bool)
		filtered out elements
	"""
	abs_dist = np.array(abs_dist)
	filtered_array_new = []
	last_mark=None
	for i, row in enumerate(array):
		#if np.any(filtered_array[i] > 0):
			#filtered_row[:] = True
		if last_mark is None or not turn_on:
			last_mark=row
			filtered_row = np.array(np.absolute(row) != np.absolute(row))
		else:
			filtered_row = np.absolute(row-last_mark) < abs_dist
			last_mark = filtered_row*last_mark+(np.logical_not(filtered_row))*row

		filtered_array_new.append(filtered_row)
	return np.array(filtered_array_new)

def filter_indices_after_index(mtbop_data):
	filtered_array = mtbop_data['filtered_array']
	index = mtbop_data['max_current_index']
	filtered_array[index:,:] = 1
	
def index_union(filtered_array, selected_axes_dict):
	inx_sum = np.sum(filtered_array[:,selected_axes_dict], axis=1)
	inx_union = np.where(inx_sum == 0)[0]
	return inx_union

def index_cut(filtered_array, selected_axes_dict):
	inx_product = np.prod(filtered_array[:,selected_axes_dict], axis=1)
	inx_cut = np.where(inx_product == 0)[0]
	return inx_cut

@timeit
def select_mtbop_data(mtbop_data, selected_axes):
	filtered_array = mtbop_data['filtered_array']
	initial_values = mtbop_data['initial_values']
	raw_data = mtbop_data['raw_data']
	inxs = index_union(filtered_array, selected_axes)
	selected_data = raw_data[:,selected_axes][inxs,:]
	if mtbop_data['plot_delta']:
		selected_data[:,:] -= initial_values[[selected_axes]]
	return selected_data

@timeit
def plot_selected(mtbop_data):
	selected_axes_dict = mtbop_data['selected_axes_dict']
	raw_data = mtbop_data['raw_data']
	filtered_array = mtbop_data['filtered_array']
	max_current = mtbop_data['max_current']
	for axes in selected_axes_dict:
		plot_data = select_mtbop_data(mtbop_data, axes['selected_axes'])
		x = plot_data[:,0]
		y = plot_data[:,1]
		plt.xlabel(axes['0 name'] + '(' + axes['0 unit'] + ')')
		plt.ylabel(axes['1 name'] + '(' + axes['1 unit'] + ')')
		if axes['0 scaling'] == 'current': 
			print "Using current scaling for axis 1"
			x = (x/max_current)**2. 
			plt.xlabel('$(I/I_q)^2$')
		if axes['1 scaling'] == 'current': 
			print "Using current scaling for axis 2"
			y = (y/max_current)**2. 
			plt.ylabel('$(I/I_q)**2$')
		plt.plot(x,y)

@timeit
def plot_selected_list(mtbop_data_list):
	for mtbop_data in mtbop_data_list:
		plot_selected(mtbop_data)

def plot_raw_current(mtbop_data):
	raw_data = mtbop_data['raw_data']
	current_index = mtbop_data['current_index']
	plt.plot(raw_data[:, current_index], label=mtbop_data['filepath'])

@timeit
def plot_selected_list_raw_current(mtbop_data_list):
	for mtbop_data in mtbop_data_list:
		plot_raw_current(mtbop_data)
	plt.legend()

@timeit
def load_mtbop_data(filepath, override_pickle=False, header_lines=38):
	"""
	Loads mtbop data from txt file.
	
	Parameters
	__________
	filepath : string
		path of the mtbop data file
	header_lines : int
		number of header lines in mtbop file
	
	Returns
	-------
	mtbop_data : json?
		{
		 'filepath' = string,
		 channel_names = [string, .., string],
		 channel_units = [string, .., string],
		 raw_data = numpy.ndarray}
	"""
	basename = os.path.splitext(filepath)[0]
	pickle_path = basename + '.pickle'
	if os.path.isfile(pickle_path) and not override_pickle:
		print "Found a pickle:", pickle_path, "for", filepath, "using that instead."
		return cpickle_load(pickle_path)

	mtbop_data = {}
	with codecs.open(basename + '.txt') as f:
		head = [next(f).decode('utf-8', 'ignore') for x in xrange(header_lines)]
	channel_names = head[8].strip('\r\n').split('\t')
	channel_names = [name.split()[0].replace('_cryo','') for name in channel_names]
	channel_units = head[9].strip('\r\n').split('\t')
	#raw_data_df = pd.read_csv(filepath, sep = '\t', skiprows=header_lines, names=channel_names)

	with codecs.open(filepath) as f:
		raw_data = np.loadtxt(f, skiprows=header_lines)
	mtbop_data['filepath'] = filepath
	mtbop_data['channel_names'] = channel_names
	mtbop_data['channel_units'] = channel_units
	#mtbop_data['raw_data_df'] = raw_data_df
	#mtbop_data['raw_data'] = raw_data_df.values
	mtbop_data['raw_data'] = raw_data
	mtbop_data['filtered_array'] = np.zeros(np.shape(raw_data))
	mtbop_data['channel_tolerance'] = tol_cols(mtbop_data['raw_data'])
	detect_max_current(mtbop_data)

	filter_indices_after_index(mtbop_data)
	detect_first_values(mtbop_data)
	
	cpickle_dump(mtbop_data, pickle_path)
	return mtbop_data

def load_mtbop_list(filepath_list, override_pickle=False, header_lines=38):
	mtbop_data_list = []
	for filepath in filepath_list:
		mtbop_data = load_mtbop_data(filepath, override_pickle)
		mtbop_data_list.append(mtbop_data)
	return mtbop_data_list

def find_channel_index(channel_names, channel_name):
	return channel_names.index(channel_name)

def find_axes(mtbop_data):
	channel_names = mtbop_data['channel_names']
	name_pairs = mtbop_data['selected_channel_names']
	channel_units = mtbop_data['channel_units']
	mtbop_data['current_index'] = channel_names.index('Current')
	axes_list = []
	for pair in name_pairs:
		if pair[0] in channel_names and pair[1] in channel_names: 
			first = channel_names.index(pair[0])
			second = channel_names.index(pair[1])
			axes = {'selected_axes': (first,second)}
			if pair[0].lower() == 'current':
				axes['0 scaling']='current'
			else: 
				axes['0 scaling']=''
			if pair[1].lower() == 'current':
				axes['1 scaling']='current'
			else: 
				axes['1 scaling']=''
			axes['0 name'] = pair[0]
			axes['1 name'] = pair[1]
			axes['0 unit'] = channel_units[first]
			axes['1 unit'] = channel_units[second]
			axes_list.append(axes)
	return axes_list

@timeit
def test_pandas_load_data(filepath, header_lines=38):
	raw_data = pd.read_csv(filepath, sep = '\t', skiprows=header_lines,usecols=range(21))
	print raw_data.values[2]

@timeit
def test_numpy_load_data(filepath, header_lines=38):
	raw_data = np.loadtxt(filepath, skiprows=header_lines)
	print raw_data[2]

def test_load_data():
	filepath = 'Data/test/01QuenchDC_161013_16236.txt'
	test_pandas_load_data(filepath)
	test_numpy_load_data(filepath)

@timeit
def test_cpickle_dump():
	filepath = 'Data/test/01QuenchDC_161013_16236.txt'
	mtbop_data = load_mtbop_data(filepath)
	with open('test.pickle', 'wb') as handle:
		pickle.dump(mtbop_data, handle, protocol=pickle.HIGHEST_PROTOCOL)
@timeit
def cpickle_dump(data, pickle_path):
	with open(pickle_path, 'wb') as handle:
		pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)
	
@timeit
def cpickle_load(pickle_path):
	with open(pickle_path, 'rb') as handle:
		mtbop_data = pickle.load(handle)
	return mtbop_data

def print_max_current(mtbop_data):
	print mtbop_data['filepath'], mtbop_data['max_current']

def make_info_table(mtbop_data_list):
	info_table = [['File path', 'Max Current']]
	for mtbop_data in mtbop_data_list:
		info_table.append([mtbop_data['filepath'], mtbop_data['max_current']])
	return info_table

def test_plot_selected():
	#use_dist_filter = False
	#nof_points = 1e4
	filepath = 'Data/test/01QuenchDC_161013_16236.txt'
	#filepath = 'Data/test/02QuenchDC_161013_18055.txt'
	selected_channel_names = [('Current','CO105Z_cryo'), ('Current','CO105T_cryo')]
	#selected_channel_names = [('Current','CO105T_cryo')]
	plot_delta = True
	#selected_channel_names = [('CO105T_cryo','Current')]
	#selected_channel_names = [('Current','SHTT_cryo')]
	#selected_channel_names = [('CO105T_cryo','SHTT_cryo')]

	mtbop_data = load_mtbop_data(filepath)
	mtbop_data['selected_channel_names'] = selected_channel_names
	mtbop_data['selected_axes_dict'] = find_axes(mtbop_data)
	mtbop_data['plot_delta'] = plot_delta

	#mtbop_data['filter_approx_points'] = nof_points
	#filtered_array = dist_filter(raw_data, tol/nof_points, filtered_array, turn_on=use_dist_filter)

	plot_selected(mtbop_data)

	plt.show()

def test_plot_selected_list():
	#datadir = 'Data/MQXFS4' 
	#datadir = 'Data/HCMQXFS004-CR000001' 
	#filepath_list = [datadir + '/' + f for f in os.listdir(datadir) if ".txt" in f]
	max_current_low_limit = 15
	
	#filepath_list = ['Data/test/01QuenchDC_161013_16236.txt', 
		         #'Data/test/02QuenchDC_161013_18055.txt']
	#selected_channel_names = [('CO105T_cryo','Current')]

	filepath_list = ['Data/HCMQXFS004-CR000001/180704_12_15_23_Quench.txt', 
                         'Data/HCMQXFS004-CR000001/180709_14_40_50_Quench.txt', 
                         'Data/HCMQXFS004-CR000001/180730_18_07_11_Quench.txt']

	print filepath_list

	mtbop_data_list = []
	selected_channel_names = [('Current','CO109T')]
	plot_delta = True
	for filepath in filepath_list:
		mtbop_data = load_mtbop_data(filepath)
		if mtbop_data['max_current'] > max_current_low_limit:
			mtbop_data_list.append(mtbop_data)
			mtbop_data = mtbop_data_list[-1]
			mtbop_data['selected_channel_names'] = selected_channel_names
			mtbop_data['selected_axes_dict'] = find_axes(mtbop_data)
			mtbop_data['plot_delta'] = plot_delta
			print_max_current(mtbop_data)

	plot_selected_list(mtbop_data_list)
	#plot_selected_list_raw_current(mtbop_data_list)

	plt.show()

def plot_selected_files(filepath_list, selected_channel_names):
	max_current_low_limit = 15

	mtbop_data_list = []
	plot_delta = True
	for filepath in filepath_list:
		mtbop_data = load_mtbop_data(filepath)
		if mtbop_data['max_current'] > max_current_low_limit:
			mtbop_data_list.append(mtbop_data)
			mtbop_data = mtbop_data_list[-1]
			mtbop_data['selected_channel_names'] = selected_channel_names
			mtbop_data['selected_axes_dict'] = find_axes(mtbop_data)
			mtbop_data['plot_delta'] = plot_delta
			print_max_current(mtbop_data)

	plot_selected_list(mtbop_data_list)
	#plot_selected_list_raw_current(mtbop_data_list)

	plt.show()

def testing():
	test_plot_selected()
	#test_plot_selected_list()
	#test_load_data()
	#test_cpickle_dump()
	#cpickle_load('test.pickle')


def print_mtbop_info(filepaths):
	info_table = make_info_table(mtbop_data_list)
	print tabulate.tabulate(info_table, headers="firstrow")

def find_channel_names(filepaths):
	channel_names_list = []
	for mtbop_data in mtbop_data_list:
		channel_names_list += mtbop_data['channel_names']
	return list(set(channel_names_list))

#def find_strain_pairs(filepaths):
#	channel_names_list = find_channel_names(filepaths)
#	strain_pairs_list = []
#	for channel_name in channel_names_list:
#		strain_pairs_list.append((channel_name, difflib.get_close_matches(channel_name, channel_names_list)))
#	return strain_pairs_list

def print_channel_names(filepaths):
	print "Channel names:"
	print find_channel_names(filepaths)
	#print find_strain_pairs(filepaths)
	
def find_files(paths):
	filepaths=[]
	for path in paths:
		if os.path.isdir(path):
			filepaths += [path + '/' + f for f in os.listdir(path) if ".txt" in f]
		elif os.path.isfile(path) and ".txt" in f:
			filepaths.append(path)
	return filepaths

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Plot mtbop strain data')
	parser.add_argument('paths', nargs='+', type=str)
	parser.add_argument('-i', '--info', action='store_true', default=False)
	parser.add_argument('-c', '--channel-names', action='store_true', default=False)
	parser.add_argument('-p', '--plot', action='store_true', default=False) 
	parser.add_argument('-x', '--x-axis', type=str)
	parser.add_argument('-y', '--y-axis', nargs='+', type=str)
	parser.add_argument('-s', '--plot-stress', nargs='+', type=str)
	parser.add_argument('-o', '--override-pickle', action='store_true', default=False) 

	args = parser.parse_args()
	paths = args.paths
	info = args.info
	channel_names = args.channel_names
	plot = args.plot
	x_axis = args.x_axis
	y_axis = args.y_axis
	plot_stress = args.plot_stress
	override_pickle = args.override_pickle

	filepaths = find_files(paths)

	mtbop_data_list = load_mtbop_list(filepaths, override_pickle)

	if info: print_mtbop_info(filepaths)
	if channel_names: print_channel_names(filepaths)
	if plot:
		selected_channel_names = [(x_axis,y) for y in y_axis]
		plot_selected_files(filepaths, selected_channel_names, plot_stress)

	




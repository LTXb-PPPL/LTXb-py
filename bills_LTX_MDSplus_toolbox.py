# Copied from \u\phughes\python\pauls_LTX_MDSplus_toolbox.py on 28Oct2019

from __future__ import print_function
from math import *
# from scipy import *
# from scipy.integrate import cumtrapz
from scipy.interpolate.interpolate import interp1d
from pylab import *
import time as timetools
import MDSplus as mds
import inspect, os

diagnosticsyes = 0  # general flag for diagnostic reporting

''' ### REMOTE vs IN-LAB differentiation based on run path happens here ### '''
pathdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
# print( os.path.exists('/u/phughes/python'),pathdir )
if os.path.exists('/u/wcapecch/python'):
	runfromhome = 0
	print("In the lab eh? Carry on")
else:
	runfromhome = 1
	print(f"how'd I get here?")
	# print("I see you're working from home... enjoy your coffee!!")

# elif 'Paul\\Documents\\Laboratory' in pathdir:

# ************************************ LTX_B CONSTANTS **************************************

nodebase_magnetic = '.DIAGNOSTICS.MAGNETIC'
nodebase_mirnovs = nodebase_magnetic + '.MIRNOVS.'
nodebase_SC = nodebase_magnetic + '.SADDLE_COILS:'
nodes_FL_C = ['\FL_CS_{0:02}'.format(flstep) for flstep in range(1, 11 + 1)]
nodes_FL_U = ['\FL_US_{0}'.format(flstep) for flstep in range(1, 8 + 1)]
nodes_FL_L = ['\FL_LS_{0}'.format(flstep) for flstep in range(1, 8 + 1)]
nodes_TA = ['\MIRNOV_TA_S{0:02}'.format(tastep) for tastep in range(1, 10 + 1)]
TAdict = {'A': 1, 'D': 2, 'E': 3, 'F': 4, 'H': 5, 'I': 6, 'L': 7, 'M': 8, 'N': 9, 'P': 10}
SESAsen = ['I_EDGE', 'JK_UPPER', 'JK_MIDPLANE', 'LM_MIDPLANE', 'NO_MIDPLANE', 'NO_UPPER', 'P_EDGE']
SESAdir = ['PERP', 'TAN']
nodes_SESA = ['\MIRNOV_{0}_{1}'.format(SESAsen[sstep], SESAdir[dstep]) for dstep in range(0, 2) for sstep in
              range(0, 7)]
nodes_REPA_a = ['\MIRNOV_RE_PA{0:02}A'.format(pastep) for pastep in range(1, 20 + 1)]
nodes_REPA_p = ['\MIRNOV_RE_PA{0:02}P'.format(pastep) for pastep in range(1, 20 + 1)]
nodes_REPA_t = ['\MIRNOV_RE_PA{0:02}T'.format(pastep) for pastep in range(1, 20 + 1)]
nodes_REPA = [nodes_REPA_a, nodes_REPA_p, nodes_REPA_t]
# node_IP='\IP'
node_IP = '.DIAGNOSTICS.MAGNETIC:IP_ROG'
node_Vloop = nodebase_magnetic + ':V_LOOP'
node_DML = nodebase_magnetic + ':DIAMAG_LOOP'
nodes_I_coil = ['\I_COIL_' + thiscoil for thiscoil in
                ['TF', 'YELLOW', 'GREEN', 'ORANGE_T', 'ORANGE_B', 'BLUE', 'RED', 'INTERNAL', 'OH']]
dict_nodes_I_coil = {'tf': '\I_COIL_TF',
                     'yellow': '\I_COIL_YELLOW',
                     'green': '\I_COIL_GREEN',
                     'orange_t': '\I_COIL_ORANGE_T',
                     'orange_b': '\I_COIL_ORANGE_B',
                     'blue': '\I_COIL_BLUE',
                     'red': '\I_COIL_RED',
                     'internal': '\I_COIL_INTERNAL',
                     'oh': '\I_COIL_OH'}
nodes_RTD_BD = ['.OPER_DIAGS.RTDS.BEAM_DUMP:' + 'SENSOR_{0:02}'.format(RTDstep) for RTDstep in range(1, 11 + 1)]


# ********************************** FUNCTIONS FROM HERE ON DOWN *************************************

def get_tree_conn(temp_shot, treename='ltx_b'):
	"""
	Establishes a tree pointer to pass to commands like "get_data."  Automatically differentiates between remote or in-lab access based on the above run path criterion.
	Arguments:
		<temp_shot>: shot number target for the tree pointer (int)
	Returns:
		<temp_conn>: tree pointer (class 'MDSplus.tree.Tree')
	"""
	if runfromhome:
		try:
			temp_conn = mds.Connection('lithos')
		except:
			temp_conn = mds.Connection('lithos.pppl.gov:8000::')  # "lithos.pppl.gov:8000::" in .cshrc file
		temp_conn.openTree(treename, temp_shot)
	else:
		# print(treename, temp_shot)
		temp_conn = mds.Tree(treename, temp_shot)
	return temp_conn


def get_transp_tree_conn(temp_shot, treename='transp_ltx'):
	"""
	Establishes a tree pointer to pass to commands like "get_data."  Automatically differentiates between remote or in-lab access based on the above run path criterion.
	Arguments:
		<temp_shot>: shot number target for the tree pointer (int)
	Returns:
		<temp_conn>: tree pointer (class 'MDSplus.tree.Tree')
	"""
	if runfromhome:
		temp_conn = mds.Connection('transpgrid')
		temp_conn.openTree(treename, temp_shot)
	else:
		temp_conn = mds.Tree(treename, temp_shot)
	return temp_conn


def get_data(mytree, node, t_start=None, t_end=None, times=None):

	try:
		if runfromhome:
			rawdata = mytree.get(node).data()
		else:
			mydatarec = mytree.getNode(node).record
			rawdata = mydatarec.data()
		
		if type(rawdata) == type(array(list())):
			if runfromhome:
				mytimes = mytree.get('dim_of(' + node + ')').data()
				rawdata = mytree.get(node).data()
			# print(mytimes,data)
			else:
				mytimes = mydatarec.dim_of().data()
				rawdata = mydatarec.data()
			# print(mean(data),data)
			# print(mytimes)
			# print(type(mytimes),type(data))
			start_idx = 0
			stop_idx = len(mytimes)
			if t_start is not None:
				start_idx = abs(mytimes - t_start).argsort()[0]
			if t_end is not None:
				stop_idx = abs(mytimes - t_end).argsort()[0]
			# print('test4')
			
			if times is None:
				# print(len(mytimes),len(data),start_idx,stop_idx)
				# print(len(mytimes[start_idx:stop_idx]),len(data[start_idx:stop_idx]),stop_idx-start_idx)
				return array(mytimes[start_idx:stop_idx]), array(rawdata[start_idx:stop_idx])
			else:
				idx1 = abs((mytimes - times[0])).argsort()[0]
				idx2 = abs((mytimes - times[-1])).argsort()[0] + 1
				if idx2 - idx1 == len(times) and allclose(mytimes[idx1:idx2], times):
					data = rawdata[idx1:idx2]
				elif min(times) >= min(mytimes) and max(times) <= max(mytimes):
					data = interp1d(mytimes, rawdata)(times)
					return array(data[start_idx:stop_idx])
				else:
					data = zeros(len(times))  # *NaN
					try:
						subtimes = (times[times >= mytimes[0]])[times[times >= mytimes[0]] <= mytimes[-1]]
					except:
						print(sys.exc_info())
						print(shape(times), shape(mytimes))
						print(min(times), min(mytimes), max(times), max(mytimes))
						print(shape(times[mytimes[0] >= times]))
					st_idx1 = (where(times >= mytimes[0]))[0][0]
					st_idx2 = (where(times <= mytimes[-1]))[0][-1]
					print(min(subtimes), min(mytimes), max(subtimes), max(mytimes))
					interpdata = interp1d(mytimes, rawdata)(subtimes)
					data[st_idx1:st_idx2 + 1] = interpdata
					return array(data)
		else:
			return rawdata
	except:
		print(f'problem occurred retrieving data for node {node}')
		return (np.array([np.nan]), np.array([np.nan]))
	

def get_nodedata(temp_node, temp_shot, temp_tStart=None, temp_tStop=None, temp_times=None):
	'''
	Branching wrapper that sometimes works better than "get_data" when accessing a node which contains only an integer or float... but try "get_data" first!
	Arguments:
		<temp_node>: tree node pointer/name (string)
		<temp_shot>: shot number target (int)
	Optional Arguments: if none of the following is used, data is returned without a timebase
		<temp_tStart>: time IN SECONDS according to node's timebase, from which to start reading data; usually used together with <temp_tStop> (float)
		<temp_tStop>: time IN SECONDS according to node's timebase, at which to stop reading data; usually used together with <temp_tStart> (float)
		<temp_times>: range of times IN SECONDS, according to node's timebase, at which to read data; overrides <temp_tStart> and <temp_tStop> (iterable of float)
	Returns:
		Contents of node
	'''
	if runfromhome:
		conn = mds.Connection('spitzer.ap.columbia.edu:8003')
		conn.openTree('hbtep2', shot)
		data = conn.get(node).data()
		return data
	else:
		return get_data_fromlab(Tree('hbtep2', temp_shot), temp_node, temp_tStart, temp_tStop, temp_times)


def get_data_old(temp_node, temp_shot, temp_tStart=None, temp_tStop=None, temp_times=None):
	'''
	Older version of "get_data"; deprecated but not erased, just in case some old code needs its format.  Interfaces more like Niko's old "get_data."
	Arguments:
		<temp_node>: tree node pointer/name (string)
		<temp_shot>: shot number target (int)
	Optional Arguments: if none of the following is used, data is returned without a timebase
		<temp_tStart>: time IN SECONDS according to node's timebase, from which to start reading data; usually used together with <temp_tStop> (float)
		<temp_tStop>: time IN SECONDS according to node's timebase, at which to stop reading data; usually used together with <temp_tStart> (float)
		<temp_times>: range of times IN SECONDS, according to node's timebase, at which to read data; overrides <temp_tStart> and <temp_tStop> (iterable of float)
	Returns:
		Return of proper read method
	'''
	# print('>>> DIAGNOSTIC: Getting data from'+temp_node)
	if runfromhome:
		return get_data_fromhome(temp_node, temp_shot, temp_tStart, temp_tStop)
	else:
		return read_data(Tree('hbtep2', temp_shot), temp_node, temp_tStart, temp_tStop, temp_times)


def signals_from_shot(tree, nodenames, t_start=None, t_end=None, times=None):
	'''Reads data in <nodenames> (iterable of strings/node pointers) from <tree> (tree pointer)

	This is just a thin wrapper to iterate over `get_data`.
	----- code by Niko Rath, modified for use with `get_data` by Paul Hughes
	'''
	
	if (times is not None and
			(t_start is not None or t_end is not None)):
		raise ValueError('t_start/t_end and times parameters are mutually exclusive')
	
	ret_times = times is None
	for (i, name) in enumerate(nodenames):
		if times is None:
			(times, data) = get_data(tree, name, t_start, t_end)
		else:
			data = get_data(tree, name, times=times)
		
		if i == 0:
			alldata = np.empty((len(nodenames), len(times)))
		
		alldata[i] = data
	
	if ret_times:
		return (times, alldata)
	else:
		return alldata


def GetPickupSubtraction(data, pickupdata, window_inds=None, t_st_sub=None, t_end_sub=None, t_start=None, t_end=None,
                         timebase=None):
	if window_inds == None:
		if t_st_sub == None:
			substartindex = 0
		else:
			substartindex = argmin(abs(timebase - t_st_sub))
		if t_end_sub == None:
			subendindex = len(data)
		else:
			subendindex = argmin(abs(timebase - t_end_sub))
		window_inds = [substartindex, subendindex]
	ratio = mean(data[window_inds[0]:window_inds[1]] / pickupdata[window_inds[0]:window_inds[1]])
	#	print(ratio)
	if t_start == None:
		startindex = 0
	else:
		startindex = argmin(abs(timebase - t_start))
	if t_end == None:
		endindex = len(data)
	else:
		endindex = argmin(abs(timebase - t_end))
	if 1:  # timebase==None:
		return (data - ratio * pickupdata)
		#	else:
		return timebase[startindex:endindex], (data - ratio * pickupdata)[startindex:endindex]


def DoZeroSubtract(time, data, myzerosubmethod='l'):
	if 'c' in myzerosubmethod:
		data -= mean(data[argmin(abs(time + 0.02)):argmin(abs(time - 0.100))])
	elif 'p' in myzerosubmethod:
		wind1ind1 = argmin(abs(time - 0.0))
		wind1ind2 = argmin(abs(time - 0.350))
		wind2ind1 = argmin(abs(time - 0.750))
		wind2ind2 = argmin(abs(time - 0.900))
		samp_times = append(mean(time[wind1ind1:wind1ind2]), mean(time[wind2ind1:wind2ind2]))
		samp_data = append(mean(data[wind1ind1:wind1ind2]), mean(data[wind2ind1:wind2ind2]))
		fittime, fit_y, y_rmserr = GetPolyfitCurve(samp_times, samp_data, 1, fit_x=time, xpoints=len(time))
		data -= fit_y
	elif 'l' in myzerosubmethod:
		wind1ind1 = argmin(abs(time - 0.00))
		wind1ind2 = argmin(abs(time - 0.10))
		wind2ind1 = argmin(abs(time - 0.90))
		wind2ind2 = argmin(abs(time - 0.95))
		samp_times = append(mean(time[wind1ind1:wind1ind2]), mean(time[wind2ind1:wind2ind2]))
		samp_data = append(mean(data[wind1ind1:wind1ind2]), mean(data[wind2ind1:wind2ind2]))
		fittime = time
		fit_slope = (samp_data[1] - samp_data[0]) / (samp_times[1] - samp_times[0])
		fit_y = mean(data[argmin(abs(time - 0.)):argmin(abs(time - 0.01))]) * 1. + fit_slope * time
		data -= fit_y
	else:
		print('WARNING: No valid zero-subtraction (Const/Linear/Polyfit) selected!')
	return data


def GetCorrectedSignal(thistree, thisnode, thiscalib=1., zerosubmethod='l', smooth=True, subfactslist=[], I_coils=[],
                       t_start=None, t_end=None):
	time, data = get_data(thistree, thisnode)
	data *= thiscalib
	data = DoZeroSubtract(time, data, myzerosubmethod=zerosubmethod)
	data = NewSmoothing(data, winhalfwid=20, wtmethod='g')
	# print(shape(data))
	if t_start == None:
		startindex = 0
	else:
		startindex = argmin(abs(time - t_start))
	if t_end == None:
		endindex = len(time)
	else:
		endindex = argmin(abs(time - t_end))
	time = time[startindex:endindex + 1]
	data = data[startindex:endindex + 1]
	# print(startindex,endindex,shape(data))
	for substep in range(0, len(subfactslist)):
		# print(shape(data),shape(I_coils[substep]),shape(I_coils))
		data -= I_coils[substep] * subfactslist[substep]
	#	if t_start==None:
	#		startindex = 0
	#	else:
	#		startindex = argmin(abs(time-t_start))
	#	if t_end==None:
	#		endindex = len(time_IP)
	#	else:
	#		endindex = argmin(abs(time-t_end))
	#	return time[startindex:endindex+1],data[startindex:endindex+1]
	return time, data


def GetPlasmaCurrent(mytree, coilcurrents=None, t_start=None, t_end=None):
	'''
	Calculates plasma current according to legacy code from Jeff, Daisuke, and Niko
	Arguments:
		<mytree>: MDSPlus tree pointer for shot in question (MDSPlus tree pointer)
		<myshot>: shot number (integer)
		<t_start>: time IN SECONDS from which to start reading data (float)
		<t_end>: time IN SECONDS at which to stop reading data (float)
	Optional Arguments:
		<fixedpickupparams>: identical to the output of "GetFixedPickups()"; can be supplied by calling routine to avoid recalculating fixed pickups (list of lists of float)
					if neglected, will rerun "GetFixedPickups()" and add up to a few seconds to runtime-per-shot
	Returns:
		<timebase_IP>: plasma current time base IN SECONDS (array of float)
		<plasmacurrent>: plasma current IN AMPERES (array of float)
	'''
	if coilcurrents == None:
		coilcurrents = []
		for coilstep in range(0, len(nodes_I_coil)):
			newtime, newcurrent = GetCorrectedSignal(mytree, nodes_I_coil[coilstep])
			coilcurrents.append(newcurrent)
	# print('gimme coil currents bruh')
	coilpickupfacts = [1.568, 0., 0.405, 0., 0., 0.35, 0.20, 0., 0.]
	# for coilstep in range(0,len(nodes_I_coil)):
	# print('For plasma current: ',shape(coilcurrents[coilstep]))
	time_IP, data_IP = GetCorrectedSignal(mytree, node_IP, I_coils=coilcurrents, subfactslist=coilpickupfacts,
	                                      t_start=t_start, t_end=t_end)
	
	if t_start == None:
		startindex = 0
	else:
		startindex = argmin(abs(time_IP - t_start))
	if t_end == None:
		endindex = len(time_IP)
	else:
		endindex = argmin(abs(time_IP - t_end))
	return time_IP[startindex:endindex + 1], data_IP[startindex:endindex + 1]


def GetTimescales(shot_tree):
	'''
	Returns bank/diagnostic timebases
	Arguments:
		<shot_tree>: tree pointer (class 'MDSplus.tree.Tree')
	Returns:
		<timescale>: timebase for most digitizers (array of float)
		<timescale_TF>: timebase for basement rack (array of float)
	'''
	timescale = get_data(shot_tree, VF_node_loc)[0]
	last_index = len(timescale) - 10
	timescale = timescale[:last_index + 1]
	
	timescale_TF = get_data(shot_tree, TFP_node_loc)[0]
	last_index_TF = len(timescale_TF) - 10
	timescale_TF = timescale_TF[:last_index_TF + 1]
	
	return timescale, timescale_TF


def GetPickups(shotnum, fixedparams, mydigitimes, myTFdigitimes):
	'''
	Returns standard coil pickup calibration values for shot <shotnum>, important for major radius and plasma current calculations
		--- this code simplified and condensed from equivalent code in hbtep.misc
	Arguments:
		<shotnum>: shot number (integer)
		<fixedparams>: standard calibration values, see output of "GetFixedPickups" below
			Allows parent code to run GetFixedPickups in header section and calculate these constant terms only once per run
		<mydigitimes>: standard timebase parameters for non-TF banks and diagnostics (list of integer and float)
			Organized as number of time points, time step length IN SECONDS, and time value of first point
		<myTFdigitimes>: standard timebase parameters for TF bank (list of integer and float)
			Organized as number of time points, time step length IN SECONDS, and time value of first point
	Returns:
		(list of lists of float): calculated pickup from OH, VF, and TF on sin1 and cos1 rogowskis, and from TF on IP rogowski, on shot <shotnum>
	'''
	[datalength, deltat, tnaught] = [len(mydigitimes), mydigitimes[1] - mydigitimes[0], mydigitimes[0]]
	[datalength_TF, deltat_TF, tnaught_TF] = [len(myTFdigitimes), myTFdigitimes[1] - myTFdigitimes[0], myTFdigitimes[0]]
	[IP_tf_pickup, sin1_tf_pickup, cos1_tf_pickup, t_ohbias_trig, oh_start_trig, sin1_oh_raw, cos1_oh_raw,
	 OHcurrent_raw, rog_tb, sin1_vf_raw, cos1_vf_raw, VFcurrent_raw] = fixedparams
	
	# TF pickup
	# IP_tf_pickup = get_data(IP_raw_loc_old,tf_pickup_shot)[0]
	# sin1_tf_pickup = get_data(IP_raw_loc_old,tf_pickup_shot)[0]
	# cos1_tf_pickup = get_data(IP_raw_loc_old,tf_pickup_shot)[0]
	#	sin1_tf_pickup = sum(sin1_tf/TorMagField)/datalength
	#	cos1_tf_pickup = sum(cos1_tf/TorMagField)/datalength
	
	# OH pickup
	# t_ohbias_trig=get_data('.timing.banks:oh_bias_st',oh_pickup_shot)*u
	# oh_start_trig=get_data('.timing.banks:oh_st',oh_pickup_shot)*u
	
	offset_index = TtI(t_ohbias_trig, deltat, tnaught)  # time lt oh_bias_trig)
	index_0 = TtI(-m, deltat, tnaught)
	index_1 = TtI(20. * m, deltat, tnaught)
	if diagnosticsyes: print(deltat, tnaught)
	rog_tb_trim = rog_tb[index_0:index_1]
	if (oh_pickup_shot > 70000):
		sin1_oh = get_data(sin1_loc, oh_pickup_shot)[1]
		cos1_oh = get_data(cos1_loc, oh_pickup_shot)[1]
		OHcurrent = get_data(OH_node_loc, oh_pickup_shot)[1]
	else:
		#		CallError('Need to fix pauls_toolbox.GetPickups to analyze shots before 70000')
		#		exit()
		# sin1_oh_raw = -get_data(sin1_raw_loc_old,oh_pickup_shot)[0]
		sin1_oh_raw -= sum(sin1_oh_raw[0:offset_index]) / offset_index
		sin1_oh_raw = sin1_oh_raw[index_0:index_1]
		if diagnosticsyes: print(sin1_oh_raw, index_0, index_1)
		sin1_oh = SimpleIntegral(sin1_oh_raw, deltat, RC=RC_rog_sin1)
		# cos1_oh_raw = get_data(cos1_raw_loc_old,oh_pickup_shot)[0]
		cos1_oh_raw -= sum(cos1_oh_raw[0:offset_index]) / offset_index
		cos1_oh_raw = cos1_oh_raw[index_0:index_1]
		cos1_oh = SimpleIntegral(cos1_oh_raw, deltat, RC=RC_rog_cos1)
		# OHcurrent_raw = get_data(OH_raw_loc_old,oh_pickup_shot)[0]
		OHcurrent_raw -= sum(OHcurrent_raw[0:offset_index]) / offset_index
		OHcurrent_raw = OHcurrent_raw[index_0:index_1]
		OHcurrent = SimpleIntegral(OHcurrent_raw, deltat, RC=RC_rog_OH) * OH_mult_factor
	sin1_oh_pickup = norm(sin1_oh[:TtI(m, deltat, rog_tb_trim[0])]) / norm(OHcurrent[:TtI(m, deltat, rog_tb_trim[0])])
	cos1_oh_pickup = norm(cos1_oh[:TtI(m, deltat, rog_tb_trim[0])]) / norm(OHcurrent[:TtI(m, deltat, rog_tb_trim[0])])
	#	sin1_oh_pickup = sum(sin1_oh/OHcurrent)/len(OHcurrent)
	#	cos1_oh_pickup = sum(cos1_oh/OHcurrent)/len(OHcurrent)
	
	# VF pickup
	if (vf_pickup_shot > 70000):
		sin1_vf = get_data(sin1_loc + ':raw', vf_pickup_shot)[1]
		cos1_vf = get_data(cos1_loc + ':raw', vf_pickup_shot)[1]
		VFcurrent = get_data(VF_node_loc + ':raw', vf_pickup_shot)[1]
	else:
		#		CallError('Need to fix pauls_toolbox.GetPickups to analyze shots before 70000')
		#		exit()
		# sin1_vf_raw = -get_data(sin1_raw_loc_old,vf_pickup_shot)[0]
		sin1_vf_raw -= sum(sin1_vf_raw[0:offset_index]) / offset_index
		sin1_vf_raw = sin1_vf_raw[index_0:index_1]
		sin1_vf = SimpleIntegral(sin1_vf_raw, deltat, RC=RC_rog_sin1)
		# cos1_vf_raw = get_data(cos1_raw_loc_old,vf_pickup_shot)[0]
		cos1_vf_raw -= sum(cos1_vf_raw[0:offset_index]) / offset_index
		cos1_vf_raw = cos1_vf_raw[index_0:index_1]
		cos1_vf = SimpleIntegral(cos1_vf_raw, deltat, RC=RC_rog_cos1)
		# VFcurrent_raw = get_data(VF_raw_loc_old,vf_pickup_shot)[0]
		VFcurrent_raw -= sum(VFcurrent_raw[0:offset_index]) / offset_index
		VFcurrent = VFcurrent_raw[index_0:index_1]
	sin1_vf_pickup = norm(sin1_vf[:TtI(m, deltat, rog_tb_trim[0])]) / norm(VFcurrent[:TtI(m, deltat, rog_tb_trim[0])])
	cos1_vf_pickup = norm(cos1_vf[:TtI(m, deltat, rog_tb_trim[0])]) / norm(VFcurrent[:TtI(m, deltat, rog_tb_trim[0])])
	#	sin1_vf_pickup = sum(sin1_vf/VFcurrent)/len(VFcurrent)
	#	cos1_vf_pickup = sum(cos1_vf/VFcurrent)/len(VFcurrent)
	cos1_vf_pickup = 0.0046315133 * (-m)
	cos1_oh_pickup = 7.0723416e-08
	return [sin1_oh_pickup, cos1_oh_pickup], [sin1_vf_pickup, cos1_vf_pickup], [IP_tf_pickup, sin1_tf_pickup,
	                                                                            cos1_tf_pickup]


def GetFixedPickups(verbose=True):
	'''
	Returns reference values for calibrating sin1, cos1, and IP rogowski pickup to coil currents
		--- this code simplified and condensed from equivalent code in hbtep.misc
	Arguments:
		<verbose>: when True, prints out status of fixed pickup readings (boolean)
	Returns:
		<IP_tf_pickup>: TF bank pickup on IP rogowski for TF pickup reference shot
		<sin1_tf_pickup>: TF bank pickup on sin1 rogowski for TF pickup reference shot
		<cos1_tf_pickup>: TF bank pickup on cos1 rogowski for TF pickup reference shot
		<t_ohbias_trig>: time at which OHB bank fires for OH pickup reference shot
		<oh_start_trig>: time at which OHS bank fires for OH pickup reference shot
		<sin1_oh_raw>: raw OH bank pickup on sin1 rogowski for OH pickup reference shot
		<cos1_oh_raw>: raw OH bank pickup on cos1 rogowski for OH pickup reference shot
		<OHcurrent_raw>: raw OH current for OH pickup reference shot
		<rog_tb>: timebase for <OHcurrent_raw>
		<sin1_vf_raw>: raw VF bank pickup on sin1 rogowski for VF pickup reference shot
		<cos1_vf_raw>: raw VF bank pickup on sin1 rogowski for VF pickup reference shot
		<VFcurrent_raw>: raw VF current for VF pickup reference shot
	NOTE:
		Best used by running at beginning of code and assigning return values to a constant; streamlines code by getting picking values only one time per run
		Mainly used for calculation of plasma current
	'''
	fixedpickups_st = timetools.clock()
	if verbose: print('>>> Initializing coil pickup parameters... ')
	# TF pickup
	tftree = get_tree_conn(tf_pickup_shot)
	IP_tf_pickup = get_data(tftree, IP_raw_loc_old)[1]
	if verbose: print(' > > > IP-TF')
	sin1_tf_pickup = get_data(tftree, IP_raw_loc_old)[1]
	if verbose: print(' > > > sin1-TF')
	cos1_tf_pickup = get_data(tftree, IP_raw_loc_old)[1]
	if verbose: print(' > > > cos1-TF')
	
	ohtree = get_tree_conn(oh_pickup_shot)
	t_ohbias_trig = get_data(ohtree, '.timing.banks:oh_bias_st') * u
	if verbose: print(' > > > t_OHB')
	oh_start_trig = get_data(ohtree, '.timing.banks:oh_st') * u
	if verbose: print(' > > > t_OHS')
	
	rog_tb = get_data(ohtree, OH_raw_loc_old)[0]
	if verbose: print(' > > > OH timebase')
	
	sin1_oh_raw = -get_data(ohtree, sin1_raw_loc_old)[1]
	if verbose: print(' > > > raw sin1-OH')
	cos1_oh_raw = get_data(ohtree, cos1_raw_loc_old)[1]
	if verbose: print(' > > > raw cos1-OH')
	OHcurrent_raw = get_data(ohtree, OH_raw_loc_old)[1]
	if verbose: print(' > > > raw I_OH')
	
	#	rog_tb = get_data(ohtree,OH_raw_loc_old)[0]
	#	print('    >>> rog_tb')
	
	vftree = get_tree_conn(vf_pickup_shot)
	sin1_vf_raw = -get_data(vftree, sin1_raw_loc_old)[1]
	if verbose: print(' > > > raw sin1-VF')
	cos1_vf_raw = get_data(vftree, cos1_raw_loc_old)[1]
	if verbose: print(' > > > raw cos1-VF')
	VFcurrent_raw = get_data(vftree, VF_raw_loc_old)[1]
	if verbose: print(' > > > raw I_VF')
	
	if verbose: print(' > > > Took {0} seconds'.format(timetools.clock() - fixedpickups_st))
	return IP_tf_pickup, sin1_tf_pickup, cos1_tf_pickup, t_ohbias_trig, oh_start_trig, sin1_oh_raw, cos1_oh_raw, OHcurrent_raw, rog_tb, sin1_vf_raw, cos1_vf_raw, VFcurrent_raw


def calc_r_major(tree, t_start=None, t_end=None, times=None, tf_sub=False):
	'''
	Calculate major radius ----- code by Niko Rath
		(based in part on IDL code and major radius analysis by Jeff Levesque and Daisuke Shiraki); Jeff is the resident expert on this code!
	Arguments:
		<tree>: MDSPlus tree instance (class 'MDSplus.tree.Tree')
	Optional Arguments: if none of the following is used, data is returned without a timebase
		<t_start>: time IN SECONDS according to node's timebase, from which to start reading data; usually used together with t_end (float)
		<t_end>: time IN SECONDS according to node's timebase, at which to stop reading data; usually used together with t_start (float)
		<times>: range of times IN SECONDS, according to node's timebase, at which to read data; overrides t_start and t_end (iterable of float)
	Returns:  format depends on timebase arguments...  if <times> is used, only <r_major> is returned; else (<times>, <r_major>) is returned
		<r_major>: vector of major radius values IN METERS at requested time range or time points, or all points if not specified (array of float)
		<times>: vector of time values IN SECONDS at requested time range, only returned if timebase arguments are used (array of float)
	'''
	
	# Determined by Daisuke during copper plasma calibration
	a = .00643005
	b = -1.10423
	c = 48.2567
	
	# Calculated by Jeff, but still has errors
	vf_pickup = 0.0046315133 * -1e-3
	oh_pickup = 7.0723416e-08
	
	return_times = times is None
	
	# Read raw signal
	if times is None:
		raw_end = t_end
	else:
		raw_end = times[-1] + 2 * (times[1] - times[0])
	(cos1_times, cos1_raw) = get_data(tree, '.sensors.rogowskis:cos_1:raw',
	                                  t_end=raw_end)
	
	# Subtract offset
	offset_time = 0
	offset_mask = cos1_times > offset_time
	cos1_raw -= cos1_raw[~offset_mask].mean()
	
	if tf_sub:
		# tree2 = MDSplus.Tree('hbtep2', 74453) # 6.1kV reference shot
		tree2 = get_tree_conn(74453)
		
		# Pickup
		cos1_p_raw = get_data(tree2, '.devices.west_rack:cpci:input_02',
		                      times=cos1_times[offset_mask])
		
		# Offset in pickup
		#       pre_trig_cos1 = get_data(tree2, '.devices.west_rack:cpci:input_02',
		#                                  t_end=tree2.getNode('.timing.banks:tf_st').data() * 1e-6)[1]
		ptc1_endtime = get_data(tree2, '.timing.banks:tf_st') * 1e-6
		pre_trig_cos1 = get_data(tree2, '.devices.west_rack:cpci:input_02', t_end=ptc1_endtime)[1]
		if len(pre_trig_cos1) > 0:
			cos1_p_raw -= pre_trig_cos1.mean()
		else:
			raise RuntimeError('No pre-TF data available for offset subtraction in pickup')
		
		# Subtract corrected offset
		cos1_raw[offset_mask] -= cos1_p_raw
	
	# Integrate
	cos1RC = get_data(tree, '.sensors.rogowskis:cos_1:rc_time')
	cos1 = (cumtrapz(cos1_raw, cos1_times)
	        + cos1_raw[:-1] * cos1RC)
	cos1_times = cos1_times[:-1]
	
	if times is not None:
		cos1 = interp1d(cos1_times, cos1)(times)
	elif t_start is not None:
		idx = abs(cos1_times - t_start).argsort()[0]
		times = cos1_times[idx:]
		cos1 = cos1[idx:]
	else:
		times = cos1_times
	
	vf = get_data(tree, '.sensors.vf_current', times=times)
	oh = get_data(tree, '.sensors.oh_current', times=times)
	ip = get_data(tree, '.sensors.rogowskis:ip', times=times)
	# ip *= tree.getNode('.sensors.rogowskis:ip:gain').data()
	ip *= get_data(tree, '.sensors.rogowskis:ip:gain')
	
	if len(vf) > len(oh):
		vf = vf[:len(oh)]
	else:
		oh = oh[:len(vf)]
	
	# print('DIAGNOSTIC: vf = ',len(vf))
	# print('DIAGNOSTIC: vf_pickup = ',vf_pickup)
	# print('DIAGNOSTIC: oh = ',len(oh))
	# print('DIAGNOSTIC: oh_pickup = ',oh_pickup)
	# print('DIAGNOSTIC: ip = ',len(ip))
	pickup = vf * vf_pickup + oh * oh_pickup
	# print('DIAGNOSTIC: pickup = ',len(pickup))
	# print('DIAGNOSTIC: cos1 = ',len(cos1))
	ratio = ip / (cos1 - pickup)
	
	arg = b ** 2 - 4 * a * (c - ratio)
	arg[arg < 0] = 0
	r_major = (-b + sqrt(arg)) / (2 * a)
	r_major /= 100  # Convert to meters
	
	if return_times:
		return (times, r_major)
	else:
		return r_major


def calc_r_minor(r_major, vpos=None):
	'''Calculate minor radius ----- code by Niko Rath'''
	'''vertical position sensitivity added by Paul Hughes'''
	'''
	Arguments:
		<r_major>: major radius IN METERS (array of float)
	Optional Arguments: if <vpos>=None, vertical position is assumed to be zero at all time points
		<vpos>: vertical position relative to chamber centerline IN METERS (array of float, or None)
	Returns:
	<r_minor>: minor radius IN METERS; piggybacks on r_major timebase (array of float)
	'''
	
	if vpos is None:
		vpos = zeros(len(r_major))
	
	r_minor = empty_like(r_major)
	r_minor[:] = 0.15 - abs(vpos)  # Up/down limited
	
	mask = r_major > (0.92) + abs(vpos)
	r_minor[mask] = 1.07 - r_major[mask]  # Outboard limited
	
	mask = r_major < (0.92 - 0.01704) - abs(vpos)
	r_minor[mask] = r_major[mask] - 0.75296  # Inboard limited
	
	mask = r_minor < 0.
	r_minor[mask] = zeros(len(r_minor))[mask]
	
	return r_minor


def CalcVertPos(this_shot, this_IP, this_timebase, method='foursensor'):
	'''
	Calculates vertical position using top- and bottom-most chamber-mounted PA poloidal sensors (4 sensors in total)
	Arguments:
		<this_shot>: shot number (integer)
		<this_IP>: plasma current trace (array of float)
		<this_timebase>: timebase to use for vertical position, typically IP timebase (array of float)
	Optional Arguments:
		<method>: currently has two options; however, only 'foursensor' has actually been implemented properly!
			'foursensor': uses top- and bottom-most chamber-mounted PA1 and PA2 sensors to measure vertical imbalance of plasma current
			'halfsine': uses full chamber-mounted set of PA1p and PA2p sensors to measure vertical imbalance of plasma current; NOT IMPLEMENTED YET
	Returns:
		<vertpos>: vertical position IN METERS relative to machine midplane; piggybacks on timebase provided (array of float, of length len(this_IP))
	'''
	this_tree = get_tree_conn(this_shot)
	this_IP = this_IP[:len(this_timebase)]
	if method == 'foursensor':
		PA1Bnode = 'PA1_S08P'
		PA1Tnode = 'PA1_S25P'
		PA2Bnode = 'PA2_S08P'
		PA2Tnode = 'PA2_S25P'
		fullsensorlist = ['.sensors.magnetic:%s' % x for x in [PA1Bnode, PA1Tnode, PA2Bnode, PA2Tnode]]
		#	(garbagetimes, thissignalset) = signals_from_shot(this_tree, fullsensorlist, t_start=starttime, t_end=endtime )
		PA1Bsig, PA1Tsig, PA2Bsig, PA2Tsig = signals_from_shot(this_tree, fullsensorlist, times=this_timebase)
		# B = mu0*IP/2piR  ...  R = 2x10^-7 IP / B
		#	print(len(this_IP),len(PA1Tsig),len(PA1Bsig),len(PA2Tsig),len(PA2Bsig))
		vertpos = (1. / (2. * pi)) * mu_0 * this_IP * (
				(1. / abs(PA1Tsig) - 1. / abs(PA1Bsig)) + (1. / abs(PA2Tsig) - 1. / abs(PA2Bsig))) / 4.
	elif method == 'halfsine':
		# *** get PA1 signals
		PA1set_full = ['PA1_S{0:02}P'.format(thissensor) for thissensor in range(1, 8) + range(26, 33)]
		PA1set = [x for x in PA1set_full if x not in sensorblacklist]
		PA1nodes = ['.sensors.magnetic:%s' % x for x in PA1set]
		thetaPA1p = []
		PA1p_skips = 0
		for PA1step in range(0, len(PA1set)):
			if not (PA1set[PA1step] == PA1set_full[PA1step + PA1p_skips]):
				PA1p_skips += 1
			thetaPA1p.append((PA1step + PA1p_skips) * 2. * pi / 32. + 0.0001)
		PA1sigs = signals_from_shot(this_tree, PA1nodes, times=this_timebase)
		# *** get PA2 signals
		PA2set_full = ['PA1_S{0:02}P'.format(thissensor) for thissensor in range(1, 8) + range(26, 33)]
		PA2set = [x for x in PA2set_full if x not in sensorblacklist]
		PA2nodes = ['.sensors.magnetic:%s' % x for x in PA2set]
		thetaPA2p = []
		PA2p_skips = 0
		for PA2step in range(0, len(PA2set)):
			if not (PA2set[PA2step] == PA2set_full[PA2step + PA2p_skips]):
				PA2p_skips += 1
			thetaPA2p.append((PA2step + PA2p_skips) * 2. * pi / 32. + 0.0001)
		PA2sigs = signals_from_shot(this_tree, PA2nodes, times=this_timebase)
		PA1sineproj = array([sum(PA1sigs[:, tstep] * sin(thetaPA1p)) / sum(sin(thetaPA1p) ** 2.) for tstep in
		                     range(0, len(this_timebase))])
		PA2sineproj = array([sum(PA2sigs[:, tstep] * sin(thetaPA2p)) / sum(sin(thetaPA2p) ** 2.) for tstep in
		                     range(0, len(this_timebase))])
		vertpos = (1. / (2. * pi)) * mu_0 * (PA1sineproj / 2. + PA2sineproj / 2.)
	return vertpos


def calc_q(tree, times, r_major=None, r_minor=None, I_P=None):
	'''Calculate safety factor ----- code by Niko Rath'''
	'''
	Arguments:
		<tree>: tree pointer (class 'MDSplus.tree.Tree')
		<times>: timebase to use IN SECONDS (array of float)
	Optional Arguments:
		<r_major>: major radius IN METERS (array of float) - if neglected, is called by calc_q
		<r_minor>: minor radius IN METERS (array of float) - if neglected, is called by calc_q WITHOUT vertical position correction
		<I_P>: plasma current IN AMPERES (array of float) - if neglected, is called by calc_q
	Returns:
	<q>: edge safety factor (array of float)
	Note:
		TF has been measured at roughly 15% higher than TF probe data; can correct for this by multiplying return value by 1.15, if desired
	'''
	
	if r_major is None:
		r_major = calc_r_major(tree, times=times)
	
	if r_minor is None:
		r_minor = calc_r_minor(r_major)
	
	ip = get_data(tree, '.sensors.rogowskis:ip', times=times)
	
	# FIXME: TF probe has pickup from VF and OH
	tf = get_data(tree, '.sensors.tf_probe', times=times)
	tf *= get_data(tree, '.sensors:tf_probe:mr') / r_major
	
	q = r_minor ** 2 * tf / (2e-7 * ip * r_major)
	
	if I_P is None:
		I_P = 50000. * ones(len(q))
	mask = I_P[:len(q)] < 6000.
	q[mask] = k * ones(len(mask))
	
	return q


def GetCouplingToTA(R_0, r_a, Z_0, mnum=3.):
	'''
	Calculates coupling coefficient (i.e. cylindrical harmonic dropoff) of observed mode amplitude based on plasma major radius, minor radius, and vertical position
	Arguments:
		<R_0>: major radius IN METERS (array of float)
		<r_a>: minor radius IN METERS (array of float)
		<Z_0>: vertical position IN METERS (array of float)
	Returns:
		<TAcoupling>: Coupling coefficient (a/b)**2m for the toroidal sensor array; corrected for a given m-number by taking <TAcoupling>**mnum (array of float)
	'''
	R_TA = 0.90 - 0.155
	r_TA = sqrt((R_0 - R_TA) ** 2. + Z_0 ** 2.)
	TAcoupling = (r_a / r_TA) ** 2.
	return TAcoupling


def get_puff_pressures(thistree, returntimetrace=False):
	'''
	Calculates base pressure (i.e. pressure prior to gas injection) and fill pressure (i.e. pressure following gas injection) in chamber from ion gauge signal
		MAY work remotely, has not been tested as of Mar. 21st 2016
	Arguments:
		<thistree>: MDSPlus tree pointer for shot in question (MDSPlus tree pointer)
	Returns:
		(float): digitizer output puff length in microseconds
		<meanbasepress>: base pressure in nTorr (float)
		<meanfillpress>: fill pressure in nTorr (float)
	'''
	#	nodeloc_pufftrig = '.devices.basement:a14_4:input_6'
	nodeloc_iongauge = '.devices.basement:a14_04:input_5'
	voltagedivider = 2.0286  # ********************************** INCLUDES CORRECTION FOR A14 INPUT IMPEDANCE!
	t_puff_1 = (get_data(thistree, '.timing:gas_puff')[1][0]) * 10 ** (-6)
	t_puff_2 = (get_data(thistree, '.timing:gas_puff')[1][1]) * 10 ** (-6)
	#	t_puff_1,t_puff_2 = ((thistree.getNode('.timing:gas_puff')).data())*10**(-6)
	t_basementA14start = (get_data(thistree, '.timing.digitizer:basement_tf')) * 10 ** (-6)
	#	t_basementA14start = ((thistree.getNode('.timing.digitizer:basement_tf')).data())*10**(-6)
	t_tfbankstart = (get_data(thistree, '.timing.banks:TF_ST')) * 10 ** (-6)
	#	t_tfbankstart = ((thistree.getNode('.timing.banks:TF_ST')).data())*10**(-6)
	p_base_times = [t_puff_1 - 0.001 - 1. / 60., t_puff_1 - 0.001]
	p_fill_times = [t_tfbankstart - 0.006 - 1. / 60., t_tfbankstart - 0.006]
	p_IG_rawsig = (get_data(thistree, nodeloc_iongauge))[1]
	t_p_IG = (get_data(thistree, nodeloc_iongauge))[0]
	#	p_IG_rawsig = (thistree.getNode(nodeloc_iongauge)).data()
	#	t_p_IG = (thistree.getNode(nodeloc_iongauge)).dim_of().data()
	#	DGN_note(lineno(),['t_p_IG','p_base_times','p_fill_times'],[t_p_IG,p_base_times,p_fill_times])
	dt = t_p_IG[1] - t_p_IG[0]
	p_IG_smth = NewSmoothing(p_IG_rawsig[int(round((p_base_times[0] - t_p_IG[0]) / dt)) // 2:int(
		round((p_fill_times[1] - t_p_IG[0]) / dt)) * 2], 250, 'g')
	p_IG_proc = 10 ** (abs(p_IG_smth) * (voltagedivider) - 2.)
	minindex = int(round((p_base_times[0] - t_p_IG[0]) / dt)) // 2
	meanbasepress = mean(p_IG_proc[int(round((p_base_times[0] - t_p_IG[0]) / dt)) - minindex:int(
		round((p_base_times[1] - t_p_IG[0]) / dt)) - minindex]) / 0.35
	meanfillpress = mean(p_IG_proc[int(round((p_fill_times[0] - t_p_IG[0]) / dt)) - minindex:int(
		round((p_fill_times[1] - t_p_IG[0]) / dt)) - minindex]) / 0.35
	if returntimetrace:
		return (t_puff_2 - t_puff_1) * 10 ** 6, meanbasepress, meanfillpress, p_IG_smth  # in nTorr
	else:
		return (t_puff_2 - t_puff_1) * 10 ** 6, meanbasepress, meanfillpress  # in nTorr


def get_disruption_time(t_IP, IP):
	'''
	Calculates the time of disruption based on properties of the plasma current trace.  Caveats: 1) if the last value in the plasma current vector is too high (i.e. the time base ends before the plasma has clearly disrupted), it returns a time 100us after the end of the time base; 2) if the plasma current spike is very weak, the function can become confused.
	Arguments:
		<t_IP>: plasma current time base IN SECONDS (array of float)
		<IP>: plasma current IN AMPERES (array of float)
	Returns:
		<currentpeaktime>: time of disruption (specifically, the current spike) IN SECONDS (float)
	'''
	dIPdt = NewSmoothing(SimpleDerivative(NewSmoothing(IP, 20, 'g', peakwid=6, edgetreat='s'), t_IP), 8, 'g', peakwid=2,
	                     edgetreat='s')
	if IP[-1] < 2000.:
		if max(dIPdt[argmin(abs(t_IP - 1.25 * m)):]) > 10 * M:
			riseindex = argmax(dIPdt)
			currentpeakindex = argmin(abs(dIPdt)[riseindex:riseindex + 10]) + riseindex
			currentpeaktime = t_IP[currentpeakindex]
		else:
			lastplasmaindex = argmin(dIPdt)
			riseindex = argmax(dIPdt[lastplasmaindex - 500:lastplasmaindex]) + (lastplasmaindex - 500)
			currentpeakindex = argmin(abs(dIPdt)[riseindex:riseindex + 10]) + riseindex
			currentpeaktime = t_IP[currentpeakindex]
	else:
		currentpeaktime = t_IP.max() + 0.1 * m
	return currentpeaktime


''' ****** no more code after this ****** '''

''' No, seriously, it's just a text buffer, from here on. '''

import matplotlib.pyplot as plt
import numpy as np
import os
from helpful_stuff import is_nbi_shot, get_tree_conn, get_data

def validate_nubeamify_inputs(dir, trdat, **kwargs):
	if len(trdat) != 15:
		return 0, f'TRDAT filename incorrect format. Got {trdat} but expected 105795A01TR.DAT'
	if 'output_file_letter' not in kwargs.keys():
		ofl = 'A'
		while os.path.isfile(f'{dir}{trdat[0:6]}{ofl}{trdat[7:]}'):
			if ofl == 'Z':
				return 0, 'No output file letter available'
			ofl = chr(ord(ofl)+1)  # increment letter
		trdat_out = f'{trdat[0:6]}{ofl}{trdat[7:]}'
		resp = input(f'No output file letter given; output file will be- {trdat_out}, is this ok? (Y/N)')
		if resp.lower() != 'y':
			return 0, 'User denied output_file_letter, terminating.'
	if not ('twin' and 'ltx_shots' in kwargs.keys() or 'ebeam' and 'pbeam' in kwargs.keys() or 'ebeam' and 'ibeam' in kwargs.keys() or 'pbeam' and 'ibeam' in kwargs.keys()):
		return 0, 'missing info: include either ltx_shots and twin, or two of three- ebeam, pbeam, ibeam'
	return 0, 'inputs ok...'


def get_nbi_inputs_for_trdat(ltx_shots, twin):
	einj, pinj = np.array([]), np.array([])
	for shot in ltx_shots:
		try:
			t = get_tree_conn(shot, treename='ltx_b')
			print('gathering data for shot {} occurred on {}'.format(shot, get_data(t, '.metadata:timestamp')))
			if is_nbi_shot(shot, t):
				(tibeam, ibeam) = get_data(t, '.oper_diags.ltx_nbi.source_diags.i_hvps')
				(tbeam, vbeam) = get_data(t, '.oper_diags.ltx_nbi.source_diags.v_hvps')
				ibeam = np.interp(tbeam, tibeam, ibeam)
				ii = np.where((tbeam >= twin[0]) & (tbeam <= twin[1]))
				einj = np.append(einj, np.mean(vbeam[ii]))
				pinj = np.append(pinj, np.mean(vbeam[ii]*ibeam[ii]))
			else:
				print(f'shot {shot} is not an NBI shot- omitting')
		except:
			print(f'problem occurred gathering data for shot {shot}: returning')
			return [np.nan]
	return [np.mean(einj), np.mean(pinj)]


def ltx_nubeamify(dir, trdat, **kwargs):
	# output_file_letter, ltx_shots, twin, ebeam, pbeam, ibeam
	ok, mess = validate_nubeamify_inputs(dir, trdat, **kwargs)
	if not ok:
		print(f"Problem encountered with inputs: {mess}")
		return
	else:
		print(mess)
	
	if 'ltx_shots' in kwargs.keys():
		ep_inj = get_nbi_inputs_for_trdat(kwargs['ltx_shots'], kwargs['twin'])
		if len(ep_inj) != 2:
			return
		else:
			kwargs['ebeam'], kwargs['pbeam'] = ep_inj[0], ep_inj[1]
	if 'ibeam' in kwargs.keys():
		kwargs['pbeam'] = kwargs['ebeam'] * kwargs['ibeam']
		
	# we should have ebeam and pbeam defined
	# go through original tr.dat and look for existing NBI fields, ask if ok to overwrite if present
	# copy original to new file with new output_file_letter and append NBI fields to end
	# if it works here, look to move over into bot (Marvin maybe?) and relocate all necessary subroutines
	
	
if __name__ == '__main__':
	dir = 'Z:/transp/t105795/'
	trdat = '105795A01TR.DAT'
	ltx_nubeamify(dir, trdat, ltx_shots=[105795], twin=[.465, .472])

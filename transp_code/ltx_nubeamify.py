import datetime
import warnings

import matplotlib.pyplot as plt
import numpy as np
import os
from helpful_stuff import is_nbi_shot, get_tree_conn, get_data


# def validate_nubeamify_inputs(dir, trdat, **kwargs):


def get_nbi_inputs_for_trdat(ltx_shots):
	einj, pinj = np.array([]), np.array([])
	pad = 0.25e-3  # distance to step twin away from beam turnon/turnoff
	ton, toff = np.array([]), np.array([])
	for shot in ltx_shots:
		try:
			t = get_tree_conn(shot, treename='ltx_b')
			print('gathering data for shot {} occurred on {}'.format(shot, get_data(t, '.metadata:timestamp')))
			if is_nbi_shot(shot, t):
				(tibeam, ibeam) = get_data(t, '.oper_diags.ltx_nbi.source_diags.i_hvps')
				(tbeam, vbeam) = get_data(t, '.oper_diags.ltx_nbi.source_diags.v_hvps')
				ibeam = np.interp(tbeam, tibeam, ibeam)
				tbeamon, tbeamoff = tbeam[np.where(vbeam > 1000)][0] + pad, tbeam[np.where(vbeam > 1000)][-1] - pad
				ton, toff = np.append(ton, tbeamon), np.append(toff, tbeamoff)
				ii = np.where((tbeam >= tbeamon + pad) & (tbeam <= tbeamoff - pad))
				einj = np.append(einj, np.mean(vbeam[ii]))
				pinj = np.append(pinj, np.mean(vbeam[ii] * ibeam[ii]))
			else:
				print(f'shot {shot} is not an NBI shot- omitting')
		except:
			print(f'problem occurred gathering data for shot {shot}: returning')
			return [np.nan]
	if max(ton) - min(ton) > 1.e-3:
		print('discrepancy in beam turnon times, returning.')
		return [np.nan]
	return [np.mean(einj), np.mean(pinj), np.mean(ton), np.mean(toff)]


def copy_content(fromfile, tofile):
	print(f'copying {fromfile} to {tofile}')
	content = ''
	readin = open(fromfile, 'r')
	for line in readin:
		content += line
	readin.close()
	writeto = open(tofile, 'w')
	writeto.write(content)
	writeto.close()


def get_beam_dict(eb, pb, twin):
	ts = datetime.datetime.today()
	bhead = '!==================================================\n! Neutral Beam\n! http://w3.pppl.gov/~pshare/help/body_transp_hlp.html#outfile196.html\n! Written by ltx_neubeamify.py on {ts.strftime("%m/%d/%Y at %I:%M%p")}\n!==================================================\n'
	bdict = {'NLBCCW': '.T       ! .T for CCW BT',
	         'NLJCCW': '.F       ! .T for CCW Ip',
	         'SELAVG': "'FBM BDENS2 BDENSS3 BMVOL EBA2PL EBA2PP' !quantities	to average",
	         'OUTTIM': ', '.join(
		         [str(round(t, 5)) for t in np.linspace(twin[0], 8 / 5 * (twin[1] - twin[0]) + twin[0], num=9)]),
	         'FRAC_DEP_LIM': '30.0',
	         'NLBGFLR': '.T',
	         'NLNBI_TEST_OUT': '.F',
	         'DN0OUT': '5.E10       ! external neutral density in cm ** -3',
	         'NPTCLS': '10000      ! no.of MC ptcls',
	         'NMCURB': '3      ! Use Lin - Liu & Hinton model',
	         'NPRAD': '2',
	         'NLBEAM': '.T',
	         'DTBEAM': '0.0001',
	         'NDEP0': '5000        !DE no.of deposition tracks',
	         'GOOMAX': '5000.0     ! max goosing factor, DE uncommented',
	         'REDGE': '4.126         !aperture dimensions, DE uncommented "halfwidth of aperture" for beam exit',
	         'NBSHAP': '2          ! 1=rectangular source grid, 2=circular source grid',
	         'FOCLR': '180.        ! focal lengths DE " ** R" is horizontal focal length',
	         'FOCLZ': '180.        !DE this is the vertical focal length',
	         'DIVR': '.02          ! divergences DE this is the horizontal divergence in radians',
	         'DIVZ': '.02        !DE vertical divergence of beam in radians',
	         'BMWIDR': '.19        ! source grid dimensions DE this is the half width',
	         'BMWIDZ': '.19        !DE this is the half height of the source grid',
	         'NCIRLM': '0          ! circle limiters',
	         'CRLMR1': '0.',
	         'CRLMY1': '0.',
	         'CRLMRD': '0.',
	         'NLINLM': '14      ! line limiters',
	         'ALNLMR': '18.52, 18.52, 27.9, 57.1, 104.33, 131.91, 149.51, 157.1, 149.51, 131.91, 104.33, 57.1, 27.9, 18.52',
	         'ALNLMY': '0.0, 100.81, 160.34, 160.34, 146.01, 103.97, 54.5, 0.0, -54.5, -103.97, -146.01, -160.34, -160.34, -100.81',
	         'ALNLMT': '90.0, -81.0, 0.0, 16.89, 56.73, 71.5, 81.04, 90.0, -81.04, -71.5, -56.73, -16.89, 0.0, 81.00',
	         'NBEAM': '1',
	         'NLCO(1)': '.T         !This means the beam is co-injected with the plasma current',
	         'EINJA(1)': f'{round(eb)}  ! the beam energy in eV',
	         'PINJA(1)': f'{round(pb)}  ! The beam power in Watts',
	         'TBONA(1)': f'{round(twin[0], 5)} !The beam turn on time',
	         'TBOFFA(1)': f'{round(twin[1], 5)}',
	         'ABEAMA(1)': '1.         !The atomic weight of the injected species (ie Hydrogen)',
	         'XZBEAMA(1)': '1.         !The charge of the injected species (seems like it would be 0?)',
	         'PDELTA(1)': '0.001    !Has no meaning, retained for compatibility',
	         'RTCENA(1)': '21.3     !Beam tangency radius in centimeters',
	         'XLBAPA(1)': '189.1    !Distance of ion source to beam aperture',
	         'XLBTNA(1)': '257.0    !The distance between the bean ion source and the tangecy radius',
	         'FFULLA(1)': '0.8     !The percent fraction of full energy particles accelerated ie H ^ +',
	         'FHALFA(1)': '0.1     !The percent fraction of half energy perticles accelerated ie H_2 ^ + !The percentage of the 1 / 3rd energy is just 1-FFULA-FHALF',
	         'XYBAPA(1)': '0.0         !Elevation of beam centerline above / below plasma midplane',
	         'XYBSCA(1)': '0.0         !Elevation of beam ion source above / below the plasma midplane',
	         'XBZETA(1)': '0.0         !Toroidal angle of the beam source in (R, zeta, Z) coordinates',
	         'nbi_pserve': '1         !transp told me to do this.This is not from Deyon Liu',
	         'EBDMAX': '2.50E+04 ! WJC max beam energy, sets min to 2 / 3 of 1 % of EBDMAX'}
	return bhead, bdict


def update_with_nbi(file, updict, header):
	print(f'Adding NBI to: {file}')
	used_keys = []
	new_content = header
	readin = open(file, 'r')
	for line in readin:
		for k, v in updict.items():
			if line.split('=')[0].strip() == k:
				line = f'{k} = {v}  !set by robot\n'
				used_keys.append(k)
		new_content += line
	for k in updict.keys():
		if k not in used_keys:
			print(f"Adding '{k}' to {file}")
			line = f'{k} = {updict[k]} ! added by robot\n'
			new_content += line
		ntimes = sum([1 for used in used_keys if used == k])  # count up num times key is found in TR.DAT
		if ntimes > 1:
			warnings.warn(f'Found multiple lines in {file} for {k}')
	readin.close()
	
	writeto = open(file, 'w')
	writeto.write(new_content)
	writeto.close()


def ltx_nubeamify(dir, trdat, **kwargs):
	# output_file_letter, ltx_shots, twin, ebeam, pbeam, ibeam
	if len(trdat) != 15:
		print(f'TRDAT filename incorrect format. Got {trdat} but expected 105795A01TR.DAT')
		return
	if 'ltx_shots' not in kwargs.keys():
		if not 'twin' in kwargs.keys():
			print('missing info: include either ltx_shots or (twin and two of three- ebeam, pbeam, ibeam)')
			return
		elif not (
				'ebeam' and 'pbeam' in kwargs.keys() or 'ebeam' and 'ibeam' in kwargs.keys() or 'pbeam' and 'ibeam' in kwargs.keys()):
			print('missing info: include either ltx_shots or (twin and two of three- ebeam, pbeam, ibeam)')
			return
	
	# ensure ebeam, pbeam defined
	if 'ltx_shots' in kwargs.keys():
		ep_inj = get_nbi_inputs_for_trdat(kwargs['ltx_shots'])
		if len(ep_inj) != 4:
			return
		else:
			kwargs['ebeam'], kwargs['pbeam'] = ep_inj[0], ep_inj[1]
		if 'twin' not in kwargs.keys():  # if not set by user, use ton/toff from nbi signals
			kwargs['twin'] = [ep_inj[2], ep_inj[3]]
	if 'ibeam' in kwargs.keys():
		kwargs['pbeam'] = kwargs['ebeam'] * kwargs['ibeam']
	
	# ensure output_file_letter defined
	if 'output_file_letter' not in kwargs.keys():
		ofl = 'A'
		while os.path.isfile(f'{dir}{trdat[0:6]}{ofl}{trdat[7:]}'):
			if ofl == 'Z':
				print('No output file letter available')
				return
			ofl = chr(ord(ofl) + 1)  # increment letter
		trdat_out = f'{trdat[0:6]}{ofl}{trdat[7:]}'
		resp = input(f'No output file letter given; output file will be- {trdat_out}, is this ok? (Y/N)\n')
		if resp.lower() != 'y':
			print('User denied output_file_letter, terminating.')
			return
	else:
		trdat_out = f"{trdat[0:6]}{kwargs['output_file_letter'].upper()}{trdat[7:]}"
		if os.path.isfile(f'{dir}{trdat_out}'):
			yn = input(f'output file {trdat_out} already exists: overwrite? (y/n)\n')
			if yn.lower() != 'y':
				print("fine I won't")
				return
	
	# create new file
	copy_content(f'{dir}{trdat}', f'{dir}{trdat_out}')
	
	bhead, bdict = get_beam_dict(kwargs['ebeam'], kwargs['pbeam'], kwargs['twin'])
	update_with_nbi(f'{dir}{trdat_out}', bdict, bhead)
	print(f'\nTR.DAT FILE WRITTEN:\n{dir}{trdat_out} is ready to be submitted...\n')


if __name__ == '__main__':
	dir = 'Z:/transp/t105952/'
	trdat = '105952A01TR.DAT'
	ltx_nubeamify(dir, trdat, ltx_shots=[105952], output_file_letter='c')

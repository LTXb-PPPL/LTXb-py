import io
import os
import stat
import warnings

import matplotlib.pyplot as plt
import numpy as np


# SOOOOO tired of going through countless files resetting stuff

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


def update_trdat(files, updict=None, comment=[], uncomment=[], remove=[]):
	for file in files:
		print(f'updating: {file}')
		used_keys = []
		new_content = ''
		readin = open(file, 'r')
		for line in readin:
			for k, v in updict.items():
				if line.split('=')[0].strip() == k:
					line = f'{k} = {v}  !set by robot\n'
					used_keys.append(k)
			for k in comment:
				if line.strip().startswith(k):
					line = f'!{line}'
			for k in uncomment:
				if line.startswith(f'!{k}'):
					line = line.replace('!', '')
			if not any([rmv in line for rmv in remove]):
				new_content += line
		for k in updict.keys():
			if k not in used_keys:
				print(f"Adding '{k}' to {file}")
				line = f'{k} = {updict[k]} ! added by robot\n'
				new_content = line + new_content  # add it at the TOP so it occurs before &END in TR.DAT
			ntimes = sum([1 for used in used_keys if used == k])  # count up num times key is found in TR.DAT
			if ntimes > 1:
				warnings.warn(f'Found multiple lines in {file} for {k}')
		
		readin.close()
		
		writeto = open(file, 'w')
		writeto.write(new_content)
		writeto.close()


def create_new(copyfrom=None, newnames=None, setdict=None):
	filelist = [copyfrom[:-9] + nn + 'TR.DAT' for nn in newnames]
	go = 1
	for file in filelist:
		if os.path.exists(file):
			yn, go = '', 0
			print(f'file {file} already exists')
	if not go:
		while yn != 'y' and yn != 'n':
			yn = input('OVERWRITE ALL? (y/n)')
		if yn == 'y':
			go = 1
	if go:
		for i, file in enumerate(filelist):
			copy_content(copyfrom, file)  # create new files copied from [copyfrom]
			update = {}
			for k, v in setdict.items():
				update[k] = v[i]
			update_trdat([file], updict=update)
	write_sh_file2(filelist)


def write_sh_file(filelist):
	nums = [int(fn[-8:-6]) for fn in filelist]
	if not np.unique([nums[i + 1] - nums[i] for i in np.arange(len(nums) - 1)])[0] == 1:
		print('non-sequential TR.DAT filenames found, write your own shell script. returning...')
		return
	else:
		dir = '/'.join(filelist[0].split('/')[:-1])
		tofile = f'{dir}/run_transp_submissions.sh'
		nprocs = np.floor(256 / len(filelist))
		shotletter = filelist[0].split('/')[-1][-15:-8]
		content = f'#!/bin/sh\n#created by robot to run TR.DAT files\n' \
		          f'echo "256 core user limit, only submit using <={int(nprocs)} CORES"\n' \
		          f'for i in $(seq -f "%02g" {min(nums)} {max(nums)})\ndo\n' \
		          '	echo " "\n	echo " "\n	echo " "\n' \
		          f'	echo "starting {shotletter}$i"\n' \
		          f'	tr_start {shotletter}$i\n' \
		          f'	echo "submitting {shotletter}$i"\n' \
		          f'	tr_send {shotletter}$i\n' \
		          'done\necho " "\necho " "\necho " "\n' \
		          'echo "SUBMISSION COMPLETE"'
		writeto = io.open(tofile, 'w', newline='\n')
		writeto.write(content)
		writeto.close()
		os.chmod(tofile, stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
		# os.chmod(tofile, 0o777)
		print(f'wrote file: {tofile}')


def write_sh_file2(filelist):
	nums = [int(fn[-8:-6]) for fn in filelist]
	if not np.unique([nums[i + 1] - nums[i] for i in np.arange(len(nums) - 1)])[0] == 1:
		print('non-sequential TR.DAT filenames found, write your own shell script. returning...')
		return
	else:
		dir = '/'.join(filelist[0].split('/')[:-1])
		tofile = f'{dir}/run_transp_submissions.sh'
		nprocs = np.floor(256 / len(filelist))
		shotletter = filelist[0].split('/')[-1][-15:-8]
		
		# f'echo "256 core user limit, only submit using <={int(nprocs)} CORES"\n' \
		content = f'#!/usr/bin/expect --\n#created by robot to run TR.DAT files\n' \
		          f'for i in $(seq -f "%02g" {min(nums)} {max(nums)})\ndo\n' \
		          '	echo " "\n	echo " "\n	echo " "\n' \
		          f'	echo "starting {shotletter}$i"\n' \
		          f'	spawn tr_start {shotletter}$i\n' \
		          '	expect "tokyy_id"\n	send "ltx"\r\n	expect "tokYY_id"\n	send "y"\r\n' \
		          '	expect "Enter Comments for TR.INF File:"\n	send "x"\r\n' \
		          '	expect "verify your email address = wcapecch@pppl.gov (Y/N):"\n	send "y"\r\n' \
		          '	expect eof\n' \
		          f'	echo "submitting {shotletter}$i"\n' \
		          f'	tr_send {shotletter}$i\n' \
		          'done\necho " "\necho " "\necho " "\n' \
		          'echo "SUBMISSION COMPLETE"'
		writeto = io.open(tofile, 'w', newline='\n')
		writeto.write(content)
		writeto.close()
		os.chmod(tofile, stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
		# os.chmod(tofile, 0o777)
		print(f'wrote file: {tofile}')


def ip_n0_eb_scan(cleancopy=False):
	# need to standardize all the runs here
	# ip_arr = [81, 100, 121, 138, 155]  # kA
	ip_scan = [f'Z:/transp/t100002/100002C{str(i).zfill(2)}TR.dat' for i in [1, 2, 3, 4, 5]]
	ip_arr = ["'HANS81KA'", "'HANS100KA'", "'HANS121KA'", "'HANS138KA'", "'HANS155KA'"]
	
	n0_scan_155ka = [f'Z:/transp/t100002/100002C{str(i).zfill(2)}TR.dat' for i in [6, 7, 8, 9, 10, 11, 12, 13, 14, 15]]
	n0_scan_81ka = [f'Z:/transp/t100002/100002C{str(i).zfill(2)}TR.dat' for i in
	                [16, 17, 18, 19, 20, 21, 22, 23, 24, 25]]
	# n0_arr = [1.e18, 2.3e18, 5.e18, 1.e19, 2.3e19, 5.e19, 1.e20, 2.3e20, 5.e20, 8.e20]
	n0_arr = ["'CAP1e18'", "'CAP2P3e18'", "'CAP5e18'", "'CAP1e19'", "'CAP2P3e19'", "'CAP5e19'", "'CAP1e20'",
	          "'CAP2P3e20'", "'CAP5e20'", "'CAP8e20'"]
	
	eb_scan_155ka = [f'Z:/transp/t100002/100002C{str(i).zfill(2)}TR.dat' for i in [26, 27, 28, 29, 30, 31]]
	eb_scan_81ka = [f'Z:/transp/t100002/100002C{str(i).zfill(2)}TR.dat' for i in [32, 33, 34, 35, 36, 37]]
	# e_nbi = [10, 12, 14, 16, 18, 20]  # keV
	eb_arr = ['1.0d4', '1.2d4', '1.4d4', '1.6d4', '1.8d4', '2.0d4']
	
	optim = [f'Z:/transp/t100002/100002C{str(i).zfill(2)}TR.dat' for i in [38, 39]]
	opt_arr = ['1.0d4', '2.0d4']
	# start with template from good_template_TR.DAT - new limiter, 3us beam, good timesteps, etc.
	
	if cleancopy:
		copy_from = f'Z:/transp/t100002/good_template_TR.dat'
		for filelist in [ip_scan, n0_scan_155ka, n0_scan_81ka, eb_scan_155ka, eb_scan_81ka, optim]:
			for file in filelist:
				copy_content(copy_from, file)
	
	nominal_n0 = "'CAP1e19'"
	nominal_eb = '2.0d4'
	for i, file in enumerate(ip_scan):
		eq = ip_arr[i]
		update = {'PREGRB': eq, 'PREPRS': eq, 'PRETRF': eq, 'PREPLF': eq, 'PREQPR': eq, 'PREMMX': eq, 'PRECUR': eq,
		          'PRELIM': eq, 'PRENER': nominal_n0, 'EINJA(1)': nominal_eb}
		update_trdat([file], updict=update)
	for i, file in enumerate(n0_scan_155ka):
		eq = "'HANS155KA'"
		update = {'PREGRB': eq, 'PREPRS': eq, 'PRETRF': eq, 'PREPLF': eq, 'PREQPR': eq, 'PREMMX': eq, 'PRECUR': eq,
		          'PRELIM': eq, 'PRENER': n0_arr[i], 'EINJA(1)': nominal_eb}
		update_trdat([file], updict=update)
	for i, file in enumerate(n0_scan_81ka):
		eq = "'HANS81KA'"
		update = {'PREGRB': eq, 'PREPRS': eq, 'PRETRF': eq, 'PREPLF': eq, 'PREQPR': eq, 'PREMMX': eq, 'PRECUR': eq,
		          'PRELIM': eq, 'PRENER': n0_arr[i], 'EINJA(1)': nominal_eb}
		update_trdat([file], updict=update)
	for i, file in enumerate(eb_scan_155ka):
		eq = "'HANS155KA'"
		update = {'PREGRB': eq, 'PREPRS': eq, 'PRETRF': eq, 'PREPLF': eq, 'PREQPR': eq, 'PREMMX': eq, 'PRECUR': eq,
		          'PRELIM': eq, 'PRENER': nominal_n0, 'EINJA(1)': eb_arr[i]}
		update_trdat([file], updict=update)
	for i, file in enumerate(eb_scan_81ka):
		eq = "'HANS81KA'"
		update = {'PREGRB': eq, 'PREPRS': eq, 'PRETRF': eq, 'PREPLF': eq, 'PREQPR': eq, 'PREMMX': eq, 'PRECUR': eq,
		          'PRELIM': eq, 'PRENER': nominal_n0, 'EINJA(1)': eb_arr[i]}
		update_trdat([file], updict=update)
	for i, file in enumerate(optim):
		eq = "'HANS155KA'"
		update = {'PREGRB': eq, 'PREPRS': eq, 'PRETRF': eq, 'PREPLF': eq, 'PREQPR': eq, 'PREMMX': eq, 'PRECUR': eq,
		          'PRELIM': eq, 'PRENER': "'CAP5e19'", 'EINJA(1)': opt_arr[i]}
		update_trdat([file], updict=update)


if __name__ == '__main__':
	# create tangency scans for 106536R02 at all 4 orientations (normal ops, rev-Ip, rev-Bt, rev-Ip*Bt)
	rtans = [19., 28., 36., 44., 53.]
	nlbccw = ['.T', '.T', '.F', '.F']
	nljccw = ['.F', '.T', '.F', '.T']
	nlco1 = ['.T', '.F', '.T', '.F']
	_, bb = np.meshgrid(rtans, nlbccw)
	_, jj = np.meshgrid(rtans, nljccw)
	rr, cc = np.meshgrid(rtans, nlco1)
	a=1
	# create_new('Z:/transp/t106536/106536R02TR.DAT', [f'T{i + 1:02}' for i in np.arange(len(rr.flatten()))],
	#            {'RTCENA(1)': [f'{rt}' for rt in rr.flatten()], 'NLCO(1)': cc.flatten(), 'NLBCCW': bb.flatten(),
	#             'NLJCCW': jj.flatten()})

# focal length scan (180 nominal). range 160, 165, ... 200
# foc = [str(f) for f in np.linspace(160, 200, endpoint=True, num=9)]
# create_new('Z:/transp/t103617/103617C01TR.DAT', [f'F{i:02}' for i in np.arange(1, 10)],
#            {'FOCLR': foc, 'FOCLZ': foc})

# divergence scan (.02 nominal) range .02,.03, ... .10
# div = [str(d) for d in np.linspace(.02, .1, endpoint=True, num=9)]
# create_new('Z:/transp/t103617/103617C01TR.DAT', [f'D{i:02}' for i in np.arange(1, 10)], {'DIVR': div, 'DIVZ': div})
# a = 1

# create_new('Z:/transp/t103617/103617C01TR.DAT',
#            ['V01', 'V02', 'V03', 'V04', 'V05', 'V06', 'V07', 'V08', 'V09', 'V10', 'V11'], {
# 	           'EINJA(1)': ['1.0d4', '1.1d4', '1.2d4', '1.3d4', '1.4d4', '1.5d4', '1.6d4', '1.7d4', '1.8d4',
# 	                        '1.9d4', '2.0d4']})
# create_new('Z:/transp/t103617/103617C01TR.DAT',
#            ['N01', 'N02', 'N03', 'N04', 'N05', 'N06', 'N07', 'N08', 'N09', 'N10', 'N11'], {
# 	           'DN0OUT': ['1.e8', '5.e8', '1.e9', '5.e9', '1.e10', '5.e10', '1.e11', '5.e11', '1.e12', '5.e12',
# 	                      '1.e13']})
# perv = 15.e-6
# vb = np.arange(10, 21) * 1000.  # 10-20kV
# ib = perv * vb ** 1.5  # A
# pb = ib * vb  # W
# vbarr = [f'{v:1.1e}' for v in vb]
# pbarr = [f'{p:1.1e}' for p in pb]
# newnames = [f'P{i:02}' for i in np.arange(1, 12)]
# create_new('Z:/transp/t103617/103617C01TR.DAT', newnames, {'EINJA(1)': vbarr, 'PINJA(1)': pbarr})

# ip_n0_eb_scan(cleancopy=True)

# fix_stuff_with_current_density_energy_scan = False
# if fix_stuff_with_current_density_energy_scan:
# 	files = [f'Z:/transp/t100002/100002D{str(i).zfill(2)}TR.dat' for i in np.arange(15) + 1]
# 	update = {'EINJA(1)': '1.6d4'}
# 	update_trdat(files, updict=update)

# create tangency scans for 105795C02 and 105952C02
# rtans = [18.0, 20.0, 22., 24., 26., 28., 30., 32., 34., 36., 38., 40.]
# create_new('Z:/transp/t105795/105795C02TR.DAT', [f'T{i + 1:02}' for i in np.arange(len(rtans))],
#            {'RTCENA(1)': [f'{rt}' for rt in rtans]})
# create_new('Z:/transp/t105952/105952C02TR.DAT', [f'T{i + 1:02}' for i in np.arange(len(rtans))],
#            {'RTCENA(1)': [f'{rt}' for rt in rtans]})

# neutral density scan at various energies 3Sept22
# n0_arr = ['1.E8', '1.E9', '1.E10', '1.E11', '1.E12']  # num/cm^3
# ebeam_arr = ['10.E3', '13.E3', '16.E3']  # eV
# nn, ee = np.meshgrid(n0_arr, ebeam_arr)
# create_new('Z:/transp/t106536/106536R01TR.DAT', [f'N{i+1:02}' for i in np.arange(len(nn.flatten()))],
#            {'EINJA(1)': ee.flatten(), 'DN0OUT': nn.flatten()})
# create_new('Z:/transp/t105795/105795R01TR.DAT', [f'N{i+1:02}' for i in np.arange(len(nn.flatten()))],
#            {'EINJA(1)': ee.flatten(), 'DN0OUT': nn.flatten()})

import subprocess
import os

cwd = os.getcwd()
reco_dir = 'Y:/reconstructions/runday/'
matlab_dir = 'Z:/matlab/conbeam/'
fig_dir = "Z:/matlab/botdata/"


def create_map(eqdsk, shot, runnum, ebeam, nel):
	os.chdir(matlab_dir)
	# eqdsk = eqdsk.replace('\\', '/')
	print(f'making map for: {eqdsk}')
	command_string = ''.join(
		f'shots={[shot]};Ecsts={[ebeam]};runnums={[runnum]};times={[-1]};ips={[-1]};nels={[nel]};bot_create_map;'.split(
			'\n'))
	msgs = subprocess.run(['matlab', '-batch', command_string], capture_output=True)
	# print('exist status:', msgs.returncode)
	# print(f'{msgs.stdout.decode()}')
	# print(f'{msgs.stderr.decode()}')
	os.chdir(cwd)

	
if __name__ == '__main__':
	pass

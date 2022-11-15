import io

from toolbox.helpful_stuff import *

nprocs = 5.2
askprocs=2
while askprocs <= nprocs:
	if 2*askprocs <= nprocs:
		askprocs *= 2
	else:
		break
content = '#!/usr/bin/expect --\n#created by robot to run TR.DAT files\n'
for i in np.arange(2):
	content += f'spawn tr_start 105795F{i + 1:02}\n' \
	           f'expect "Enter number of processes"\nsend "{askprocs}\\r"\n' \
	           'expect "tokyy_id"\nsend "ltx\\r"\nexpect "tokYY_id"\nsend "y\\r"\n' \
	           'expect "Enter Comments for TR.INF File:"\nsend "x\\r"\n' \
	           'expect "verify your email address = wcapecch@pppl.gov (Y/N):"\nsend "y\\r"\n' \
	           'expect eof\n' \
	           f'spawn tr_send 105795F{i + 1:02}\n'
tofile = 'Z:/transp/t105795/test_batch_submit.sh'
writeto = io.open(tofile, 'w', newline='\n')
writeto.write(content)
writeto.close()

if __name__ == '__main__':
	# shot105124()
	# gauss2d_integral()
	do_nothing = 1

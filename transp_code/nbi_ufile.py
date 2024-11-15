"""
This code is to replicate OMFit generation of Ufiles for transp runs with time-dependent NBI data
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import re
from mpl_toolkits.mplot3d import axes3d
import io
import datetime

matplotlib.use('TkAgg')  # allows plotting in debug mode
clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']


def read_nbi_ufile(fn=None, lbls=None):
    if fn is None:
        direc = 'C:/Users/willi/PycharmProjects/LTXb-py/transp_code/data/'
        fn = f'{direc}OMF176440.NBI'
        lbls = ['30LT', '30RT', '150LT_UP', '150_LT_MU', '150LT_ML', '150LT_LO', '150RT_UP', '150_RT_MU', '150RT_ML',
                '150RT_LO', '210LT', '210RT', '330LT', '330RT']

    f = open(fn, 'r')
    form1 = r'-?\d+.\d+E[+-]\d+'
    form2 = r'-?\d+.\d+e[+-]\d+'
    for i in np.arange(3):  # skip stuff in first 24 lines
        f.readline()
    nbeams = int(float((re.findall(form1, f.readline())[0])))
    for i in np.arange(5):  # skip a few more lines
        f.readline()
    nx0 = int(re.findall(r'\d+', f.readline())[0])
    nx1 = int(re.findall(r'\d+', f.readline())[0])
    x0, x1, mtrx = [], [], []
    form = None
    for i in np.arange(np.ceil(nx0 / 6.)):
        line = f.readline()
        if form is None:  # perform check to see which form to use
            if len(re.findall(form1, line)) > len(re.findall(form2, line)):
                form = form1
            else:
                form = form2
        x0.extend([float(n) for n in re.findall(form, line)])
    for j in np.arange(np.ceil(nx1 / 6.)):
        line = f.readline()
        x1.extend([float(n) for n in re.findall(form, line)])
    for k in np.arange(np.ceil(nx0 * nx1 / 6.)):
        line = f.readline()
        mtrx.extend([float(n) for n in re.findall(form, line)])
    f.close()

    xx0, xx1 = np.meshgrid(x0, x1)
    m = np.array(mtrx).reshape((nx1, nx0))
    # ax = plt.figure().add_subplot(projection='3d')
    # ax.plot_surface(xx0, xx1, m, edgecolor='royalblue', lw=0.5, rstride=8, cstride=8, alpha=0.3)
    # plt.show()

    fig, ((axa, axb), (axc, axd)) = plt.subplots(nrows=2, ncols=2, sharex=True)
    for nb in np.arange(nbeams):
        if nb > len(clrs):
            ls = '--'
        else:
            ls = '-'
        if lbls is not None:
            axa.plot(x0, m[nb, :], c=clrs[nb % len(clrs)], label=lbls[nb], ls=ls)
        else:
            axa.plot(x0, m[nb, :], c=clrs[nb % len(clrs)], ls=ls)
        axb.plot(x0, m[nb + nbeams, :], c=clrs[nb % len(clrs)], ls=ls)
        axc.plot(x0, m[nb + 2*nbeams, :], c=clrs[nb % len(clrs)], ls=ls)
        axd.plot(x0, m[nb + 3*nbeams, :], c=clrs[nb % len(clrs)], ls=ls)
    if lbls is not None:
        axa.legend()
    axa.set_ylabel('pbeam (W)')
    axb.set_ylabel('ebeam (V)')
    axc.set_ylabel('ffull')
    axd.set_ylabel('fhalf')
    plt.tight_layout()
    plt.show()


def write_nbi_ufile(shot=None):
    if shot is None:  # make up some data to write for testing
        t = np.linspace(0, 300, num=1000)
        pb = np.zeros_like(t)
        eb = np.zeros_like(t)
        itb = np.where((t >= 200) & (t < 210))
        pb[itb] = 700.e3  # 700 kW
        eb[itb] = 20.e3  # 20 keV
        ffull = np.ones_like(pb)*.8
        fhalf = np.ones_like(pb)*.1
        shot = 123456
        shotdate = '20July-1969'  # copying ufile form, just using random date of moon landing
    else:
        pass  # get data from mdsplus for this shot

    direc = 'C:/Users/willi/PycharmProjects/LTXb-py/transp_code/data/'
    tofile = f'{direc}testfile{shot}.NBI'
    nbeams = 1  # FOR NOW!! Will LTX ever see a 2nd beam installed? Stay tuned!!
    nt = len(t)
    nchan = 4 * nbeams  # not sure what each of the 4 are right now, power and energy are two...

    content = f'{shot}LTX 0 0 0             ; SHOT\n{shotdate}               ;SHOT DATE\n'
    content += '1                           ;Number of associated scalar quantities\n'
    content += f'{nbeams:.4E}                         ;Scalar, label follows:\n'  # format to match test Ufile
    content += f'NBEAM:  # of beams\nTime Seconds               ;Indep var label: X0\n'
    content += f'Channel No                 ;Indep var label X1\n'
    content += f'Mixed Beam Data            ;Dep var label\n'
    content += f'0                          ;Proc Code 0:Raw 1:Avg 2:Sm 3:Avg+Sm\n'
    content += f'{nt}                       ;# of X0 pts\n{nchan}                           ;# of X1 pts\n'

    # write time data
    for irow in np.arange(int(np.ceil(len(t) / 6))):  # 6 columns
        content += '  '.join([f'{val:.5e}' for val in t[irow * 6:irow * 6 + 6]])+'\n'
    # write channel array
    for irow in np.arange(int(np.ceil(nchan / 6))):  # single row for us with only 1 beam, but leave it general
        content += '  '.join([f'{val+1:.5e}' for val in np.arange(nchan)[irow * 6:irow * 6 + 6]])+'\n'
    # write data arrays
    data = np.append(pb, np.append(eb, np.append(ffull, fhalf)))
    for irow in np.arange(int(np.ceil(len(data) / 6))):
        content += '  '.join([f'{val:.5e}' for val in data[irow * 6:irow * 6 + 6]])+'\n'

    content += ';-------END OF DATA-------------\n'
    content += f'Written by python routine on {datetime.datetime.today().strftime("%d %b %Y")}'
    content += f'Author: Bill Capecchi wcapecch@pppl.gov'

    writeto = io.open(tofile, 'w', newline='\n')
    writeto.write(content)
    writeto.close()


if __name__ == '__main__':
    # read_nbi_ufile()
    # write_nbi_ufile()
    direc = 'C:/Users/willi/PycharmProjects/LTXb-py/transp_code/data/'
    tofile = f'{direc}testfile{123456}.NBI'
    read_nbi_ufile(tofile)

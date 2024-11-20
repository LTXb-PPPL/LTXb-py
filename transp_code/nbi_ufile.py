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
from bills_LTX_MDSplus_toolbox import get_tree_conn, get_data
import os

matplotlib.use('TkAgg')  # allows plotting in debug mode
clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']


def retrieve_nbi_shot_data(shot, plot=False):
    if shot > 200000:
        onlynbi = True
    else:
        onlynbi = False
    if onlynbi:
        print(f'shot {shot} is NBI only, not suitable for TRANSP run')
        return
    tree = 'ltx_b'
    prefix = '.oper_diags.ltx_nbi'
    t = get_tree_conn(shot, treename=tree)
    ts = get_data(t, '.metadata:timestamp')
    print(f'gathering data for shot {shot} occurred on {ts}')
    (tibeam, ibeam) = get_data(t, f'{prefix}.source_diags.i_hvps')
    (tbeam, vbeam) = get_data(t, f'{prefix}.source_diags.v_hvps')
    ibeam = np.interp(tbeam, tibeam, ibeam)
    pbeam = ibeam * vbeam / 1000.  # kW
    if plot:
        fig, ax = plt.subplots()
        axr = ax.twinx()
        ax.plot(tbeam, vbeam)
        axr.plot(tbeam, pbeam)
        ax.set_ylabel('ebeam (V)')
        axr.set_ylabel('pbeam (W)')
        ax.set_xlabel('time (s)')
    return tbeam, pbeam, vbeam, ts


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
        axc.plot(x0, m[nb + 2 * nbeams, :], c=clrs[nb % len(clrs)], ls=ls)
        axd.plot(x0, m[nb + 3 * nbeams, :], c=clrs[nb % len(clrs)], ls=ls)
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
        shot = 123456
        shotdate = '20July-1969'  # copying ufile form, just using random date of moon landing
    else:
        t, pb, eb, shotdate = retrieve_nbi_shot_data(shot, plot=False)

    # pb[np.where(pb < 0)] = 0.
    # eb[np.where(eb < 0)] = 0.

    # until we know better, set full/half fractions here
    ffull = np.ones_like(pb) * .8
    fhalf = np.ones_like(pb) * .1

    direc = f'{os.getcwd()}/data/'
    tofile = f'{direc}AVN{shot}.NBI'  # "time varying neutral beam"
    nbeams = 1  # FOR NOW!! Will LTX ever see a 2nd beam installed? Stay tuned!!
    nt = len(t)
    nchan = 4 * nbeams  # pbeam, ebeam, ffull, fhalf

    content = f'{shot}LTX 2 0 6             ; SHOT\n{shotdate}               ;SHOT DATE\n'
    content += '1                           ;Number of associated scalar quantities\n'
    content += f'{nbeams:.4E}                         ;Scalar, label follows:\n'  # format to match test Ufile
    content += f'NBEAM:  # of beams\n'
    content += 'TIME                SECONDS ;-INDEPENDENT VARIABLE LABEL: X-\n'
    content += f'CHANNEL NO                    ;-INDEPENDENT VARIABLE LABEL: Y-\n'
    content += f'MIXED BEAM DATA  MIXED        ;-DEPENDENT VARIABLE LABEL-\n'
    content += f'0                             ;-PROC CODE- 0:RAW 1:AVG 2:SM 3:AVG+SM\n'
    content += f'        {nt}                    ;-# OF X PTS-\n'
    content += f'          {nchan}                    ;-# OF Y PTS- X,Y,F(X,Y) DATA FOLLOW:\n'

    # write time data
    for irow in np.arange(int(np.ceil(len(t) / 6))):  # 6 columns
        content += '  '.join([f'{val:.5e}' for val in t[irow * 6:irow * 6 + 6]]) + '\n'
    # write channel array
    for irow in np.arange(int(np.ceil(nchan / 6))):  # single row for us with only 1 beam, but leave it general
        content += '  '.join([f'{val + 1:.5e}' for val in np.arange(nchan)[irow * 6:irow * 6 + 6]]) + '\n'
    # write data arrays
    data = np.append(pb, np.append(eb, np.append(ffull, fhalf)))
    for irow in np.arange(int(np.ceil(len(data) / 6))):
        content += '  '.join([f'{val:.5e}' for val in data[irow * 6:irow * 6 + 6]]) + '\n'

    content += ';-------END OF DATA-------------\n'
    content += f'Written by python routine on {datetime.datetime.today().strftime("%d %b %Y")}'
    content += f'Author: Bill Capecchi wcapecch@pppl.gov'

    writeto = io.open(tofile, 'w', newline='\n')
    writeto.write(content)
    writeto.close()
    print(f'wrote {tofile}')


if __name__ == '__main__':
    # read_nbi_ufile()
    # write_nbi_ufile()
    # direc = 'C:/Users/willi/PycharmProjects/LTXb-py/transp_code/data/'
    # tofile = f'{direc}testfile{123456}.NBI'
    # read_nbi_ufile(tofile)
    write_nbi_ufile(108890)
    # shot = 108890
    # read_nbi_ufile(f'{os.getcwd()}/data/TVN{shot}.NBI')

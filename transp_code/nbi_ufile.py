"""
This code is to replicate OMFit generation of Ufiles for transp runs with time-dependent NBI data
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import re
import io
import datetime
import os

if os.getcwd().startswith('Z:'):  # working on laptop
    from bills_LTX_MDSplus_toolbox import get_tree_conn, get_data

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


def read_nbi_ufile(shot, pre):  # fn=None, lbls=None):
    direc = f'Z:/transp/t{shot}/'
    fn = f'{direc}{pre}{shot}.NBI'  # "time varying neutral beam"
    fnshrt = f'{pre}{shot}.NBI'

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

    m = np.array(mtrx).reshape((nx1, nx0))

    fig, ((axa, axb), (axc, axd)) = plt.subplots(nrows=2, ncols=2, sharex=True)
    for nb in np.arange(nbeams):
        if nb > len(clrs):
            ls = '--'
        else:
            ls = '-'
        axa.plot(x0, m[nb, :], c=clrs[nb % len(clrs)], ls=ls)
        axb.plot(x0, m[nb + nbeams, :], c=clrs[nb % len(clrs)], ls=ls)
        axc.plot(x0, m[nb + 2 * nbeams, :], c=clrs[nb % len(clrs)], ls=ls)
        axd.plot(x0, m[nb + 3 * nbeams, :], c=clrs[nb % len(clrs)], ls=ls)
    axa.set_ylabel('pbeam (W)')
    axb.set_ylabel('ebeam (V)')
    axc.set_ylabel('ffull')
    axd.set_ylabel('fhalf')
    plt.suptitle(fnshrt)
    plt.tight_layout()
    plt.show()


def get_trdat_timing(trdat=None):
    if trdat is None:
        trdat = 'C:/Users/willi/Dropbox/work_stuff/108890V01TR.DAT'
    tstrt, tend = None, None
    f = open(trdat, 'r')
    lines = f.readlines()
    for l in lines:
        if 'TINIT' in l:
            tstrt = float(l.split('=')[-1].split(' ')[0])
        if 'FTIME' in l and 'FTIME_OK' not in l:
            tend = float(l.split('=')[-1].split(' ')[0])
    f.close()
    if tstrt is None:
        print(f'no start time found in {trdat}')
        return None, None
    if tend is None:
        print(f'no end time found in {trdat}')
        return None, None
    return tstrt, tend


def write_nbi_ufile(shot=None, alphanum='A01'):
    # assume working from work laptop (so we can get beam data)
    direc = f'Z:/transp/t{shot}/'
    trdat = f'{direc}{shot}{alphanum}TR.DAT'
    tstrt, tend = get_trdat_timing(trdat)
    tdat, pdat, edat, shotdate = retrieve_nbi_shot_data(shot, plot=False)
    t = np.arange(tstrt, tend, tdat[1] - tdat[0])  # new time axis spanning sim time
    pb, eb = np.interp(t, tdat, pdat), np.interp(t, tdat, edat)

    pb[np.where(pb < 0)] = 0.
    eb[np.where(eb < 0)] = 0.

    # until we know better, set full/half fractions here
    ffull = np.ones_like(pb) * .8
    fhalf = np.ones_like(pb) * .1

    tofile = f'{direc}PVN{shot}.NBI'  # "time varying neutral beam"
    nbeams = 1  # FOR NOW!! Will LTX ever see a 2nd beam installed? Stay tuned!!
    nt = len(t)
    nchan = 4 * nbeams  # pbeam, ebeam, ffull, fhalf

    content = f' {shot}LTX  2 0 6             ; SHOT\n'
    content += f' {shotdate}               ;SHOT DATE\n'
    content += '   1                           ;Number of associated scalar quantities\n'
    content += f'  {nbeams:.4E}                         ;Scalar, label follows:\n'  # format to match test Ufile
    content += f' NBEAM:    # OF BEAMS          \n'
    content += ' TIME                SECONDS ;-INDEPENDENT VARIABLE LABEL: X-\n'
    content += f' CHANNEL NO                    ;-INDEPENDENT VARIABLE LABEL: Y-\n'
    content += f' MIXED BEAM DATA  MIXED        ;-DEPENDENT VARIABLE LABEL-\n'
    content += f' 0                             ;-PROC CODE- 0:RAW 1:AVG 2:SM 3:AVG+SM\n'
    content += f'       {nt}                    ;-# OF X PTS-\n'
    content += f'          {nchan}                    ;-# OF Y PTS- X,Y,F(X,Y) DATA FOLLOW:\n'

    # write time data
    for irow in np.arange(int(np.ceil(len(t) / 6))):  # 6 columns
        content += '  ' + '  '.join([f'{val:.5e}' for val in t[irow * 6:irow * 6 + 6]]) + '\n'
    # write channel array
    for irow in np.arange(int(np.ceil(nchan / 6))):  # single row for us with only 1 beam, but leave it general
        content += '  ' + '  '.join([f'{val + 1:.5e}' for val in np.arange(nchan)[irow * 6:irow * 6 + 6]]) + '\n'
    # write data arrays
    data = np.append(pb, np.append(eb, np.append(ffull, fhalf)))
    for irow in np.arange(int(np.ceil(len(data) / 6))):
        content += '  ' + '  '.join([f'{val:.5e}' for val in data[irow * 6:irow * 6 + 6]]) + '\n'

    content += ';-------END OF DATA-------------\n'
    content += f'Written by python routine on {datetime.datetime.today().strftime("%d %b %Y")}\n'
    content += f'Author: Bill Capecchi wcapecch@pppl.gov'

    writeto = io.open(tofile, 'w', newline='\n')
    writeto.write(content)
    writeto.close()
    print(f'wrote {tofile}')


def average_multiple_shots(ltx_shots):
    einj, pinj = np.array([]), np.array([])
    pad = 0.25e-3  # distance to step twin away from beam turnon/turnoff
    ton, toff = np.array([]), np.array([])
    for shot in ltx_shots:
        try:
            t = get_tree_conn(shot, treename='ltx_b')
            ts = get_data(t, '.metadata:timestamp')
            # print('gathering data for shot {} occurred on {}'.format(shot, get_data(t, '.metadata:timestamp')))
            (tibeam, ibeam) = get_data(t, '.oper_diags.ltx_nbi.source_diags.i_hvps')
            (tbeam, vbeam) = get_data(t, '.oper_diags.ltx_nbi.source_diags.v_hvps')
            ibeam = np.interp(tbeam, tibeam, ibeam)
            tbeamon, tbeamoff = tbeam[np.where(vbeam > 1000)][0] + pad, tbeam[np.where(vbeam > 1000)][-1] - pad
            ton, toff = np.append(ton, tbeamon), np.append(toff, tbeamoff)
            ii = np.where((tbeam >= tbeamon + pad) & (tbeam <= tbeamoff - pad))
            einj = np.append(einj, np.mean(vbeam[ii]))
            pinj = np.append(pinj, np.mean(vbeam[ii] * ibeam[ii]))
            print(f'shot {shot} <Einj> = {np.mean(einj) / 1000.:.2f} [keV] occurred {ts}')
        except:
            print(f'problem occurred gathering data for shot {shot}: returning')
            return [np.nan]
    checkshots = False
    if max(ton) - min(ton) > 1.e-3:
        print('discrepancy in beam turnon times, check shots')
        checkshots = True
    transmission_efficiency = 0.9  # fraction not lost along beam path into vessel
    neutralization_efficiency = 0.8  # fraction neutralized in neutralizer
    power_into_vessel = np.mean(pinj) * transmission_efficiency * neutralization_efficiency

    tbeam, pinj, einj, shotdate = 1, 1, 1, 1
    return tbeam, pinj, einj, shotdate


def write_nbi_ufile_multishot(direc=None, trshot=None, tralphanum=None, shots=None):
    # MUST work from work laptop (so we can get beam data)
    if trshot is None:
        trshot = 108890
    if direc is None:
        direc = f'Z:/transp/t{trshot}/'
    if tralphanum is None:
        tralphanum = 'C01'
    if shots is None:
        shots = [108886, 108887, 108890, 108880, 108874, 108881, 108895, 108891, 108876, 108875, 108894]
    trdat = f'{direc}{trshot}{tralphanum}TR.DAT'

    # tstrt, tend = get_trdat_timing(trdat)
    tdat, pdat, edat, shotdate = average_multiple_shots(shots)

    tdat, pdat, edat, shotdate = retrieve_nbi_shot_data(shot, plot=False)
    t = np.arange(tstrt, tend, tdat[1] - tdat[0])  # new time axis spanning sim time
    pb, eb = np.interp(t, tdat, pdat), np.interp(t, tdat, edat)

    pb[np.where(pb < 0)] = 0.
    eb[np.where(eb < 0)] = 0.

    # until we know better, set full/half fractions here
    ffull = np.ones_like(pb) * .8
    fhalf = np.ones_like(pb) * .1

    tofile = f'{direc}PVN{shot}.NBI'  # "time varying neutral beam"
    nbeams = 1  # FOR NOW!! Will LTX ever see a 2nd beam installed? Stay tuned!!
    nt = len(t)
    nchan = 4 * nbeams  # pbeam, ebeam, ffull, fhalf

    content = f' {shot}LTX  2 0 6             ; SHOT\n'
    content += f' {shotdate}               ;SHOT DATE\n'
    content += '   1                           ;Number of associated scalar quantities\n'
    content += f'  {nbeams:.4E}                         ;Scalar, label follows:\n'  # format to match test Ufile
    content += f' NBEAM:    # OF BEAMS          \n'
    content += ' TIME                SECONDS ;-INDEPENDENT VARIABLE LABEL: X-\n'
    content += f' CHANNEL NO                    ;-INDEPENDENT VARIABLE LABEL: Y-\n'
    content += f' MIXED BEAM DATA  MIXED        ;-DEPENDENT VARIABLE LABEL-\n'
    content += f' 0                             ;-PROC CODE- 0:RAW 1:AVG 2:SM 3:AVG+SM\n'
    content += f'       {nt}                    ;-# OF X PTS-\n'
    content += f'          {nchan}                    ;-# OF Y PTS- X,Y,F(X,Y) DATA FOLLOW:\n'

    # write time data
    for irow in np.arange(int(np.ceil(len(t) / 6))):  # 6 columns
        content += '  ' + '  '.join([f'{val:.5e}' for val in t[irow * 6:irow * 6 + 6]]) + '\n'
    # write channel array
    for irow in np.arange(int(np.ceil(nchan / 6))):  # single row for us with only 1 beam, but leave it general
        content += '  ' + '  '.join([f'{val + 1:.5e}' for val in np.arange(nchan)[irow * 6:irow * 6 + 6]]) + '\n'
    # write data arrays
    data = np.append(pb, np.append(eb, np.append(ffull, fhalf)))
    for irow in np.arange(int(np.ceil(len(data) / 6))):
        content += '  ' + '  '.join([f'{val:.5e}' for val in data[irow * 6:irow * 6 + 6]]) + '\n'

    content += ';-------END OF DATA-------------\n'
    content += f'Written by python routine on {datetime.datetime.today().strftime("%d %b %Y")}\n'
    content += f'Author: Bill Capecchi wcapecch@pppl.gov'

    writeto = io.open(tofile, 'w', newline='\n')
    writeto.write(content)
    writeto.close()
    print(f'wrote {tofile}')


if __name__ == '__main__':
    shot = 108890
    # read_nbi_ufile(f'{os.getcwd()}/data/TVN{shot}.NBI')

    # write_nbi_ufile(108890, alphanum='V01')
    read_nbi_ufile(shot, pre='PVN')

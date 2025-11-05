import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from MDSplus import Tree

matplotlib.use('TkAgg')  # allows plotting in debug mode
clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']


def get_anpa_rawdata(shot, tree_name="ltxb", npt=None, cal=False):
    """
    Load ANPA raw detector signals from MDSplus.

    todo: The original code digitizes current and voltage in det1, det2, target, lens, capacitor. We should do that

    Parameters
    ----------
    shot : int
        Shot number to load
    tree_name : str, optional
        Name of the MDSplus tree (default: "ltxb")
    npt : int, optional
        Number of time points to load from each channel (default was: 128000)

    Returns
    -------
    anpa : ndarray, shape (10, npt)
        Raw detector waveforms.
        Channels 0–9 correspond to H detectors (H1–H10)
        Channels 10–19 correspond to D detectors (D1–D10)  # OMITTED FOR LTX
    t_anpa : ndarray, shape (npt)

    Notes
    -----
    This version does NOT apply bias-dependent gain correction.
    A placeholder function is included below to add it later.
    """

    noise_signal = 'det'
    anpa_noise_coeff_file = '/home/eilerman/anpa_proc/noise/ec.dat'
    rel_cal = [0.315, 0.389, 0.56, 1.000, 0.955, 1.177, 1.212, 1.065, 1.106, 0.502]  # H channels only

    tree = Tree(tree_name, shot)

    # Output: 10 detectors × npt samples (H only)
    if npt is None:
        npt = len(tree.getNode('\\anpa_h_1').data())
    print(f'data length: {npt} pts')
    anpa = np.zeros((10, npt), dtype=float)

    # --- Load Hydrogen channels (\anpa_h_1 ... \anpa_h_10)
    for i in range(10):
        node = f"\\anpa_h_{i + 1}"
        sig = tree.getNode(node).data()
        anpa[i, :] = sig[:npt]
    t_anpa = tree.get('dim_of(\\anpa_h_1)').data()

    if cal:
        for i in np.arange(10):
            anpa[i, :] /= rel_cal[i]

    return anpa, t_anpa


# Optional: place-holder to add bias compensation later
def apply_bias_compensation(anpa, bias_voltages, calibration_curve):
    """
    Placeholder for bias compensation logic.

    Parameters
    ----------
    anpa : ndarray
        Raw detector data (output of get_anpa_rawdata)
    bias_voltages : array-like, length 20
        Bias voltage applied to each detector during the shot
    calibration_curve : callable or lookup table
        Function G(V) mapping bias voltage to detector gain

    Returns
    -------
    corrected : ndarray
        Bias-corrected detector signals

    Notes
    -----
    You will replace this once you locate your gain calibration function.
    """
    # Example: match all detectors to a reference bias V_ref
    V_ref = np.mean(bias_voltages)

    gain_factors = calibration_curve(V_ref) / calibration_curve(np.array(bias_voltages))[:, None]
    return anpa * gain_factors


def energy_calibration():
    """
    :return:
      energy: peaked energy for each channel (keV)
      low_energy: lower limit of energy for each channel (keV)
      high_energy: upper limit of energy for each channel (keV)
    """
    energy = np.array([11.48, 14.17, 17.11, 20.69, 23.93, 27.23, 30.96, 34.82, 38.70, 43.56])  # H
    three_sigma = np.array([1.85, 2.62, 3.61, 4.73, 5.05, 5.29, 5.61, 5.89, 5.98, 6.92])  # H only
    low_energy = energy - three_sigma
    high_energy = energy + three_sigma
    return energy, low_energy, high_energy


def apply_filter(anpa):
    # optional place to apply Savortsky filter or whatever
    return anpa


def plot_anpa_rawdata(shot, logscale=False, cmap='inferno'):
    anpa, t_anpa = get_anpa_rawdata(shot)  # anpa shape = [10 channels, time]
    anpa = apply_filter(anpa)
    energy, low_energy, high_energy = energy_calibration()

    # todo: need to digitize target, lens, capacitor, and detector voltages!
    v_target = 0.  # foil voltage (for shifting energy sensitivity range)
    energy = energy - v_target
    # low_energy = low_energy - v_target
    # high_energy = high_energy - v_target

    Z = np.where(anpa>0, anpa, np.nan)
    if logscale:
        Z = np.log10(Z)

    plt.figure(figsize=(9,5))
    t_edges = np.concatenate([t_anpa, [t_anpa[-1] + (t_anpa[-1]-t_anpa[-2])]])
    e_edges = np.concatenate([energy, [energy[-1] + (energy[-1]-energy[-2])]])

    plt.pcolormesh(t_edges, e_edges, Z, shading='auto', cmap=cmap)
    plt.xlabel("Time (s)")
    plt.ylabel("Energy (eV)")
    plt.title("ANPA Energy Spectrum")
    cb = plt.colorbar()
    cb.set_label("log10(counts)" if logscale else "counts")
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    plot_anpa_rawdata(112562)

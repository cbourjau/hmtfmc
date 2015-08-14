"""
This file provides functions to extract data from the "Sums" histogram. The functions should, return (lists of)
plottables.
"""

import sys, os, string, random

from rootpy import asrootpy
from rootpy.plotting import HistStack, Canvas, Legend, Pad, Hist1D, Hist2D, Graph

from post_utils import gen_random_name

import ROOT


def get_dNdeta_binned_in_mult(h2d, event_counter, nch_max, with_mb=True):
    """
    Create dN/deta stack for various multiplicity bins from given 2D histogram. If `with_mb` is `True`,
    also add dNdeta for minimum bias to the stack.
    Parameters
    ----------
    h2d : Hist2D
          feta_Nch from the Sums list
    event_counter : Hist1D
         Event counter histogram used for scaling
    Returns
    -------
    list
          List of plottables
    """
    stack = list()
    nbins = h2d.yaxis.GetNbins()
    nch_step = 10
    last_mult_bin = False
    for mult_bin in range(1, nbins, nch_step):
        if mult_bin + nch_step > nch_max:
            last_mult_bin = True
            mult_bin_upper = nbins
        else:
            mult_bin_upper = mult_bin + nch_step - 1
        h2d.yaxis.set_range(mult_bin, mult_bin_upper)
        h = asrootpy(h2d.projection_x())
        h.name = str(mult_bin)
        h.title = (str(int(h2d.yaxis.get_bin_low_edge(mult_bin)))
                   + ' #leq N_{ch}^{est} #leq ' +
                   str(int(h2d.yaxis.get_bin_up_edge(mult_bin_upper))))
        h.yaxis.title = '1/N dN_{ch}^{est}/d#eta'
        # scale by the number of events in this mult_bin
        try:
            h.Scale(1.0 / float(event_counter.Integral(mult_bin, mult_bin_upper)))
        except ZeroDivisionError:
            raise ZeroDivisionError("Consider lowering the cutoff factor to avoid this")
        stack.append(h)
        if last_mult_bin:
            break
    if with_mb:
        h2d.yaxis.set_range(0, -1)  # full Nch range
        h = asrootpy(h2d.projection_x())
        h.title = 'MB'
        h.name = 'mb_dNdeta'
        h.yaxis.title = '1/N dN_{ch}^{est}/d#eta'
        h.Scale(1.0 / float(event_counter.Integral(0, -1)))
        stack.append(h)
    return stack


def get_identified_vs_mult(h3d, pdg):
    """
    Return 1D counter histogram of identified particles vs N_ch^est
    Parameters
    ----------
    h3d: Hist3D
         x: est_mult; y: pT; z: pids1
    pdg: str
         pdg code as string
    Return
    ------
    Hist1D:
         x: Nch_est y: counts
    """
    pid_bin = h3d.zaxis.find_bin(pdg)
    if pid_bin == 0:
        raise ValueError("given pdg ({}) does not exist in histogram".format(pdg))

    h3d.zaxis.SetRange(pid_bin, pid_bin)
    h = asrootpy(h3d.Project3D("yx").ProjectionX())
    h.SetName(gen_random_name())
    return h


def get_NchEst1_vs_NchEst2(results_post, est1, est2):
    """
    Returns the correlation histogram between estimator 1 and estimator two
    Parameters
    ----------
    est1 : str
           Name of estimator 1
    est2 : str
           Name of estimator 2
    Returns
    -------
    Hist2D :
            Hist2D with Nch est1 on the x- and Nch est2 on the y-axis_title
    """
    corr_hist_dir = results_post.Get("correlations")
    hist_name = "corr_hist_{}_vs_{}".format(est1, est2)
    corr_hist = asrootpy(corr_hist_dir.Get(hist_name))
    corr_hist.name = gen_random_name()
    if not isinstance(corr_hist, ROOT.TH2):
        raise ValueError("Correlation histogram with name {} not found".format(hist_name))
    return corr_hist


def get_PNch_vs_estmult(results_post, est):
    """
    Parameters
    ----------
    results_post : TDirectory
                   Results_post directory
    est : str
          Estimator name
    Returns
    -------
    Hist1D :
            Counter Histogram for Number of events with Nch in the estimator region
    """
    # nasty hardcoded:
    ref_est = "EtaLt05"
    corr_hist = get_NchEst1_vs_NchEst2(results_post, ref_est, est)
    return asrootpy(corr_hist.ProjectionX(gen_random_name()))

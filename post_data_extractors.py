"""
This file provides functions to extract data from the "Sums" histogram. The functions should, return (lists of)
plottables.
"""

import sys, os, string, random

from rootpy import asrootpy
from rootpy.plotting import HistStack, Canvas, Legend, Pad, Hist1D, Hist2D, Graph

from post_utils import gen_random_name

import ROOT


def get_Nch_edges_for_percentile_edges(percentile_edges, event_counter):
    """
    Returns the bin edges so that each bin reprecents the given percentile edges
    Parameters
    ----------
    percentile_edges : list
                       Edges in percentiles; high percentile comes first
    event_counter : Hist1D
                    Event counter histogram with Nch (in estimator region, binwidth=1) on the x-axis
    Returns
    -------
    list :
           Edges in Nch of the estimator, first edge is lowest in Nch, last edge the highest in Nch
    """
    nch_bins = event_counter.GetXaxis().GetNbins()
    n_total = event_counter.Integral(1, nch_bins)
    accu_counts = [0] + [event_counter.Integral(1, binidx) for binidx in range(1, nch_bins)]
    # percentage of events having nch=this_idx or more Nch; hence starts at 1.0 at idx=0 and decreases
    # but the floats might not equal 0.0!!!
    decumulative_percents = [1 - (accu_count / float(n_total)) for accu_count in accu_counts]
    nch_edges = []  # int(round(n_total - n_total * perc)) for perc in percentile_edges]
    for perc in percentile_edges:
        # find the idx where the given percentile edge would fit into the decumulative list
        for idx, decum_perc in enumerate(decumulative_percents):
            if perc > decum_perc:
                # this should pick up 1.0
                nch_edges.append(idx - 1)
                break
        else:
            # the loop did not break, then the last bin is _probably_ 0.0?
            if perc == 0.0:
                nch_edges.append(nch_bins - 1)
            else:
                # something is fishy!
                raise IndexError("Something is fishy here!!!")
    return nch_edges


def get_meanpt_vs_estmult(resutlts_est_dir, pids):
    """
    Create a 1Dprofile for the given pids and the given estimator name
    """
    # find the mult vs pt histograms for the given pids
    mult_vs_pts = []
    for pid in pids:
        mult_vs_pts.append(asrootpy(getattr(resutlts_est_dir.mult_pt, pid)))
    profx = sum(mult_vs_pts).ProfileX()
    profx.name = gen_random_name()
    return profx


def get_dNdeta_in_mult_interval(h2d, event_counter, mult_interval):
    """
    Create dN/deta stack for various multiplicity bins from given 2D histogram. If `with_mb` is `True`,
    also add dNdeta for minimum bias to the stack.
    Parameters
    ----------
    h2d : Hist2D
          feta_Nch from the Sums list
    event_counter : Hist1D
                    Event counter histogram used for scaling
    mult_interval : list
               Nch bin edges
    Returns
    -------
    Hist1D
    """
    h2d.yaxis.set_range(mult_interval[0] + 1, mult_interval[1] + 1)
    h = asrootpy(h2d.projection_x(gen_random_name()))
    h.title = "{} - {} %".format(100 * mult_interval[0], 100 * mult_interval[1])
    # scale by the number of events in this mult_interval and bin width
    try:
        h.Scale(1.0 / float(event_counter.Integral(mult_interval[0] + 1, mult_interval[1] + 1)), "width")
    except ZeroDivisionError:
        raise ZeroDivisionError("Consider increasing the nch interval to avoid this")
    return h


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


def get_NchEst1_vs_NchEst2(sums, est1, est2):
    """
    Returns the correlation histogram between estimator 1 and estimator two
    Parameters
    ----------
    sums : TList
           Sums directory
    est1 : str
           Name of estimator 1
    est2 : str
           Name of estimator 2
    Returns
    -------
    Hist2D :
            Hist2D with Nch est1 on the x- and Nch est2 on the y-axis_title
    """
    if not isinstance(sums, ROOT.TList):
        raise TypeError("{} is not of type ROOT.TList".format(sums))
    hist_name = "fcorr_thisNch_vs_refNch"
    corr_hist = asrootpy(sums.FindObject(est2).FindObject(hist_name))
    if not isinstance(corr_hist, ROOT.TH2):
        raise ValueError("Correlation histogram with name {} not found".format(hist_name))
    return corr_hist


def get_PNch_vs_estmult(sums, est):
    """
    Parameters
    ----------
    sums : TList
           Sums directory
    est : str
          Estimator name
    Returns
    -------
    Hist1D :
            Counter Histogram for Number of events with Nch in the estimator region
    """
    if not isinstance(sums, ROOT.TList):
        raise TypeError("{} is not of type ROOT.TList".format(sums))
    # nasty hardcoded:
    ref_est = "EtaLt05"
    corr_hist = get_NchEst1_vs_NchEst2(sums, ref_est, est)
    return asrootpy(corr_hist.ProjectionX(gen_random_name()))


def get_nMPI_vs_Nch(sums_est_dir):
    """
    Parameters
    ----------
    sums_est_dir : TList
        Estimator directory from Sums

    Returns
    -------
    Profile :
        Profile with N_ch on the xaxis and <nMPI> on the yaxis
    """
    nch_nmpi = sums_est_dir.FindObject("fNch_vs_nMPI")
    profx = nch_nmpi.ProfileX()
    return profx


def get_pT_distribution(results_est_dir, pids, nch_low, nch_up, normalized):
    """
    Parameters
    ----------
    results_est_dir : TDirectory
               Directory of a given estimator
    pids : list
           List of strings denoting requested pids
    nch_low, nch_up : int
           Lower and upper limit of Nch for which the p_T distribution should be made
    normalized : Boolean
           Should the distribution be normalized to yield P(p_T)?
    Returns
    -------
    Hist1D :
            Histogram P(p_T)
    """
    mult_pt_hists = []
    for pid in pids:
        mult_pt_hists.append(getattr(results_est_dir.mult_pt, pid))
    summed_mult_pt = sum(mult_pt_hists)
    summed_mult_pt.xaxis.SetRangeUser(nch_low, nch_up)
    projy = asrootpy(summed_mult_pt.ProjectionY())
    projy.name = gen_random_name()
    event_counter = asrootpy(results_est_dir.event_counter)
    # Scale by the number of events in the interval; N_ch == bin - 1
    projy.Scale(1.0 / event_counter.Integral(nch_low, nch_up))
    if normalized:
        projy.Scale(1.0 / projy.Integral())
    return projy


def get_mean_nMPI(sums_est_dir, nch_low, nch_up):
    """
    Get the mean nMPI of events in a given N_ch interval
    Parameters
    ----------
    sums_est_dir : TList
                   List for a given estimator
    nch_low, nch_up : int
           Lower and upper limit of Nch for which <nMPI> should be calculated
    Returns
    -------
    Float :
           <nMPI>
    """
    nch_vs_nmpi = asrootpy(sums_est_dir.FindObject("fNch_vs_nMPI"))
    nch_vs_nmpi.xaxis.SetRangeUser(nch_low, nch_up)
    return nch_vs_nmpi.GetMean(2)

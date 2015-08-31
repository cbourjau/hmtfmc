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
                       Edges in percentiles; has to start at 1 and end at 0
    event_counter : Hist1D
                    Event counter histogram with Nch (in estimator region) on the x-axis
    Returns
    -------
    list :
           Edges in Nch of the estimator, first edge is 0 last edge the highest available bin
    """
    nch_bins = event_counter.GetXaxis().GetNbins()
    n_total = event_counter.Integral(1, nch_bins)
    accu_counts = [0] + [event_counter.Integral(1, binidx) for binidx in range(1, nch_bins)]
    accu_percents = [1 - (accu_count / float(n_total)) for accu_count in accu_counts]

    nch_edges = []  # int(round(n_total - n_total * perc)) for perc in percentile_edges]
    iter_wanted_perc_edges = iter(percentile_edges)
    perc_edge = next(iter_wanted_perc_edges)
    for idx, (accu_count, accu_percent) in enumerate(zip(accu_counts, accu_percents)):
        # the idx is the nch bin we are looking for
        try:
            if accu_percent <= perc_edge:  # this should pick up the first bin with 1.0 = 1.0
                nch_edges.append(idx)
                perc_edge = next(iter_wanted_perc_edges)
                # check if we are skipping a bin because the percentile intervals are too small:
                if accu_percent <= perc_edge:
                    raise ValueError("The given percentile edges narrower than one Nch bin!")
        except StopIteration:
            # We found the last bin
            break
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


def get_dNdeta_binned_in_mult(h2d, event_counter, percent_bins=None, nch_max=None, with_mb=True):
    """
    Create dN/deta stack for various multiplicity bins from given 2D histogram. If `with_mb` is `True`,
    also add dNdeta for minimum bias to the stack.
    Parameters
    ----------
    h2d : Hist2D
          feta_Nch from the Sums list
    event_counter : Hist1D
                    Event counter histogram used for scaling
    percent_bins : list
                   Percentile bin edges; either this or nch_max has to be given
    nch_max : int
              The cutoff value above which, everything is packed in one bin, if percent_bins is given,
              this value is ignored
    with_mb : bool
              If true, the last plottable in the return list is the minimum bias distribution

    Returns
    -------
    list
          List of plottables
    """
    stack = list()
    nbins = h2d.yaxis.GetNbins()
    nch_step = 10
    last_mult_bin = False
    if percent_bins:
        nch_edges = get_Nch_edges_for_percentile_edges(percent_bins, event_counter)
    else:
        nch_edges = range(0, nbins, nch_step)

    for i, (nch_low, nch_up) in enumerate(zip(nch_edges[:-1], nch_edges[1:])):
        nch_bin_low = nch_low + 1  # root's crack pot binning...
        nch_bin_up = nch_up + 1    # root's crack pot binning...
        if percent_bins is None:
            if (nch_bin_up > nch_max):
                last_mult_bin = True
                nch_bin_up = nbins
        h2d.yaxis.set_range(nch_bin_low, nch_bin_up)
        h = asrootpy(h2d.projection_x())
        h.name = str(nch_bin_low)
        if percent_bins is None:
            h.title = (str(int(nch_low))
                       + '#leqN_{ch}^{est}#leq' +
                       str(int(nch_up)))
        else:
            h.title = "{} - {} %".format(100 * percent_bins[i], 100 * percent_bins[i + 1])
        h.yaxis.title = '1/N dN_{ch}^{est}/d#eta'
        # scale by the number of events in this mult_bin
        try:
            h.Scale(1.0 / float(event_counter.Integral(nch_bin_low, nch_bin_up)))
        except ZeroDivisionError:
            raise ZeroDivisionError("Consider lowering the cutoff factor / percentil binning to avoid this")
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
        import ipdb; ipdb.set_trace()
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

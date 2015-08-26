"""
This is the entry point to run the post analysis. This file contains the steering part.
Ie. it does define which plots to make, but not how. The logic which extracts the data from the
primary generated "Sums" is in the post_data_extractor.py. These functions provide (lists of) plottables.
These plottables are then plotted with the plotting functions in plotting_util.py
"""

import sys
import os
import string
import random
import ipdb

if len(sys.argv) != 2:
    print "Usage: python ./post.py path_to_root_file.root"
    quit()

from rootpy.io import root_open
from rootpy import asrootpy, ROOT, log
from rootpy.plotting import Hist1D, Hist2D, Canvas

from post_data_extractors import get_dNdeta_binned_in_mult, get_identified_vs_mult,\
    get_Nch_edges_for_percentile_edges,\
    get_NchEst1_vs_NchEst2, get_PNch_vs_estmult
from post_utils import create_stack_pid_ratio_over_pt,\
    remap_x_values,\
    remove_zero_value_points, remove_non_mutual_points,\
    remove_points_with_equal_x, remove_points_with_x_err_gt_1NchRef

from roofi import Figure


def gen_random_name():
    """Generate a random name for temp hists"""
    return ''.join(random.choice(string.ascii_letters + string.digits) for _ in range(25))


def make_estimator_title(name):
    if name == 'EtaLt05':
        return '|#eta|#leq0.5'
    elif name == 'EtaLt08':
        return '|#eta|#leq0.8'
    elif name == 'EtaLt15':
        return '|#eta|#leq1.5'
    elif name == 'Eta08_15':
        return '0.8#leq|#eta|#leq1.5'
    else:
        return name


kPROTON = str(2212)
kANTIPROTON = str(-2212)
kLAMBDA = str(3122)
kANTILAMBDA = str(-3122)
kK0S = str(310)
kKPLUS = str(321)
kKMINUS = str(-321)
kPIPLUS = str(211)
kPIMINUS = str(-211)
kPI0 = str(111)
kXI = str(3312)
kANTIXI = str(-3312)
kOMEGAMINUS = str(3334)
kOMEGAPLUS = str(-3334)


# use the last mult bin starts at a multiplicity  x times larger than the mean in this estimator
mean_mult_cutoff_factor = 3


def get_est_dirs(sums):
    return (somedir for somedir in sums if somedir.GetName() in considered_ests)


def _plot_particle_ratios_vs_estmult(f, sums, results_post, pids1, pids2, scale=None, ytitle=''):
    ratio_vs_estmult_dir = 'MultEstimators/results_post/pid_ratios_vs_estmult'
    fig = Figure()
    if not ytitle:
        fig.ytitle = ", ".join(pids1) + " / " + ", ".join(pids2)
    else:
        fig.ytitle = ytitle

    for est_dir in sums:
        if est_dir.GetName() not in ['EtaLt05', 'EtaLt08', 'EtaLt15', 'Eta08_15', 'V0M']:
            continue
        h3d = asrootpy(est_dir.FindObject("fNch_pT_pid"))
        pids1hists = [get_identified_vs_mult(h3d, pdg) for pdg in pids1]
        pids2hists = [get_identified_vs_mult(h3d, pdg) for pdg in pids2]

        pids1_px = sum(pids1hists)
        pids2_px = sum(pids2hists)
        ratio1d = pids1_px / pids2_px

        fig.xtitle = "N_{ch}|_{" + make_estimator_title(est_dir.GetName()) + "}"

        if scale:
            ratio1d.Scale(scale)
        fig.add_plottable(ratio1d, legend_title=make_estimator_title(est_dir.GetName()))
    name = "_".join(pids1) + "_div_" + "_".join(pids2)
    fig.save_to_root_file(f, name, ratio_vs_estmult_dir)


def _plot_particle_ratios_vs_refmult(f, sums, results_post, pids1, pids2, scale=None, ytitle=''):
    """
    plot and write to file the ratio of the two pid-lists (pids1/pids2). Plot is vs refmult.
    This function depends on the correlation histograms to be present in f
    """
    ratio_vs_refmult_dir = 'MultEstimators/results_post/pid_ratios_vs_refmult'
    fig = Figure()

    refest = "EtaLt05"
    if not ytitle:
        fig.ytitle = ", ".join(pids1) + " / " + ", ".join(pids2)
    else:
        fig.ytitle = ytitle

    for est_dir in get_est_dirs(sums):
        h3d = asrootpy(est_dir.FindObject("fNch_pT_pid"))
        corr_hist = asrootpy(est_dir.FindObject("fcorr_thisNch_vs_refNch"))

        pids1_vs_estmult = sum([get_identified_vs_mult(h3d, pdg) for pdg in pids1])
        pids2_vs_estmult = sum([get_identified_vs_mult(h3d, pdg) for pdg in pids2])

        # remap histograms using the correlation between the current estimator and the reference one
        pids1_vs_refmult = remap_x_values(pids1_vs_estmult, corr_hist)
        pids2_vs_refmult = remap_x_values(pids2_vs_estmult, corr_hist)

        # sanitize
        for g in [pids1_vs_refmult, pids2_vs_refmult]:
            remove_zero_value_points(g)
            remove_points_with_x_err_gt_1NchRef(g)
            remove_points_with_equal_x(g)
        remove_non_mutual_points(pids1_vs_refmult, pids2_vs_refmult)

        try:
            ratio = pids1_vs_refmult / pids2_vs_refmult
        except ZeroDivisionError:
            print "ZeroDivisionError in {}".format(est_dir.GetName())
            continue
        if scale:
            ratio.Scale(scale)
        fig.xtitle = "N_{ch}|_{" + refest + "}"
        fig.add_plottable(ratio, make_estimator_title(est_dir.GetName()))

    name = "_".join(pids1) + "_div_" + "_".join(pids2)
    fig.save_to_root_file(f, name, ratio_vs_refmult_dir)


def _make_event_counters(f, sums, results_post):
    log.info("Creating event counters")
    for est_dir in get_est_dirs(sums):
        results_est_dir = results_post.FindObject(est_dir.GetName())
        corr = asrootpy(est_dir.FindObject("fcorr_thisNch_vs_refNch"))
        counter = asrootpy(corr.ProjectionX())
        counter.name = "event_counter"
        f.cd("MultEstimators/results_post/" + est_dir.GetName())
        results_est_dir.WriteTObject(counter)
        #counter.Write()


def _make_dNdeta(f, sums, results_post):
    # Loop over all estimators in the Sums list:
    log.info("Creating dN/deta bin in multiplicity")
    for est_dir in get_est_dirs(sums):
        results_est_dir = results_post.Get(est_dir.GetName())
        h3d = asrootpy(est_dir.FindObject('fNch_pT_pid'))  # use FindObject on List objects (they have no get)
        h2d = asrootpy(est_dir.FindObject('feta_Nch'))
        event_counter = asrootpy(results_est_dir.Get("event_counter"))

        mean_nch = est_dir.FindObject("feta_Nch").GetMean(2)  # mean of yaxis
        # bin in standard step size up to max_nch; from there ibs all in one bin:
        max_nch = mean_nch * mean_mult_cutoff_factor

        fig = Figure()
        perc_edges = [1, .6, .4, .2, .1, .05, .025, 0]
        hists = get_dNdeta_binned_in_mult(h2d, event_counter, percent_bins=perc_edges,  # nch_max=max_nch,
                                          with_mb=True)
        [fig.add_plottable(p, legend_title=p.title) for p in hists]
        fig.xtitle = '#eta'
        fig.ytitle = 'dN_{ch}/d#eta'
        fig.legend.position = 'seperate'
        path = results_est_dir.GetPath().split(":")[1]  # file.root:/internal/root/path
        fig.save_to_root_file(f, "dNdeta_summary", path=path)


def _make_hists_vs_pt(f, sums, results_post):
    ###########################################################
    # Category 2 on TWiki
    # create particle ratio vs pT plots

    log.info("Computing histograms vs pt")

    # Loop over all estimators in the Sums list:
    for est_dir in get_est_dirs(sums):
        dirname = 'MultEstimators/results_post/{}/pid_ratios/'.format(est_dir.GetName())
        try:
            f.mkdir(dirname, recurse=True)
        except:
            pass
        f.cd(dirname)
        h3d_orig = asrootpy(est_dir.FindObject('fNch_pT_pid'))
        event_counter = asrootpy(getattr(results_post, est_dir.GetName()).event_counter)
        perc_edges = [1, .6, .4, .2, .1, .05, .025, 0]

        nch_edges = get_Nch_edges_for_percentile_edges(perc_edges, event_counter)
        # mean_nch = est_dir.FindObject("feta_Nch").GetMean(2)  # mean of yaxis
        # # bin in standard step size up to max_nch; from there ibs all in one bin:
        # nch_cutoff = int(mean_nch * mean_mult_cutoff_factor)
        # step_size = 10
        # nch_edges = list(range(0, nch_cutoff, step_size)) + [h3d_orig.GetXaxis().GetNbins()]
        esti_title = "({0})".format(h3d_orig.title[31:])

        mult_pt_dir = results_post.FindObject(est_dir.GetName()).Get("mult_pt")
        fig = Figure()
        fig.xtitle = 'p_{T} (GeV)'
        fig.plot.ymin = 0
        fig.plot.palette = 'root'
        fig.plot.palette_ncolors = len(nch_edges) - 1
        fig.legend.position = 'bl'

        fig.delete_plottables()
        name = "proton_over_pich__vs__pt"
        hs = create_stack_pid_ratio_over_pt(mult_pt_dir, [kANTIPROTON, kPROTON], [kPIMINUS, kPIPLUS], nch_edges)
        fig.ytitle = "(p+#bar{p})/#pi^{+-}"
        fig.plot.ymax = .3
        fig.legend.title = make_estimator_title(est_dir.GetName())
        [h.SetTitle(str(int(100 * up)) + "-" + str(int(100 * low)) + "%")
         for h, up, low in zip(hs, perc_edges[:-1], perc_edges[1:])]
        [fig.add_plottable(p, p.title) for p in hs]
        fig.save_to_root_file(f, name, dirname)

        fig.delete_plottables()
        name = "Xi_over_pich__vs__pt"
        fig.plot.ymax = .06
        fig.legend.position = 'tl'
        hs = create_stack_pid_ratio_over_pt(mult_pt_dir, [kANTIXI, kXI], [kPIMINUS, kPIPLUS], nch_edges)
        fig.ytitle = "#Xi/#pi^{+-}"
        fig.legend.title = make_estimator_title(est_dir.GetName())
        [h.SetTitle(str(int(100 * up)) + "-" + str(int(100 * low)) + "%")
         for h, up, low in zip(hs, perc_edges[:-1], perc_edges[1:])]
        [fig.add_plottable(p, p.title) for p in hs]
        fig.save_to_root_file(f, name, dirname)

        fig.delete_plottables()
        name = "OmegaCh_over_pich__vs__pt"
        fig.plot.ymax = .005
        fig.legend.position = 'tl'
        hs = create_stack_pid_ratio_over_pt(mult_pt_dir, [kOMEGAMINUS, kOMEGAPLUS], [kPIMINUS, kPIPLUS], nch_edges)
        fig.ytitle = "#Omega_{ch}/#pi^{+-} "
        fig.legend.title = make_estimator_title(est_dir.GetName())
        [h.SetTitle(str(int(100 * up)) + "-" + str(int(100 * low)) + "%")
         for h, up, low in zip(hs, perc_edges[:-1], perc_edges[1:])]
        [fig.add_plottable(p, p.title) for p in hs]
        fig.save_to_root_file(f, name, dirname)

        # Ratios to pi0
        fig.delete_plottables()
        name = "pich_over_pi0__vs__pt"
        fig.plot.ymax = 2.5
        fig.legend.position = 'bl'
        hs = create_stack_pid_ratio_over_pt(mult_pt_dir, [kPIMINUS, kPIPLUS], [kPI0], nch_edges)
        fig.ytitle = "#pi^{+-}/#pi^{0}"
        fig.legend.title = make_estimator_title(est_dir.GetName())
        [h.SetTitle(str(int(100 * up)) + "-" + str(int(100 * low)) + "%")
         for h, up, low in zip(hs, perc_edges[:-1], perc_edges[1:])]
        [fig.add_plottable(p, p.title) for p in hs]
        fig.save_to_root_file(f, name, dirname)

        fig.delete_plottables()
        name = "proton_over_pi0__vs__pt"
        fig.plot.ymax = 1
        fig.legend.position = 'tr'
        hs = create_stack_pid_ratio_over_pt(mult_pt_dir, [kANTIPROTON, kPROTON], [kPI0], nch_edges)
        fig.ytitle = "p/#pi^{0}"
        fig.legend.title = make_estimator_title(est_dir.GetName())
        [h.SetTitle(str(int(100 * up)) + "-" + str(int(100 * low)) + "%")
         for h, up, low in zip(hs, perc_edges[:-1], perc_edges[1:])]
        [fig.add_plottable(p, p.title) for p in hs]
        fig.save_to_root_file(f, name, dirname)

        fig.delete_plottables()
        name = "K0S_over_pi0__vs__pt"
        fig.plot.ymax = 1.4
        fig.legend.position = 'tl'
        hs = create_stack_pid_ratio_over_pt(mult_pt_dir, [kK0S], [kPI0], nch_edges)
        fig.ytitle = "K^{0}_{S}/#pi^{0}"
        fig.legend.title = make_estimator_title(est_dir.GetName())
        [h.SetTitle(str(int(100 * up)) + "-" + str(int(100 * low)) + "%")
         for h, up, low in zip(hs, perc_edges[:-1], perc_edges[1:])]
        [fig.add_plottable(p, p.title) for p in hs]
        fig.save_to_root_file(f, name, dirname)

        fig.delete_plottables()
        name = "Lambda_over_pi0__vs__pt"
        fig.plot.ymax = .9
        fig.legend.position = 'tl'
        hs = create_stack_pid_ratio_over_pt(mult_pt_dir, [kANTILAMBDA, kLAMBDA], [kPI0], nch_edges)
        fig.ytitle = "#Lambda/#pi^{0}"
        fig.legend.title = make_estimator_title(est_dir.GetName())
        [h.SetTitle(str(int(100 * up)) + "-" + str(int(100 * low)) + "%")
         for h, up, low in zip(hs, perc_edges[:-1], perc_edges[1:])]
        [fig.add_plottable(p, p.title) for p in hs]
        fig.save_to_root_file(f, name, dirname)

        fig.delete_plottables()
        name = "Xi_over_pi0__vs__pt"
        fig.plot.ymax = .08
        fig.legend.position = 'tl'
        hs = create_stack_pid_ratio_over_pt(mult_pt_dir, [kANTIXI, kXI], [kPI0], nch_edges)
        fig.ytitle = "#Xi/#pi^{0}"
        fig.legend.title = make_estimator_title(est_dir.GetName())
        [h.SetTitle(str(int(100 * up)) + "-" + str(int(100 * low)) + "%")
         for h, up, low in zip(hs, perc_edges[:-1], perc_edges[1:])]
        [fig.add_plottable(p, p.title) for p in hs]
        fig.save_to_root_file(f, name, dirname)

        fig.delete_plottables()
        name = "OmegaCh_over_pi0__vs__pt"
        fig.plot.ymax = .005
        fig.legend.position = 'tl'
        hs = create_stack_pid_ratio_over_pt(mult_pt_dir, [kOMEGAMINUS, kOMEGAPLUS], [kPI0], nch_edges)
        fig.ytitle = "#Omega_{ch}/#pi^{0}"
        fig.legend.title = make_estimator_title(est_dir.GetName())
        [h.SetTitle(str(int(100 * up)) + "-" + str(int(100 * low)) + "%")
         for h, up, low in zip(hs, perc_edges[:-1], perc_edges[1:])]
        [fig.add_plottable(p, p.title) for p in hs]
        fig.save_to_root_file(f, name, dirname)

        # Ratios to K0S
        fig.delete_plottables()
        name = "proton_over_K0S__vs__pt"
        fig.plot.ymax = 2.6
        fig.legend.position = 'tr'
        hs = create_stack_pid_ratio_over_pt(mult_pt_dir, [kANTIPROTON, kPROTON], [kK0S], nch_edges)
        fig.ytitle = "p/K^{0}_{S}"
        fig.legend.title = make_estimator_title(est_dir.GetName())
        [h.SetTitle(str(int(100 * up)) + "-" + str(int(100 * low)) + "%")
         for h, up, low in zip(hs, perc_edges[:-1], perc_edges[1:])]
        [fig.add_plottable(p, p.title) for p in hs]
        fig.save_to_root_file(f, name, dirname)

        fig.delete_plottables()
        name = "Lambda_over_K0S__vs__pt"
        fig.plot.ymax = 1
        fig.legend.position = 'bl'
        hs = create_stack_pid_ratio_over_pt(mult_pt_dir, [kANTILAMBDA, kLAMBDA], [kK0S], nch_edges)
        fig.ytitle = "#Lambda/K^{0}_{S}"
        fig.legend.title = make_estimator_title(est_dir.GetName())
        [h.SetTitle(str(int(100 * up)) + "-" + str(int(100 * low)) + "%")
         for h, up, low in zip(hs, perc_edges[:-1], perc_edges[1:])]
        [fig.add_plottable(p, p.title) for p in hs]
        fig.save_to_root_file(f, name, dirname)

        fig.delete_plottables()
        name = "Xi_over_K0S__vs__pt"
        fig.plot.ymax = .2
        fig.legend.position = 'tl'
        hs = create_stack_pid_ratio_over_pt(mult_pt_dir, [kANTIXI, kXI], [kK0S], nch_edges)
        fig.ytitle = "#Xi/K^{0}_{S}"
        fig.legend.title = make_estimator_title(est_dir.GetName())
        [h.SetTitle(str(int(100 * up)) + "-" + str(int(100 * low)) + "%")
         for h, up, low in zip(hs, perc_edges[:-1], perc_edges[1:])]
        [fig.add_plottable(p, p.title) for p in hs]
        fig.save_to_root_file(f, name, dirname)

        fig.delete_plottables()
        name = "OmegaCh_over_K0S__vs__pt"
        fig.plot.ymax = .012
        fig.legend.position = 'tl'
        hs = create_stack_pid_ratio_over_pt(mult_pt_dir, [kOMEGAMINUS, kOMEGAPLUS], [kK0S], nch_edges)
        fig.ytitle = "#Omega_{ch}/K^{0}_{S}"
        fig.legend.title = make_estimator_title(est_dir.GetName())
        [h.SetTitle(str(int(100 * up)) + "-" + str(int(100 * low)) + "%")
         for h, up, low in zip(hs, perc_edges[:-1], perc_edges[1:])]
        [fig.add_plottable(p, p.title) for p in hs]
        fig.save_to_root_file(f, name, dirname)


def _make_PNch_plots(f, sums, results_post):
    log.info("Creating P(Nch) summary plot")
    summary_fig = Figure()
    summary_fig.xtitle = "N_{ch}^{est}"
    summary_fig.ytitle = "P(N_{ch}^{est})"
    summary_fig.legend.position = 'tr'

    for est_name in considered_ests:
        h_tmp = get_PNch_vs_estmult(sums, est_name)
        if h_tmp.Integral() > 0:
            h_tmp.Scale(1.0 / h_tmp.Integral())
            summary_fig.add_plottable(h_tmp, make_estimator_title(est_name))

    summary_fig.plot.logy = True
    path = results_post.GetPath().split(":")[1]  # file.root:/internal/root/path
    summary_fig.save_to_root_file(f, "PNch_summary", path=path)

    log.info("Creating P(Nch_est) and P(Nch_refest) histograms")
    mult_bin_size = 10
    for ref_est in ref_ests:
        for est in considered_ests:
            # Figure properties:
            fig_vs_estmult = Figure()
            fig_vs_refmult = Figure()
            fig_vs_estmult.plot.logy = True
            fig_vs_refmult.plot.logy = True

            fig_vs_estmult.legend.position = 'tr'
            fig_vs_refmult.legend.position = 'tr'

            fig_vs_estmult.xtitle = "N_{{ch}}^{{{}}}".format(est)
            fig_vs_refmult.xtitle = "N_{{ch}}^{{{}}}".format(ref_est)

            fig_vs_estmult.ytitle = "P(N_{{ch}}^{{{}}})".format(est)
            fig_vs_refmult.ytitle = "P(N_{{ch}}^{{{}}})".format(ref_est)

            corr_hist = get_NchEst1_vs_NchEst2(sums, ref_est, est)
            nch_max = corr_hist.xaxis.GetNbins()
            mean_nch_est = corr_hist.GetMean(1)  # mean of x axis
            nch_cutoff = mean_nch_est * mean_mult_cutoff_factor

            is_last_bin = False
            for nch_lower_edge in range(0, nch_max, mult_bin_size):
                if nch_lower_edge + mult_bin_size > nch_cutoff:
                    nch_upper_edge = nch_max
                    is_last_bin = True
                else:
                    nch_upper_edge = nch_lower_edge + mult_bin_size

                # vs est_mult:
                corr_hist.xaxis.SetRange(0, 0)  # reset x axis
                corr_hist.yaxis.SetRange(nch_lower_edge, nch_upper_edge)
                h_vs_est = asrootpy(corr_hist.ProjectionX(gen_random_name()))
                if h_vs_est.Integral() > 0:
                    h_vs_est.Scale(1.0 / h_vs_est.Integral())
                    fig_vs_estmult.add_plottable(h_vs_est,
                                                 "{} < N_{{ch}}^{{{}}} < {}".
                                                 format(nch_lower_edge,
                                                        make_estimator_title(ref_est),
                                                        nch_upper_edge))

                # vs ref_mult:
                corr_hist.yaxis.SetRange(0, 0)  # reset y axis
                corr_hist.xaxis.SetRange(nch_lower_edge, nch_upper_edge)
                h_vs_ref = asrootpy(corr_hist.ProjectionY(gen_random_name()))
                if h_vs_ref.Integral() > 0:
                    h_vs_ref.Scale(1.0 / h_vs_ref.Integral())
                    fig_vs_refmult.add_plottable(h_vs_ref,
                                                 "{} < N_{{ch}}^{{{}}} < {}".
                                                 format(nch_lower_edge,
                                                        make_estimator_title(est),
                                                        nch_upper_edge))

                if is_last_bin:
                    break

            path = results_post.GetPath().split(":")[1] + "/" + est  # file.root:/internal/root/path
            # vs est_mult
            fig_vs_estmult.save_to_root_file(f, "PNch{}_binned_in_Nch{}".format(est, ref_est), path)
            # vs est_mult
            fig_vs_refmult.save_to_root_file(f, "PNch{}_binned_in_Nch{}".format(ref_est, est), path)


def _make_mult_vs_pt_plots(f, sums, results_post):
    log.info("Makeing 2D mult pt plots for each particle kind")
    for est_dir in get_est_dirs(sums):
        path = (results_post.GetPath().split(":")[1]  # file.root:/internal/root/path
                + "/" + est_dir.GetName()
                + "/mult_pt")
        try:
            f.mkdir(path, recurse=True)
        except ValueError:
            pass
        f.cd(path)

        h3d = asrootpy(est_dir.FindObject('fNch_pT_pid'))
        # loop through all particle kinds:
        nPIDs = h3d.zaxis.GetNbins()
        for ibin in range(1, nPIDs + 1):
            h3d.zaxis.SetRange(ibin, ibin)
            mult_pt = asrootpy(h3d.Project3D("yx"))
            mult_pt.name = h3d.zaxis.GetBinLabel(ibin)
            mult_pt.Write()


def _make_dNdeta_mb_ratio_plots(f, sums, results_post):
    # Create ratio plots; depends on the previously created histograms
    log.info("Creating ratios of dN/deta plots for each multiplicity bin")
    for est_dir in (somedir for somedir in results_post if somedir.GetName() in considered_ests):
        res_dir_str = "MultEstimators/results_post/" + est_dir.GetName()
        # get the histograms out of the summary plot, even if it is hackish...
        plot_pad = est_dir.Get("dNdeta_summary").FindObject("plot")
        # some of the hists are just the frames, neglect those
        hists = [asrootpy(h) for h in plot_pad.GetListOfPrimitives()
                 if (isinstance(h, ROOT.TH1) and h.Integral() > 0)]
        # find the first histogram with name mb_dNdeta
        mb_hist_name = 'mb_dNdeta'
        try:
            mb_hist = next(h for h in hists if h.name == mb_hist_name)
        except StopIteration:
            raise StopIteration("no histogram named {} was found in dNdeta summary plot".format(mb_hist_name))
        ratios = [h / mb_hist for h in hists if h.name != mb_hist_name]
        fig = Figure()
        fig.xtitle = '#eta'
        fig.ytitle = 'dN/d#eta|_{mult} / dN/d#eta|_{MB}'
        [fig.add_plottable(p, p.title) for p in ratios]
        name = "dNdeta_ratio_to_mb_canvas"
        fig.save_to_root_file(f, name, res_dir_str)


def _make_correlation_plots(f, sums, results_post):
    # Make correlations between estimators
    log.info("Correlating N_ch of each estimator")
    corr_dir = 'MultEstimators/results_post/correlations'
    try:
        f.mkdir(corr_dir, recurse=True)
    except:
        pass
    # Take ntuple from the first estimator and then add friends to this one
    nt0 = sums[0].FindObject("fEventTuple")
    nt0.SetAlias(sums[0].GetName(), "fEventTuple")

    # build ntuple
    for est_dir in sums[1:]:
        nt0.AddFriend(est_dir.FindObject("fEventTuple"), est_dir.GetName())
    for ref_est in considered_ests:
        for est_dir in sums:
            log.info("Correlating {} with {}".format(ref_est, est_dir.GetName()))
            corr_hist = Hist2D(400, 0, 400,
                               400, 0, 400,
                               name="corr_hist_{}_vs_{}".format(ref_est, est_dir.GetName()))
            # Lables are deliberatly swaped, see Projection below!
            corr_hist.title = ("Correlation N_{{ch}} in {0} and {1};N_{{ch}} {1};N_{{ch}} {0}"
                               .format(ref_est, est_dir.GetName()))

            # this projects onto y:x, to make coding more adventurous
            nt0.Project(corr_hist.name, "{0}.nch:{1}.nch".format(ref_est, est_dir.GetName()),
                        "ev_weight")
            corr_hist.drawstyle = 'colz'
            f.cd(corr_dir)
            corr_hist.write()


def _make_pid_ratio_plots(f, sums, results_post):
    log.info("Creating plots vs refmult")

    # Proton / pi_ch
    _plot_particle_ratios_vs_refmult(f, sums, results_post, ['-2212', '2212'], ['-211', '211'],
                                     ytitle="p/#pi^{+-}")
    # K / pi_ch
    _plot_particle_ratios_vs_refmult(f, sums, results_post, ['310', '321', '-321'], ['-211', '211'],
                                     ytitle="K^{*}/#pi^{+-}")
    # Lambda / pi_ch
    _plot_particle_ratios_vs_refmult(f, sums, results_post, ['3122'], ['-211', '211'],
                                     ytitle="#Lambda / #pi^{+-}")
    # Xi / pi_ch
    _plot_particle_ratios_vs_refmult(f, sums, results_post, ['3312'], ['-211', '211'],
                                     ytitle="#Xi / #pi^{+-}")
    # Omega / pi_ch
    _plot_particle_ratios_vs_refmult(f, sums, results_post, ['3334', '-3334'], ['-211', '211'],
                                     ytitle="#Omega / #pi^{+-}")
    # pi_ch/pi0
    _plot_particle_ratios_vs_refmult(f, sums, results_post, ['-211', '211'], ['111'],
                                     ytitle="#pi^{+-}/#pi^{0}")
    # proton / pi0
    _plot_particle_ratios_vs_refmult(f, sums, results_post, ['-2212', '2212'], ['111'],
                                     ytitle="p/#pi^{0}")
    # K / pi0
    _plot_particle_ratios_vs_refmult(f, sums, results_post, ['310', '321', '-321'], ['111'],
                                     ytitle="K^{*}/#pi^{0}")
    # Lambda / pi0
    _plot_particle_ratios_vs_refmult(f, sums, results_post, ['3122'], ['111'],
                                     ytitle="#Lambda/#pi^{0}")
    # Xi / pi0
    _plot_particle_ratios_vs_refmult(f, sums, results_post, ['3312'], ['111'],
                                     ytitle="#Xi/#pi^{0}")
    # Omega / pi0
    _plot_particle_ratios_vs_refmult(f, sums, results_post, ['3334', '-3334'], ['111'],
                                     ytitle="#Omega/#pi^{0}")
    # K_ch / K0_S
    _plot_particle_ratios_vs_refmult(f, sums, results_post, ['321', '-321'], ['310'], scale=.5,
                                     ytitle="(K^{+}+K^{-}) / (2#timesK^{0}_{S})")
    # K0_S / Lambda
    _plot_particle_ratios_vs_refmult(f, sums, results_post, ['310'], ['-3122', '3122'],
                                     ytitle="K^{0}_{S} / #Lambda")
    # K0_S / Xi
    _plot_particle_ratios_vs_refmult(f, sums, results_post, ['310'], ['3312'],
                                     ytitle="K^{0}_{S} / #Xi")

    ######################################################################################
    # vs Est mult
    _plot_particle_ratios_vs_estmult(f, sums, results_post, ['321', '-321'], ['310'],
                                     scale=.5, ytitle="(K^{+} + K^{-}) / (2*K_{S}^{0})")


def _delete_results_dir(f, sums):
    # delete old result directory
    f.rm('MultEstimators/results_post')
    f.Write()


def _mk_results_dir(f, sums):
    f.mkdir('MultEstimators/results_post', recurse=True)
    for est_dir in sums:
        try:
            resdir = f.MultEstimators.results_post.mkdir(est_dir.GetName())
            resdir.Write()
        except:
            pass

if __name__ == "__main__":
    # go into batch mode
    ROOT.gROOT.SetBatch(True)

    log = log["/post"]  # set name of this script in logger
    log.info("IsBatch: {0}".format(ROOT.gROOT.IsBatch()))  # Results in "DEBUG:myapp] Hello"

    # Rebin multiplicity with factor:
    rebin_mult = 10
    ref_ests = ['EtaLt05', ]
    considered_ests = ['EtaLt05', 'EtaLt08', 'EtaLt15', 'Eta08_15', 'V0M', 'V0A', 'V0C']

    functions = [
        _delete_results_dir,
        _mk_results_dir,
    ]

    with root_open(sys.argv[1], 'update') as f:
        sums = f.MultEstimators.Sums
        for func in functions:
            func(f, sums)
        results_post = f.MultEstimators.results_post
        _make_event_counters(f, sums, results_post)
    with root_open(sys.argv[1], 'update') as f:
        sums = f.MultEstimators.Sums
        results_post = f.MultEstimators.results_post
        _make_dNdeta(f, sums, results_post)
        #     _make_correlation_plots(f, sums, results_post)
        _make_PNch_plots(f, sums, results_post)
        _make_mult_vs_pt_plots(f, sums, results_post)
    with root_open(sys.argv[1], 'update') as f:
        sums = f.MultEstimators.Sums
        results_post = f.MultEstimators.results_post
        _make_hists_vs_pt(f, sums, results_post)  # needs updated results_post!
        # _make_dNdeta_mb_ratio_plots(f, sums, results_post)
        _make_pid_ratio_plots(f, sums, results_post)

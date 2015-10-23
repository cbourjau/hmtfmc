"""
This is the entry point to run the post analysis. This file contains the steering part.
Ie. it does define which plots to make, but not how. The logic which extracts the data from the
primary generated "Sums" is in the post_data_extractor.py. These functions provide (lists of) plottables.
These plottables are then plotted with the plotting functions in plotting_util.py
"""

import sys
import string
import random

if len(sys.argv) != 2:
    print "Usage: python ./post.py path_to_root_file.root"
    quit()

from rootpy.io import root_open
from rootpy import asrootpy, ROOT, log, collection
from rootpy.plotting import Hist1D, Hist2D, Canvas, Legend

from post_data_extractors import get_dNdeta_in_mult_interval, get_identified_vs_mult,\
    get_Nch_edges_for_percentile_edges, get_nMPI_vs_Nch,\
    get_NchEst1_vs_NchEst2, get_PNch_vs_estmult,\
    get_meanpt_vs_estmult, get_pT_distribution, get_mean_nMPI
from post_utils import create_stack_pid_ratio_over_pt,\
    remap_x_values,\
    remove_zero_value_points, remove_non_mutual_points,\
    remove_points_with_equal_x, remove_points_with_x_err_gt_1NchRef

from roofi import Figure, Styles
from roofi.figure import get_color_generator


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
mean_mult_cutoff_factor = 4
perc_bins = [(1, 0.7), (.5, .4), (.1, .05), (0.001, 0.0)]


def get_est_dirs(sums):
    return (somedir for somedir in sums if somedir.GetName() in considered_ests)


def _plot_particle_ratios_vs_estmult(f, sums, results_post, pids1, pids2, scale=None, ytitle=''):
    ratio_vs_estmult_dir = 'MultEstimators/results_post/pid_ratios_vs_estmult'
    fig = Figure()
    if not ytitle:
        fig.ytitle = ", ".join(pids1) + " / " + ", ".join(pids2)
    else:
        fig.ytitle = ytitle

    for est_dir in get_est_dirs(sums):
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


def _get_graphs_particle_ratios_vs_refmult(f, sums, results_post, pids1, pids2, scale=None, ytitle=''):
    """
    Returns list of ratios of the two pid-lists (pids1/pids2) vs refmult.
    This function depends on the correlation histograms to be present in f
    """
    ratios = []
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
        ratio.title = make_estimator_title(est_dir.GetName())
        ratios.append(ratio)
    return ratios


def _make_event_counters(f, sums, results_post):
    log.info("Creating event counters")
    for est_dir in get_est_dirs(sums):
        results_est_dir = results_post.__getattr__(est_dir.GetName())
        corr = asrootpy(est_dir.FindObject("fcorr_thisNch_vs_refNch"))
        counter = asrootpy(corr.ProjectionX())
        counter.name = "event_counter"
        f.cd("MultEstimators/results_post/" + est_dir.GetName())
        results_est_dir.WriteTObject(counter)


def _make_dNdeta(f, sums, results_post):
    # Loop over all estimators in the Sums list:
    log.info("Creating dN/deta bin in multiplicity")
    for est_dir in get_est_dirs(sums):
        results_est_dir = results_post.Get(est_dir.GetName())
        h2d = asrootpy(est_dir.FindObject('feta_Nch'))
        event_counter = asrootpy(results_est_dir.Get("event_counter"))

        fig = Figure()
        fig_mb_ratio = Figure()
        fig.plot.palette = fig_mb_ratio.plot.palette = 'colorblind'
        fig.xtitle = fig_mb_ratio.xtitle = '#eta'
        fig.ytitle = '1/N dN_{ch}/d#eta'
        fig_mb_ratio.ytitle = 'dN/d#eta (1/MB)'
        fig.legend.title = fig_mb_ratio.legend.title = make_estimator_title(est_dir.GetName())
        fig.plot.ymin = fig_mb_ratio.plot.ymin = 0
        dNdeta_mb = get_dNdeta_in_mult_interval(h2d, event_counter, [0, 250])
        for nch_int, perc_int in zip(NCH_EDGES[est_dir.GetName()], perc_bins):
            title = "{}%-{}%".format(perc_int[0] * 100, perc_int[1] * 100)
            dNdeta_in_interval = get_dNdeta_in_mult_interval(h2d, event_counter, nch_int)
            fig.add_plottable(dNdeta_in_interval, legend_title=title)
            fig_mb_ratio.add_plottable(dNdeta_in_interval / dNdeta_mb, legend_title=title)
        # add MB as well:
        title = "MB"
        fig.add_plottable(dNdeta_mb, legend_title=title)
        path = results_est_dir.GetPath().split(":")[1]  # file.root:/internal/root/path
        fig.save_to_root_file(f, "dNdeta_summary", path=path)
        fig_mb_ratio.save_to_root_file(f, "dNdeta_MB_ratio_summary", path=path)


def _pt_distribution_ratios(f, sums, results_post):
    ###########################################################
    # Category 2 on TWiki
    # create particle ratio vs pT plots

    log.info("Computing histograms vs pt")

    # Loop over all estimators in the Sums list:
    for est_dir in get_est_dirs(results_post):
        dirname = 'MultEstimators/results_post/{}/pid_ratios/'.format(est_dir.GetName())

        # event_counter = asrootpy(getattr(results_post, est_dir.GetName()).event_counter)

        # mult_pt_dir = results_post.FindObject(est_dir.GetName()).Get("mult_pt")
        fig = Figure()
        fig.xtitle = 'p_{T} (GeV)'
        fig.plot.ymin = 0
        fig.plot.palette = 'colorblind'
        # fig.plot.palette_ncolors = len(nch_edges) - 1
        fig.legend.position = 'bl'

        mult_binned_pt_dists = {}
        mult_binned_pt_dists['proton'] = [
            get_pT_distribution(est_dir, [kANTIPROTON, kPROTON], nch_low, nch_up)
            for nch_low, nch_up in NCH_EDGES[est_dir.GetName()]
        ]
        mult_binned_pt_dists['pi_ch'] = [
            get_pT_distribution(est_dir, [kPIMINUS, kPIPLUS], nch_low, nch_up)
            for nch_low, nch_up in NCH_EDGES[est_dir.GetName()]
        ]
        mult_binned_pt_dists['xi'] = [
            get_pT_distribution(est_dir, [kANTIXI, kXI], nch_low, nch_up)
            for nch_low, nch_up in NCH_EDGES[est_dir.GetName()]
        ]
        mult_binned_pt_dists['omega'] = [
            get_pT_distribution(est_dir, [kOMEGAMINUS, kOMEGAPLUS], nch_low, nch_up)
            for nch_low, nch_up in NCH_EDGES[est_dir.GetName()]
        ]
        mult_binned_pt_dists['lambda'] = [
            get_pT_distribution(est_dir, [kANTILAMBDA, kLAMBDA], nch_low, nch_up)
            for nch_low, nch_up in NCH_EDGES[est_dir.GetName()]
        ]
        mult_binned_pt_dists['k0s'] = [
            get_pT_distribution(est_dir, [kK0S], nch_low, nch_up)
            for nch_low, nch_up in NCH_EDGES[est_dir.GetName()]
        ]
        mult_binned_pt_dists['pi0'] = [
            get_pT_distribution(est_dir, [kPI0], nch_low, nch_up)
            for nch_low, nch_up in NCH_EDGES[est_dir.GetName()]
        ]
        perc_titles = ["{}%-{}%".format(perc_bin[0] * 100, perc_bin[1] * 100) for perc_bin in perc_bins]

        fig.delete_plottables()
        name = "proton_over_pich__vs__pt"
        fig.ytitle = "(p+#bar{p})/#pi^{+-}"
        fig.plot.ymax = .3
        fig.legend.title = make_estimator_title(est_dir.GetName())
        [
            fig.add_plottable(h1 / h2, legend_title=title)
            for h1, h2, title in zip(mult_binned_pt_dists['proton'], mult_binned_pt_dists['pi_ch'], perc_titles)
        ]
        fig.save_to_root_file(f, name, dirname)

        fig.delete_plottables()
        name = "Xi_over_pich__vs__pt"
        fig.plot.ymax = .06
        fig.legend.position = 'tl'
        fig.ytitle = "#Xi/#pi^{+-}"
        fig.legend.title = make_estimator_title(est_dir.GetName())
        [
            fig.add_plottable(h1 / h2, legend_title=title)
            for h1, h2, title in zip(mult_binned_pt_dists['xi'], mult_binned_pt_dists['pi_ch'], perc_titles)
        ]
        fig.save_to_root_file(f, name, dirname)

        fig.delete_plottables()
        name = "OmegaCh_over_pich__vs__pt"
        fig.plot.ymax = .005
        fig.legend.position = 'tl'
        fig.ytitle = "#Omega_{ch}/#pi^{+-} "
        fig.legend.title = make_estimator_title(est_dir.GetName())
        [
            fig.add_plottable(h1 / h2, legend_title=title)
            for h1, h2, title in zip(mult_binned_pt_dists['omega'], mult_binned_pt_dists['pi_ch'], perc_titles)
        ]
        fig.save_to_root_file(f, name, dirname)

        # Ratios to pi0
        fig.delete_plottables()
        name = "pich_over_pi0__vs__pt"
        fig.plot.ymax = 2.5
        fig.legend.position = 'bl'
        fig.ytitle = "#pi^{+-}/#pi^{0}"
        fig.legend.title = make_estimator_title(est_dir.GetName())
        [
            fig.add_plottable(h1 / h2, legend_title=title)
            for h1, h2, title in zip(mult_binned_pt_dists['pi_ch'], mult_binned_pt_dists['pi0'], perc_titles)
        ]
        fig.save_to_root_file(f, name, dirname)

        fig.delete_plottables()
        name = "proton_over_pi0__vs__pt"
        fig.plot.ymax = 1
        fig.legend.position = 'tr'
        fig.ytitle = "p/#pi^{0}"
        fig.legend.title = make_estimator_title(est_dir.GetName())
        [
            fig.add_plottable(h1 / h2, legend_title=title)
            for h1, h2, title in zip(mult_binned_pt_dists['proton'], mult_binned_pt_dists['pi0'], perc_titles)
        ]
        fig.save_to_root_file(f, name, dirname)

        fig.delete_plottables()
        name = "K0S_over_pi0__vs__pt"
        fig.plot.ymax = 1.4
        fig.legend.position = 'tl'
        fig.ytitle = "K^{0}_{S}/#pi^{0}"
        fig.legend.title = make_estimator_title(est_dir.GetName())
        [
            fig.add_plottable(h1 / h2, legend_title=title)
            for h1, h2, title in zip(mult_binned_pt_dists['k0s'], mult_binned_pt_dists['pi0'], perc_titles)
        ]
        fig.save_to_root_file(f, name, dirname)

        fig.delete_plottables()
        name = "Lambda_over_pi0__vs__pt"
        fig.plot.ymax = .9
        fig.legend.position = 'tl'
        fig.ytitle = "#Lambda/#pi^{0}"
        fig.legend.title = make_estimator_title(est_dir.GetName())
        [
            fig.add_plottable(h1 / h2, legend_title=title)
            for h1, h2, title in zip(mult_binned_pt_dists['lambda'], mult_binned_pt_dists['pi0'], perc_titles)
        ]
        fig.save_to_root_file(f, name, dirname)

        fig.delete_plottables()
        name = "Xi_over_pi0__vs__pt"
        fig.plot.ymax = .08
        fig.legend.position = 'tl'
        fig.ytitle = "#Xi/#pi^{0}"
        fig.legend.title = make_estimator_title(est_dir.GetName())
        [
            fig.add_plottable(h1 / h2, legend_title=title)
            for h1, h2, title in zip(mult_binned_pt_dists['xi'], mult_binned_pt_dists['pi0'], perc_titles)
        ]
        fig.save_to_root_file(f, name, dirname)

        fig.delete_plottables()
        name = "OmegaCh_over_pi0__vs__pt"
        fig.plot.ymax = .005
        fig.legend.position = 'tl'
        fig.ytitle = "#Omega_{ch}/#pi^{0}"
        fig.legend.title = make_estimator_title(est_dir.GetName())
        [
            fig.add_plottable(h1 / h2, legend_title=title)
            for h1, h2, title in zip(mult_binned_pt_dists['omega'], mult_binned_pt_dists['pi0'], perc_titles)
        ]
        fig.save_to_root_file(f, name, dirname)

        # Ratios to K0S
        fig.delete_plottables()
        name = "proton_over_K0S__vs__pt"
        fig.plot.ymax = 2.6
        fig.legend.position = 'tr'
        fig.ytitle = "p/K^{0}_{S}"
        fig.legend.title = make_estimator_title(est_dir.GetName())
        [
            fig.add_plottable(h1 / h2, legend_title=title)
            for h1, h2, title in zip(mult_binned_pt_dists['proton'], mult_binned_pt_dists['k0s'], perc_titles)
        ]
        fig.save_to_root_file(f, name, dirname)

        fig.delete_plottables()
        name = "Lambda_over_K0S__vs__pt"
        fig.plot.ymax = 1
        fig.legend.position = 'bl'
        fig.ytitle = "#Lambda/K^{0}_{S}"
        fig.legend.title = make_estimator_title(est_dir.GetName())
        [
            fig.add_plottable(h1 / h2, legend_title=title)
            for h1, h2, title in zip(mult_binned_pt_dists['lambda'], mult_binned_pt_dists['k0s'], perc_titles)
        ]
        fig.save_to_root_file(f, name, dirname)

        fig.delete_plottables()
        name = "Xi_over_K0S__vs__pt"
        fig.plot.ymax = .2
        fig.legend.position = 'tl'
        fig.ytitle = "#Xi/K^{0}_{S}"
        fig.legend.title = make_estimator_title(est_dir.GetName())
        [
            fig.add_plottable(h1 / h2, legend_title=title)
            for h1, h2, title in zip(mult_binned_pt_dists['xi'], mult_binned_pt_dists['k0s'], perc_titles)
        ]
        fig.save_to_root_file(f, name, dirname)

        fig.delete_plottables()
        name = "OmegaCh_over_K0S__vs__pt"
        fig.plot.ymax = .012
        fig.legend.position = 'tl'
        fig.ytitle = "#Omega_{ch}/K^{0}_{S}"
        fig.legend.title = make_estimator_title(est_dir.GetName())
        [
            fig.add_plottable(h1 / h2, legend_title=title)
            for h1, h2, title in zip(mult_binned_pt_dists['omega'], mult_binned_pt_dists['k0s'], perc_titles)
        ]
        fig.save_to_root_file(f, name, dirname)


def _make_PNch_plots(f, sums, results_post):
    log.info("Creating P(Nch) summary plot")
    summary_fig = Figure()
    summary_fig.xtitle = "N_{ch}^{est}"
    summary_fig.ytitle = "P(N_{ch}^{est})"
    summary_fig.legend.position = 'tr'

    for est_dir in get_est_dirs(sums):
        est_name = est_dir.GetName()
        h_tmp = get_PNch_vs_estmult(sums, est_name)
        if h_tmp.Integral() > 0:
            h_tmp.Scale(1.0 / h_tmp.Integral())
            summary_fig.add_plottable(h_tmp, make_estimator_title(est_name))

    summary_fig.plot.logy = True
    path = results_post.GetPath().split(":")[1]  # file.root:/internal/root/path
    summary_fig.save_to_root_file(f, "PNch_summary", path=path)

    log.info("Creating P(Nch_est) and P(Nch_refest) histograms")
    mult_bin_size = 10
    for ref_est_name in ref_ests:
        for res_est_dir in get_est_dirs(results_post):
            est_name = res_est_dir.GetName()
            # Figure properties:
            fig_vs_estmult = Figure()
            fig_vs_refmult = Figure()
            fig_vs_estmult.plot.logy = True
            fig_vs_refmult.plot.logy = True

            fig_vs_estmult.legend.position = 'tr'
            fig_vs_refmult.legend.position = 'tr'

            fig_vs_estmult.xtitle = "N_{{ch}}^{{{}}}".format(est_name)
            fig_vs_refmult.xtitle = "N_{{ch}}^{{{}}}".format(ref_est_name)

            fig_vs_estmult.ytitle = "P(N_{{ch}}^{{{}}})".format(est_name)
            fig_vs_refmult.ytitle = "P(N_{{ch}}^{{{}}})".format(ref_est_name)

            corr_hist = get_NchEst1_vs_NchEst2(sums, ref_est_name, est_name)

            # logic when dealing with fixed bins given in Nch:
            # ------------------------------------------------
            # mean_nch_est = corr_hist.GetMean(1)  # mean of x axis
            # nch_max = corr_hist.xaxis.GetNbins()
            # nch_cutoff = mean_nch_est * mean_mult_cutoff_factor
            # nch_bins = [(low, low + mult_bin_size) for low in range(0, int(nch_cutoff), mult_bin_size)]
            # # a large last bin covering the rest:
            # nch_bins += [(nch_bins[-1][2], nch_max)]
            # legend_tmpl = "{} < N_{ch} < {}"
            # logic when dealing with percentile bins:
            # ----------------------------------------
            # event_counter_est = asrootpy(getattr(res_est_dir, "event_counter"))
            nch_bins_est = NCH_EDGES[est_name]
            nch_bins_ref = NCH_EDGES[ref_est_name]
            legend_tmpl = "{}% - {}%"
            fig_vs_estmult.legend.title = "Selected in {}".format(make_estimator_title(ref_est_name))
            fig_vs_refmult.legend.title = "Selected in {}".format(make_estimator_title(est_name))
            # WARNING: the following needs tweeking when going back to fixed N_ch bins!
            for nch_bin, perc_bin in zip(nch_bins_ref, perc_bins):
                # vs est_mult:
                corr_hist.xaxis.SetRange(0, 0)  # reset x axis
                corr_hist.yaxis.SetRange(nch_bin[0] + 1, nch_bin[1] + 1)  # mind the crackpot binning!
                h_vs_est = asrootpy(corr_hist.ProjectionX(gen_random_name()))
                if h_vs_est.Integral() > 0:
                    h_vs_est.Scale(1.0 / h_vs_est.Integral())
                    fig_vs_estmult.add_plottable(h_vs_est, legend_tmpl.format(perc_bin[0] * 100, perc_bin[1] * 100))
                else:
                    log.info("No charged particles in {}*100 percentile bin of estimator {}. This should not happen".
                             format(perc_bin, ref_est_name))
            for nch_bin, perc_bin in zip(nch_bins_est, perc_bins):
                # vs ref_mult:
                corr_hist.yaxis.SetRange(0, 0)  # reset y axis
                corr_hist.xaxis.SetRange(*nch_bin)
                h_vs_ref = asrootpy(corr_hist.ProjectionY(gen_random_name()))
                if h_vs_ref.Integral() > 0:
                    h_vs_ref.Scale(1.0 / h_vs_ref.Integral())
                    fig_vs_refmult.add_plottable(h_vs_ref, legend_tmpl.format(perc_bin[0] * 100, perc_bin[1] * 100))
                else:
                    log.info(
                        "No charged particles in {}*100 percentile bin of estimator {}. This should not happen".
                        format(perc_bin, est_name))

            path = results_post.GetPath().split(":")[1] + "/" + est_name  # file.root:/internal/root/path
            # vs est_mult
            fig_vs_estmult.save_to_root_file(f, "PNchEst_binned_in_Nch{}".format(ref_est_name), path)
            # vs est_mult
            fig_vs_refmult.save_to_root_file(f, "PNch{}_binned_in_NchEst".format(ref_est_name), path)


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
    ratios_dir = 'MultEstimators/results_post/pid_ratios_vs_refmult'
    fig = Figure()
    fig.plot.ncolors = len(considered_ests)
    fig.xtitle = "N_{ch}|_{" + make_estimator_title('EtaLt05') + "}"
    fig.plot.xmin = 0
    fig.plot.xmax = 40

    # Proton / pi_ch

    pids1, pids2 = ['-2212', '2212'], ['-211', '211']
    fig.ytitle = "p/#pi^{+-}"
    fig.plot.ymin, fig.plot.ymax = 0.04, 0.13
    fig.delete_plottables()
    graphs = _get_graphs_particle_ratios_vs_refmult(f, sums, results_post, pids1, pids2)
    [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
    name = "_".join(pids1) + "_div_" + "_".join(pids2)
    fig.save_to_root_file(f, name, ratios_dir)

    # K / pi_ch
    pids1, pids2 = ['310', '321', '-321'], ['-211', '211']
    fig.ytitle = "K^{*}/#pi^{+-}"
    fig.plot.ymin, fig.plot.ymax = 0.09, 0.30
    fig.delete_plottables()
    graphs = _get_graphs_particle_ratios_vs_refmult(f, sums, results_post, pids1, pids2)
    [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
    name = "_".join(pids1) + "_div_" + "_".join(pids2)
    fig.save_to_root_file(f, name, ratios_dir)

    # Lambda / pi_ch
    pids1, pids2 = ['3122'], ['-211', '211']
    fig.ytitle = "#Lambda / #pi^{+-}"
    fig.plot.ymin, fig.plot.ymax = 0.005, 0.035
    fig.delete_plottables()
    graphs = _get_graphs_particle_ratios_vs_refmult(f, sums, results_post, pids1, pids2)
    [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
    name = "_".join(pids1) + "_div_" + "_".join(pids2)
    fig.save_to_root_file(f, name, ratios_dir)

    # Xi / pi_ch
    pids1, pids2 = ['3312'], ['-211', '211']
    fig.ytitle = "#Xi / #pi^{+-}"
    fig.plot.ymin, fig.plot.ymax = 0.0004, 0.002
    fig.delete_plottables()
    graphs = _get_graphs_particle_ratios_vs_refmult(f, sums, results_post, pids1, pids2)
    [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
    name = "_".join(pids1) + "_div_" + "_".join(pids2)
    fig.save_to_root_file(f, name, ratios_dir)

    # Omega / pi_ch
    pids1, pids2 = ['3334', '-3334'], ['-211', '211']
    fig.ytitle = "#Omega / #pi^{+-}"
    fig.plot.ymin, fig.plot.ymax = 0.00001, 0.00011
    fig.delete_plottables()
    graphs = _get_graphs_particle_ratios_vs_refmult(f, sums, results_post, pids1, pids2)
    [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
    name = "_".join(pids1) + "_div_" + "_".join(pids2)
    fig.save_to_root_file(f, name, ratios_dir)

    # pi_ch/pi0
    pids1, pids2 = ['-211', '211'], ['111']
    fig.ytitle = "#pi^{+-}/#pi^{0}"
    fig.plot.ymin, fig.plot.ymax = 1.5, 2.2
    fig.delete_plottables()
    graphs = _get_graphs_particle_ratios_vs_refmult(f, sums, results_post, pids1, pids2)
    [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
    name = "_".join(pids1) + "_div_" + "_".join(pids2)
    fig.save_to_root_file(f, name, ratios_dir)

    # proton / pi0
    pids1, pids2 = ['-2212', '2212'], ['111']
    fig.ytitle = "p/#pi^{0}"
    fig.plot.ymin, fig.plot.ymax = 0.09, 0.13
    fig.legend.position = 'bl'
    fig.delete_plottables()
    graphs = _get_graphs_particle_ratios_vs_refmult(f, sums, results_post, pids1, pids2)
    [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
    name = "_".join(pids1) + "_div_" + "_".join(pids2)
    fig.save_to_root_file(f, name, ratios_dir)
    fig.legend.position = 'tl'

    # K / pi0
    pids1, pids2 = ['310', '321', '-321'], ['111']
    fig.ytitle = "K^{*}/#pi^{0}"
    fig.plot.ymin, fig.plot.ymax = 0.15, 0.34
    fig.legend.position = 'bl'
    fig.delete_plottables()
    graphs = _get_graphs_particle_ratios_vs_refmult(f, sums, results_post, pids1, pids2)
    [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
    name = "_".join(pids1) + "_div_" + "_".join(pids2)
    fig.save_to_root_file(f, name, ratios_dir)
    fig.legend.position = 'tl'

    # Lambda / pi0
    pids1, pids2 = ['3122'], ['111']
    fig.ytitle = "#Lambda/#pi^{0}"
    fig.plot.ymin, fig.plot.ymax = 0.014, 0.036
    fig.delete_plottables()
    graphs = _get_graphs_particle_ratios_vs_refmult(f, sums, results_post, pids1, pids2)
    [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
    name = "_".join(pids1) + "_div_" + "_".join(pids2)
    fig.save_to_root_file(f, name, ratios_dir)

    # Xi / pi0
    pids1, pids2 = ['3312'], ['111']
    fig.ytitle = "#Xi/#pi^{0}"
    fig.plot.ymin, fig.plot.ymax = 0.0010, 0.0018
    fig.legend.position = 'bl'
    fig.delete_plottables()
    graphs = _get_graphs_particle_ratios_vs_refmult(f, sums, results_post, pids1, pids2)
    [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
    name = "_".join(pids1) + "_div_" + "_".join(pids2)
    fig.save_to_root_file(f, name, ratios_dir)
    fig.legend.position = 'tl'

    # Omega / pi0
    pids1, pids2 = ['3334', '-3334'], ['111']
    fig.ytitle = "#Omega/#pi^{0}"
    fig.legend.position = 'bl'
    fig.plot.ymin, fig.plot.ymax = 0.00002, 0.00010
    fig.delete_plottables()
    graphs = _get_graphs_particle_ratios_vs_refmult(f, sums, results_post, pids1, pids2)
    [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
    name = "_".join(pids1) + "_div_" + "_".join(pids2)
    fig.save_to_root_file(f, name, ratios_dir)
    fig.legend.position = 'tl'

    # K_ch / K0_S
    pids1, pids2 = ['321', '-321'], ['310']
    fig.ytitle = "(K^{+}+K^{-}) / (2#timesK^{0}_{S})"
    fig.plot.ymin, fig.plot.ymax = 0.4, 1.3
    fig.legend.position = 'bl'
    fig.delete_plottables()
    graphs = _get_graphs_particle_ratios_vs_refmult(f, sums, results_post, pids1, pids2, scale=.5)
    [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
    name = "_".join(pids1) + "_div_" + "_".join(pids2)
    fig.save_to_root_file(f, name, ratios_dir)
    fig.legend.position = 'tl'

    # K0_S / Lambda
    pids1, pids2 = ['310'], ['-3122', '3122']
    fig.ytitle = "K^{0}_{S} / #Lambda"
    fig.plot.ymin, fig.plot.ymax = 1.9, 3.7
    fig.delete_plottables()
    graphs = _get_graphs_particle_ratios_vs_refmult(f, sums, results_post, pids1, pids2)
    [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
    name = "_".join(pids1) + "_div_" + "_".join(pids2)
    fig.save_to_root_file(f, name, ratios_dir)

    # K0_S / Xi
    pids1, pids2 = ['310'], ['3312']
    fig.ytitle = "K^{0}_{S} / #Xi"
    fig.plot.ymin, fig.plot.ymax = 57, 78
    fig.delete_plottables()
    graphs = _get_graphs_particle_ratios_vs_refmult(f, sums, results_post, pids1, pids2)
    [fig.add_plottable(g, legend_title=g.GetTitle()) for g in graphs]
    name = "_".join(pids1) + "_div_" + "_".join(pids2)
    fig.save_to_root_file(f, name, ratios_dir)


    # ######################################################################################
    # # vs Est mult
    # _plot_particle_ratios_vs_estmult(f, sums, results_post, ['321', '-321'], ['310'],
    #                                  scale=.5, fig.ytitle = "(K^{+} + K^{-}) / (2*K_{S}^{0})")


def _plot_meanpt_vs_ref_mult_for_pids(f, sums, results_post):
    log.info("Creating mean pT plots")
    for sums_est_dir, res_est_dir in zip(get_est_dirs(sums), get_est_dirs(results_post)):
        if sums_est_dir.GetName() != res_est_dir.GetName():
            raise IndexError("Order of estimator dirs is different in sums and results_post")
        res_dir_str = "MultEstimators/results_post/" + res_est_dir.GetName()
        corr_hist = asrootpy(sums_est_dir.FindObject("fcorr_thisNch_vs_refNch"))
        fig = Figure()
        fig.plot.palette = 'root'
        fig.plot.ncolors = 7
        fig.plot.xmin = 0
        fig.plot.xmax = 25
        fig.plot.ymin = 0.3
        fig.plot.ymax = 1.9
        fig.ytitle = "<p_{T}>"
        fig.xtitle = "N_{ch}|_{|#eta|<0.5}"
        fig.legend.title = make_estimator_title(sums_est_dir.GetName())
        graphs = []
        graphs.append(remap_x_values(get_meanpt_vs_estmult(res_est_dir, [kPI0, kPIMINUS, kPIPLUS]), corr_hist))
        graphs[-1].title = "#pi"
        graphs.append(remap_x_values(get_meanpt_vs_estmult(res_est_dir, [kKMINUS, kKPLUS]), corr_hist))
        graphs[-1].title = "K^{#pm}"
        graphs.append(remap_x_values(get_meanpt_vs_estmult(res_est_dir, [kPROTON, kANTIPROTON]), corr_hist))
        graphs[-1].title = "p"
        graphs.append(remap_x_values(get_meanpt_vs_estmult(res_est_dir, [kK0S]), corr_hist))
        graphs[-1].title = "K^{0}_{S}"
        graphs.append(remap_x_values(get_meanpt_vs_estmult(res_est_dir, [kLAMBDA, kANTILAMBDA]), corr_hist))
        graphs[-1].title = "#Lambda"
        graphs.append(remap_x_values(get_meanpt_vs_estmult(res_est_dir, [kXI, kANTIXI]), corr_hist))
        graphs[-1].title = "#Xi"
        graphs.append(remap_x_values(get_meanpt_vs_estmult(res_est_dir, [kOMEGAMINUS, kOMEGAPLUS]), corr_hist))
        graphs[-1].title = "#Omega"
        # sanitize graphs:
        for g in graphs:
            remove_zero_value_points(g)
            remove_points_with_x_err_gt_1NchRef(g)
            remove_points_with_equal_x(g)
        [fig.add_plottable(g, g.title) for g in graphs]
        fig.save_to_root_file(f, "mean_pt", res_dir_str)


def _plot_event_counter_with_shaded_perc_areas(f, results_post):
    log.info("Broken: Root sucks! Creating shaded event counter with percentile regions")
    return
    for est_dir in get_est_dirs(results_post):
        event_counter = asrootpy(getattr(est_dir, "event_counter"))
        perc_edges = [1, .6, .4, .2, .1, .05, .025, 0.012, 0]
        nch_edges = get_Nch_edges_for_percentile_edges(perc_edges, event_counter)
        c = Canvas(name="event_counter_with_perc")
        leg = Legend(len(nch_edges) - 1)
        copies = []
        colors = get_color_generator(ncolors=10)
        # Draw the hist once
        event_counter.Draw()
        for nch_low, nch_up in zip(nch_edges[:-1], nch_edges[1:]):
            copies.append(event_counter.Clone(gen_random_name()))
            copies[-1].xaxis.SetRangeUser(nch_low, nch_up)
            copies[-1].SetFillStyle(1001)
            copies[-1].color = next(colors)
            copies[-1].xaxis.title = "N_{ch}"
            copies[-1].yaxis.title = "counts"
            leg.AddEntry(copies[-1], "{}-{}%".format(str(nch_low), str(nch_up)))
            copies[-1].Draw('sameHist')
            break
        leg.Draw()
        est_dir.cd()
        c.Write()


def _plot_dNdpT(f, sums, results_post):
    log.info("1/N_evts  dN_ch/dpT plots")
    for sums_est_dir, res_est_dir in zip(get_est_dirs(sums), get_est_dirs(results_post)):
        if sums_est_dir.GetName() != res_est_dir.GetName():
            raise IndexError("Order of estimator dirs is different in sums and results_post")
        res_dir_str = "MultEstimators/results_post/" + res_est_dir.GetName()
        fig = Figure()
        fig.plot.palette = 'colorblind'
        # fig.plot.ncolors = 5
        fig.legend.title = make_estimator_title(sums_est_dir.GetName())
        fig.legend.position = 'tr'
        fig.ytitle = "#frac{1}{N_{evts}}dN/dp_{T} MB"
        fig.xtitle = "p_{T} (GeV)"
        fig.plot.logy = True
        hists = []
        nch_low, nch_up = 0, 250
        charged_particles = [kPIMINUS, kPIPLUS, kKMINUS, kKPLUS, kPROTON, kANTIPROTON,
                             kLAMBDA, kANTILAMBDA, kXI, kANTIXI, kOMEGAMINUS, kOMEGAPLUS]
        hists.append(get_pT_distribution(res_est_dir, charged_particles, nch_low, nch_up, normalized=False))
        hists[-1].title = "MB"

        event_counter = asrootpy(getattr(res_est_dir, "event_counter"))
        for perc_bin_up, perc_bin_low in perc_bins:
            nch_low, nch_up = get_Nch_edges_for_percentile_edges([perc_bin_up, perc_bin_low], event_counter)
            hists.append(get_pT_distribution(res_est_dir, charged_particles, nch_low, nch_up, normalized=False))
            hists[-1].title = "{}%-{}%".format(perc_bin_up * 100, perc_bin_low * 100)
        # scale by bin width
        [h.Scale(1, "width") for h in hists]

        [fig.add_plottable(p, p.title) for p in hists]
        fig.legend.title = "#pi^{#pm}, K^{#pm}, p, #Lambda, #Xi, #Omega"
        fig.save_to_root_file(f, "dNdpT_MB", res_dir_str)


def _plot_pT_HM_div_pt_MB(f, sums, results_post, scale_nMPI):
    log.info("Plot dN/dpT ratios scaled with nMPI")
    for sums_est_dir, res_est_dir in zip(get_est_dirs(sums), get_est_dirs(results_post)):
        if sums_est_dir.GetName() != res_est_dir.GetName():
            raise IndexError("Order of estimator dirs is different in sums and results_post")
        res_dir_str = "MultEstimators/results_post/" + res_est_dir.GetName()
        fig = Figure()
        fig.plot.palette = 'root'
        fig.plot.ncolors = 7
        fig.ytitle = ("#left[ #frac{dN^{HM}}{dp_{T}} / #frac{dN^{MB}}{dp_{T}} #right] "
                      "#times #left[ #frac{<N_{MPI}^{MB}>}{<N_{MPI}^{HM}>} #right]")
        fig.xtitle = "p_{T} (GeV)"
        fig.legend.title = make_estimator_title(sums_est_dir.GetName())
        event_counter = asrootpy(getattr(res_est_dir, "event_counter"))
        charged_particles = [kPIMINUS, kPIPLUS, kKMINUS, kKPLUS, kPROTON, kANTIPROTON,
                             kLAMBDA, kANTILAMBDA, kXI, kANTIXI, kOMEGAMINUS, kOMEGAPLUS]

        # get the MB distribution which will be used to devide the nch-binned distributions
        nch_low_mb, nch_up_mb = 0, 250
        pt_dist_mb = get_pT_distribution(res_est_dir, charged_particles, nch_low_mb, nch_up_mb, normalized=False)
        mean_nmpi_mb = get_mean_nMPI(sums_est_dir, nch_low_mb, nch_up_mb)

        event_counter = asrootpy(getattr(res_est_dir, "event_counter"))
        for perc_bin_up, perc_bin_low in perc_bins:
            nch_low, nch_up = get_Nch_edges_for_percentile_edges([perc_bin_up, perc_bin_low], event_counter)
            # get the pt distribution in this Nch interval
            pt_dist_in_interval = get_pT_distribution(res_est_dir, charged_particles,
                                                      nch_low, nch_up, normalized=False)
            title = "{}%-{}%".format(perc_bin_up * 100, perc_bin_low * 100)
            if scale_nMPI:
                mean_nmpi_hm = get_mean_nMPI(sums_est_dir, nch_low, nch_up)
                fig.add_plottable((pt_dist_in_interval / pt_dist_mb) * (mean_nmpi_mb / mean_nmpi_hm), title)
                name = "pt_hm_div_pt_mb_scaled_nMPI"
            else:
                fig.add_plottable((pt_dist_in_interval / pt_dist_mb), title)
                name = "pt_hm_div_pt_mb"
        fig.save_to_root_file(f, name, res_dir_str)


def _plot_nMPI_vs_Nch(f, sums, results_post):
    log.info("Creating nMPI(Nch) summary plot")
    summary_fig = Figure()
    summary_fig.xtitle = "N_{ch}^{est}"
    summary_fig.ytitle = "<N_{MPI}>"
    summary_fig.plot.palette = 'root'
    summary_fig.legend.position = 'br'
    summary_fig.plot.logy = True

    for est_dir in get_est_dirs(sums):
        h_tmp = get_nMPI_vs_Nch(est_dir)
        summary_fig.add_plottable(h_tmp, make_estimator_title(est_dir.GetName()))

    path = results_post.GetPath().split(":")[1]  # file.root:/internal/root/path
    summary_fig.save_to_root_file(f, "nMPI_summary", path=path)


def _delete_results_dir(f, sums):
    # delete old result directory
    f.rm('MultEstimators/results_post')
    f.Write()


def _mk_results_dir(f, sums):
    f.mkdir('MultEstimators/results_post', recurse=True)
    for est_dir in get_est_dirs(sums):
        try:
            resdir = f.MultEstimators.results_post.mkdir(est_dir.GetName())
            resdir.Write()
        except:
            pass

if __name__ == "__main__":
    # go into batch mode
    ROOT.gROOT.SetBatch(True)
    # set default style for all plots
    Figure.style = Styles.Presentation_half

    log = log["/post"]  # set name of this script in logger
    log.info("IsBatch: {0}".format(ROOT.gROOT.IsBatch()))  # Results in "DEBUG:myapp] Hello"

    # Rebin multiplicity with factor:
    rebin_mult = 10
    ref_ests = ['EtaLt05', ]
    considered_ests = ['EtaLt05', 'EtaLt08', 'EtaLt15', 'Eta08_15', 'V0M', 'V0A', 'V0C', 'ZDC']

    reset_functions = [
        _delete_results_dir,
        _mk_results_dir,
    ]
    plotting_functions = [
        _make_dNdeta,
        _make_PNch_plots,
        _plot_nMPI_vs_Nch,
        _make_mult_vs_pt_plots,
        _plot_meanpt_vs_ref_mult_for_pids,
        _make_hists_vs_pt,  # needs updated results_post!
        # _make_dNdeta_mb_ratio_plots,
        _make_pid_ratio_plots,
        _plot_dNdpT,
    ]

    def delete_lists(l):
        """Recursivley delete the lists to free memory"""
        for obj in l:
            if isinstance(obj, collection.List):
                obj.Delete()
        l.Delete()

    for func in reset_functions:
        with root_open(sys.argv[1], 'update') as f:
            sums = f.MultEstimators.Sums
            func(f, sums)
            delete_lists(sums)

    # Find Nch edges from percentile for all estimators:
    # --------------------------------------------------
    # Do it globally here because root is a piece of shit and it cannot open two event_counters
    # from different directories at the same time without corrupting the file. Fuck you ROOT!
    NCH_EDGES = {}
    with root_open(sys.argv[1], 'update') as f:
        sums = f.MultEstimators.Sums
        results_post = f.MultEstimators.results_post
        _make_event_counters(f, sums, results_post)
        delete_lists(sums)
    with root_open(sys.argv[1], 'update') as f:
        sums = f.MultEstimators.Sums
        results_post = f.MultEstimators.results_post
        for est_dir in get_est_dirs(results_post):
            event_counter = est_dir.event_counter
            NCH_EDGES[est_dir.GetName()] = [
                get_Nch_edges_for_percentile_edges(perc_bin, event_counter) for perc_bin in perc_bins]
        delete_lists(sums)

    for func in plotting_functions:
        with root_open(sys.argv[1], 'update') as f:
            sums = f.MultEstimators.Sums
            results_post = f.MultEstimators.results_post
            func(f, sums, results_post)
            delete_lists(sums)
    # the following two need extra parameters
    with root_open(sys.argv[1], 'update') as f:
        sums = f.MultEstimators.Sums
        results_post = f.MultEstimators.results_post
        _plot_pT_HM_div_pt_MB(f, sums, results_post, scale_nMPI=False)
    with root_open(sys.argv[1], 'update') as f:
        sums = f.MultEstimators.Sums
        results_post = f.MultEstimators.results_post
        _plot_pT_HM_div_pt_MB(f, sums, results_post, scale_nMPI=True)

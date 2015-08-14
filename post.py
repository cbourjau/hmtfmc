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
    get_NchEst1_vs_NchEst2, get_PNch_vs_estmult
from post_utils import create_stack_pid_ratio_over_pt,\
    remap_x_values,\
    plot_list_of_plottables, remove_zero_value_points, remove_non_mutual_points,\
    remove_points_with_equal_x, remove_points_with_x_err_gt_1NchRef


def gen_random_name():
    """Generate a random name for temp hists"""
    return ''.join(random.choice(string.ascii_letters + string.digits) for _ in range(25))


def get_color_generator(ncolors):
    """Returns a generator for n colors"""
    # return (int(800 + idx*100.0/ncolors) + 1 for idx in xrange(ncolors))

    # generated with sns.palplot(sns.color_palette("colorblind", 10))
    colorblind_colors = [(0.0, 0.4470588235294118, 0.6980392156862745),
                         (0.0, 0.6196078431372549, 0.45098039215686275),
                         (0.8352941176470589, 0.3686274509803922, 0.0),
                         (0.8, 0.4745098039215686, 0.6549019607843137),
                         (0.9411764705882353, 0.8941176470588236, 0.25882352941176473),
                         (0.33725490196078434, 0.7058823529411765, 0.9137254901960784),
                         (0.0, 0.4470588235294118, 0.6980392156862745),
                         (0.0, 0.6196078431372549, 0.45098039215686275),
                         (0.8352941176470589, 0.3686274509803922, 0.0),
                         (0.8, 0.4745098039215686, 0.6549019607843137)]
    # return iter(colorblind_colors)
    set2 = [(0.40000000596046448, 0.7607843279838562, 0.64705884456634521),
            (0.98131487965583808, 0.55538641635109398, 0.38740485135246722),
            (0.55432528607985565, 0.62711267120697922, 0.79595541393055635),
            (0.90311419262605563, 0.54185316071790801, 0.76495195557089413),
            (0.65371782148585622, 0.84708959004458262, 0.32827375098770734),
            (0.9986312957370983, 0.85096502233954041, 0.18488274134841617),
            (0.89573241682613591, 0.76784315109252932, 0.58182240093455595),
            (0.70196080207824707, 0.70196080207824707, 0.70196080207824707),
            (0.40000000596046448, 0.7607843279838562, 0.64705884456634521),
            (0.98131487965583808, 0.55538641635109398, 0.38740485135246722)]
    return iter(set2)


def make_estimator_title(name):
    if name == 'EtaLt05':
        return '|#eta| #leq 0.5'
    elif name == 'EtaLt08':
        return '|#eta| #leq 0.8'
    elif name == 'EtaLt15':
        return '|#eta| #leq 1.5'
    elif name == 'Eta08_15':
        return '0.8 #leq |#eta| #leq 1.5'
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


def _plot_particle_ratios_vs_estmult(f, sums, results_post, pids1, pids2, scale=None, ytitle=''):
    stack = list()
    colors = get_color_generator(10)
    if not ytitle:
        ytitle = ", ".join(pids1) + " / " + ", ".join(pids2)

    for est_dir in sums:
        if est_dir.GetName() not in ['EtaLt05', 'EtaLt08', 'EtaLt15', 'Eta08_15', 'V0M']:
            continue
        h3d = asrootpy(est_dir.FindObject("fNch_pT_pid"))
        pids1hists = [get_identified_vs_mult(h3d, pdg) for pdg in pids1]
        pids2hists = [get_identified_vs_mult(h3d, pdg) for pdg in pids2]

        pids1_px = sum(pids1hists)
        pids2_px = sum(pids2hists)
        ratio1d = pids1_px / pids2_px
        ratio1d.SetTitle(make_estimator_title(est_dir.GetName()))
        ratio1d.xaxis.title = "N_{ch}|_{" + make_estimator_title(est_dir.GetName()) + "}"
        ratio1d.yaxis.title = ytitle
        ratio1d.color = colors.next()
        if scale:
            ratio1d.Scale(scale)
        stack.append(ratio1d)
    c = plot_list_of_plottables(stack)
    c.name = "_".join(pids1) + "_div_" + "_".join(pids2)
    c.Write()


def _plot_particle_ratios_vs_refmult(f, sums, results_post, pids1, pids2, scale=None, ytitle=''):
    """
    plot and write to file the ratio of the two pid-lists (pids1/pids2). Plot is vs refmult.
    This function depends on the correlation histograms to be present in f
    """
    ratios = []
    colors = get_color_generator(10)
    refest = "EtaLt05"
    if not ytitle:
        ytitle = ", ".join(pids1) + " / " + ", ".join(pids2)

    for est_dir in sums:
        if est_dir.GetName() not in ['EtaLt05', 'EtaLt08', 'EtaLt15', 'Eta08_15', 'V0M']:
            continue
        h3d = asrootpy(est_dir.FindObject("fNch_pT_pid"))
        corr_hist = asrootpy(results_post.correlations.Get("corr_hist_{}_vs_{}".format(refest, est_dir.GetName())))

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
        ratio.xaxis.title = "N_{ch}|_{" + refest + "}"
        ratio.yaxis.title = ytitle
        ratio.title = make_estimator_title(est_dir.GetName())

        ratio.SetColor(colors.next())
        ratios.append(ratio)
    c = plot_list_of_plottables(ratios)
    c.name = "_".join(pids1) + "_div_" + "_".join(pids2)
    try:
        os.makedirs("./figures/pid_ratios_vs_ref_mult/")
    except OSError:
        pass
    # c.SaveAs("./figures/pid_ratios_vs_ref_mult/"+c.name+".pdf")
    c.Write()


def _make_dNdeta_and_event_counter(f, sums):
    # Loop over all estimators in the Sums list:
    log.info("Creating event counter and dN/deta plots")
    for est_dir in sums:
        if est_dir.GetName() == "Total":
            continue
        h3d = asrootpy(est_dir.FindObject('fNch_pT_pid'))
        h2d = asrootpy(est_dir.FindObject('feta_Nch'))
        nt = asrootpy(est_dir.FindObject("fEventTuple"))

        f.cd("results_post/" + est_dir.GetName())

        # Write out event counters vs. multiplicity
        h_event_counter = Hist1D(400, 0, 400,
                                 name=("event_counter"),
                                 title="Event counter vs. multiplicity in est region")
        nt.Project(h_event_counter.name, "nch", "ev_weight")
        h_event_counter.write()

        esti_title = "({0})".format(h3d.GetTitle()[31:])

        mean_nch = est_dir.FindObject("feta_Nch").GetMean(2)  # mean of yaxis
        # bin in standard step size up to max_nch; from there ibs all in one bin:
        max_nch = mean_nch * mean_mult_cutoff_factor
        l = get_dNdeta_binned_in_mult(h2d, h_event_counter, nch_max=max_nch, with_mb=True)
        c = plot_list_of_plottables(l, title="dN/d#eta vs. #eta " + esti_title)
        c.name = "dNdeta_summary"
        path = "./figures/results_post/{}/".format(est_dir.GetName())
        try:
            os.makedirs(path)
        except OSError:
            pass
        # c.SaveAs(path+c.name+".pdf")
        c.write()


def _make_hists_vs_pt(f, sums, results_post):
    ###########################################################
    # Category 2 on TWiki
    # create particle ratio vs pT plots

    log.info("Computing histograms vs pt")

    # Loop over all estimators in the Sums list:
    for est_dir in sums:
        # and do everything for weighted and unweighted:
        for postfix in postfixes:
            dirname = 'results_post/{}/pid_ratios'.format(est_dir.GetName() + postfix)
            try:
                f.mkdir(dirname, recurse=True)
            except:
                pass
            f.cd(dirname)
            h3d_orig = asrootpy(est_dir.FindObject('fNch_pT_pid' + postfix))
            h3d = asrootpy(h3d_orig.RebinX(rebin_mult, h3d_orig.name + "rebinned"))
            mean_nch = est_dir.FindObject("feta_Nch").GetMean(2)  # mean of yaxis
            # bin in standard step size up to max_nch; from there ibs all in one bin:
            max_nch = mean_nch * mean_mult_cutoff_factor
            esti_title = "({0})".format(h3d.title[31:])

            mult_pt_dir = results_post.FindObject(est_dir.GetName()).Get("mult_pt")

            hs = create_stack_pid_ratio_over_pt(mult_pt_dir, [kANTIPROTON, kPROTON], [kPIMINUS, kPIPLUS], max_nch)
            hs.title = "p/#pi^{+-} vs. p_{T} " + "{}".format(esti_title)
            c = plot_list_of_plottables(hs)
            c.name = "proton_over_pich__vs__pt"
            c.write()

            hs = create_stack_pid_ratio_over_pt(mult_pt_dir, [kANTIXI, kXI], [kPIMINUS, kPIPLUS], max_nch)
            hs.title = "#Xi/#pi^{+-} vs. p_{T} " + "{}".format(esti_title)
            c = plot_list_of_plottables(hs)
            c.name = "Xi_over_pich__vs__pt"
            c.write()

            hs = create_stack_pid_ratio_over_pt(mult_pt_dir, [kOMEGAMINUS, kOMEGAPLUS], [kPIMINUS, kPIPLUS], max_nch)
            hs.title = "#Omega_{ch}/#pi^{+-} vs. p_{T} " + "{}".format(esti_title)
            c = plot_list_of_plottables(hs)
            c.name = "OmegaCh_over_pich__vs__pt"
            c.write()

            # Ratios to pi0
            hs = create_stack_pid_ratio_over_pt(mult_pt_dir, [kPIMINUS, kPIPLUS], [kPI0], max_nch)
            hs.title = "#pi^{+-}/#pi^{0} vs. p_{T} " + "{}".format(esti_title)
            c = plot_list_of_plottables(hs)
            c.name = "pich_over_pi0__vs__pt"
            c.write()

            hs = create_stack_pid_ratio_over_pt(mult_pt_dir, [kANTIPROTON, kPROTON], [kPI0], max_nch)
            hs.title = "p/#pi^{0} vs. p_{T} " + "{}".format(est_dir.GetName())
            c = plot_list_of_plottables(hs)
            c.name = "proton_over_pi0__vs__pt"

            c.write()

            hs = create_stack_pid_ratio_over_pt(mult_pt_dir, [kK0S], [kPI0], max_nch)
            hs.title = "K0S/#pi^{0} vs. p_{T} " + "{}".format(est_dir.GetName())
            c = plot_list_of_plottables(hs)
            c.name = "K0S_over_pi0__vs__pt"
            c.write()

            hs = create_stack_pid_ratio_over_pt(mult_pt_dir, [kANTILAMBDA, kLAMBDA], [kPI0], max_nch)
            hs.title = "#Lambda/#pi^{0} vs. p_{T} " + "{}".format(est_dir.GetName())
            c = plot_list_of_plottables(hs)
            c.name = "Lambda_over_pi0__vs__pt"
            c.write()

            hs = create_stack_pid_ratio_over_pt(mult_pt_dir, [kANTIXI, kXI], [kPI0], max_nch)
            hs.title = "#Xi/#pi^{0} vs. p_{T} " + "{}".format(est_dir.GetName())
            c = plot_list_of_plottables(hs)
            c.name = "Xi_over_pi0__vs__pt"
            c.write()

            hs = create_stack_pid_ratio_over_pt(mult_pt_dir, [kOMEGAMINUS, kOMEGAPLUS], [kPI0], max_nch)
            hs.title = "#Omega_{ch}/#pi^{0} vs. p_{T} " + "{}".format(est_dir.GetName())
            c = plot_list_of_plottables(hs)
            c.name = "OmegaCh_over_pi0__vs__pt"
            c.write()

            # Ratios to K0S
            hs = create_stack_pid_ratio_over_pt(mult_pt_dir, [kANTIPROTON, kPROTON], [kK0S], max_nch)
            hs.title = "p/K^{0}_{S} vs. p_{T} " + "{}".format(esti_title)
            c = plot_list_of_plottables(hs)
            c.name = "proton_over_K0S__vs__pt"
            c.write()

            hs = create_stack_pid_ratio_over_pt(mult_pt_dir, [kANTILAMBDA, kLAMBDA], [kK0S], max_nch)
            hs.title = "#Lambda/K^{0}_{S} vs. p_{T} " + "{}".format(esti_title)
            c = plot_list_of_plottables(hs)
            c.name = "Lambda_over_K0S__vs__pt"
            c.write()

            hs = create_stack_pid_ratio_over_pt(mult_pt_dir, [kANTIXI, kXI], [kK0S], max_nch)
            hs.title = "#Xi/K^{0}_{S} vs. p_{T} " + "{}".format(esti_title)
            c = plot_list_of_plottables(hs)
            c.name = "Xi_over_K0S__vs__pt"
            c.write()

            hs = create_stack_pid_ratio_over_pt(mult_pt_dir, [kOMEGAMINUS, kOMEGAPLUS], [kK0S], max_nch)
            hs.title = "#Omega_{ch}/K^{0}_{S} vs. p_{T} " + "{}".format(esti_title)
            c = plot_list_of_plottables(hs)
            c.name = "OmegaCh_over_K0S__vs__pt"
            c.write()


def _make_PNch_plots(f, sums, results_post):
    log.info("Creating P(Nch) summary plot")
    hists_PNch_vs_estmult = []
    for est_name in considered_ests:
        h_tmp = get_PNch_vs_estmult(results_post, est_name)
        h_tmp.title = make_estimator_title(est_name)
        if h_tmp.Integral() > 0:
            hists_PNch_vs_estmult.append(h_tmp)

    [h.Scale(1.0 / h.Integral()) for h in hists_PNch_vs_estmult]

    c = plot_list_of_plottables(hists_PNch_vs_estmult, logy=True)
    c.name = "PNch_summary"
    f.cd("results_post")
    c.write()

    log.info("Creating P(Nch_est) and P(Nch_refest) histograms")
    mult_bin_size = 10
    for ref_est in ref_ests:
        for est in considered_ests:
            corr_hist = get_NchEst1_vs_NchEst2(results_post, ref_est, est)
            nch_max = corr_hist.xaxis.GetNbins()
            mean_nch_est = corr_hist.GetMean(1)  # mean of x axis
            nch_cutoff = mean_nch_est * mean_mult_cutoff_factor
            hists_PNch_vs_estmult_binned_in_refmult = []
            hists_PNch_vs_refmult_binned_in_estmult = []
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
                h_vs_est.title = "{} < N_{{ch}}^{{{}}} < {}".\
                                 format(nch_lower_edge, ref_est, nch_upper_edge)
                h_vs_est.xaxis.title = "N_{{ch}}^{{{}}}".format(est)
                if h_vs_est.Integral() > 0:
                    hists_PNch_vs_estmult_binned_in_refmult.append(h_vs_est)

                # vs ref_mult:
                corr_hist.yaxis.SetRange(0, 0)  # reset y axis
                corr_hist.xaxis.SetRange(nch_lower_edge, nch_upper_edge)
                h_vs_ref = asrootpy(corr_hist.ProjectionY(gen_random_name()))
                h_vs_ref.title = "{} < N_{{ch}}^{{{}}} < {}".\
                                 format(nch_lower_edge, est, nch_upper_edge)
                h_vs_ref.xaxis.title = "N_{{ch}}^{{{}}}".format(ref_est)
                if h_vs_ref.Integral() > 0:
                    hists_PNch_vs_refmult_binned_in_estmult.append(h_vs_ref)

                if is_last_bin:
                    break

            # Normalize each slice:
            [h.Scale(1.0 / h.Integral()) for h in hists_PNch_vs_estmult_binned_in_refmult]
            [h.Scale(1.0 / h.Integral()) for h in hists_PNch_vs_refmult_binned_in_estmult]

            # Plot:
            f.cd("results_post/" + est)

            # vs est_mult
            c = plot_list_of_plottables(hists_PNch_vs_estmult_binned_in_refmult, logy=True)
            c.name = "PNch{}_binned_in_Nch{}".format(est, ref_est)
            c.write()

            # vs est_mult
            c = plot_list_of_plottables(hists_PNch_vs_refmult_binned_in_estmult, logy=True)
            c.name = "PNch{}_binned_in_Nch{}".format(ref_est, est)
            c.write()


def _make_mult_vs_pt_plots(f, sums, results_post):
    log.info("Makeing 2D mult pt plots for each particle kind")
    for est_dir in (somedir for somedir in sums if somedir.GetName() in considered_ests):
        dir_name = "results_post/" + est_dir.GetName() + "/mult_pt"
        try:
            f.mkdir(dir_name, recurse=True)
        except ValueError:
            pass
        f.cd(dir_name)

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
        res_dir_str = "results_post/" + est_dir.GetName()
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
        title = 'dN/d#eta|_{mult} / dN/d#eta|_{MB} '
        c = plot_list_of_plottables(ratios, title)
        c.name = "dNdeta_ratio_to_mb_canvas"
        f.cd(res_dir_str)
        c.Write()


def _make_correlation_plots(f, sums, results_post):
    # Make correlations between estimators
    log.info("Correlating N_ch of each estimator")
    corr_dir = 'results_post/correlations'
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
    ratio_vs_ref_dir = 'results_post/pid_ratios_vs_refmult'
    ratio_vs_est_dir = 'results_post/pid_ratios_vs_estmult'
    try:
        f.mkdir(ratio_vs_ref_dir, recurse=True)
    except:
        pass
    f.cd(ratio_vs_ref_dir)

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

    try:
        f.mkdir(ratio_vs_est_dir, recurse=True)
    except:
        pass
    f.cd(ratio_vs_est_dir)
    _plot_particle_ratios_vs_estmult(f, sums, results_post, ['321', '-321'], ['310'],
                                     scale=.5, ytitle="(K^{+} + K^{-}) / (2*K_{S}^{0})")


def _delete_results_dir(f, sums):
    try:
        # delete old result directory
        f.rmdir('results_post')
    except:
        pass
    f.mkdir('results_post')


def _mk_results_dir(f, sums):
    for est_dir in sums:
        try:
            resdir = f.results_post.mkdir(est_dir.GetName())
            resdir.Write()
        except:
            pass

if __name__ == "__main__":
    # go into batch mode
    ROOT.gROOT.SetBatch(False)

    log = log["/post"]  # set name of this script in logger
    log.info("IsBatch: {0}".format(ROOT.gROOT.IsBatch()))  # Results in "DEBUG:myapp] Hello"

    # Rebin multiplicity with factor:
    rebin_mult = 10
    ref_ests = ['EtaLt05', ]
    considered_ests = ['EtaLt05', 'EtaLt08', 'EtaLt15', 'Eta08_15', 'V0M', 'V0A', 'V0C']
    postfixes = ["", ]  # "_unweighted"]

    functions = [
        _delete_results_dir,
        _mk_results_dir,
    ]

    with root_open(sys.argv[1], 'update') as f:
        sums = f.Sums
        for func in functions:
            func(f, sums)
        results_post = f.results_post
        _make_dNdeta_and_event_counter(f, sums)
        _make_correlation_plots(f, sums, results_post)
        _make_PNch_plots(f, sums, results_post)

        _make_mult_vs_pt_plots(f, sums, results_post)
    with root_open(sys.argv[1], 'update') as f:
        sums = f.Sums
        results_post = f.results_post
        _make_hists_vs_pt(f, sums, results_post)  # needs updated results_post!
        _make_dNdeta_mb_ratio_plots(f, sums, results_post)
        _make_pid_ratio_plots(f, sums, results_post)

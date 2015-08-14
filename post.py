"""Run post analysis on the given file"""
import sys, os, string, random
import ipdb

if len(sys.argv) != 2:
    print "Usage: python ./post.py path_to_root_file.root"
    quit()

from rootpy.io import root_open
from rootpy import asrootpy, ROOT, log
from rootpy.plotting import HistStack, Hist1D, Hist2D, Canvas
from post_utils import create_dNdeta_stack,\
    plot_histogram_stack, create_stack_pid_ratio_over_pt,\
    create_hist_pid_ratio_over_mult,\
    create_canonnical_avg_from_stacks,\
    divide_stacks, create_graph_pided_refest_vs_pidcount,\
    plot_list_of_plottables, remove_zero_value_points, remove_non_mutual_points,\
    remove_points_with_equal_x, remove_points_with_x_err_gt_1NchRef

def gen_random_name():
    """Generate a random name for temp hists"""
    return ''.join(random.choice(string.ascii_letters + string.digits) for _ in range(25))

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


def _plot_particle_ratios_vs_estmult(f, sums, results_post, pids1, pids2, scale=None, title=''):
    stack = list()
    nesti = len(sums)
    for i, est_dir in enumerate(sums):
        h3d = asrootpy(est_dir.FindObject("fNch_pT_pid"))
        pids1hists = []
        pids2hists = []
        for pid in pids1:
            h3d.zaxis.SetRange(h3d.zaxis.FindBin(pid), h3d.zaxis.FindBin(pid))
            h = asrootpy(h3d.Project3D("yx"))
            h.SetName(gen_random_name())
            pids1hists.append(h)

        for pid in pids2:
            h3d.zaxis.SetRange(h3d.zaxis.FindBin(pid), h3d.zaxis.FindBin(pid))
            h = asrootpy(h3d.Project3D("yx"))
            h.SetName(gen_random_name())
            pids2hists.append(h)

        # sum up each histogram
        pids1_px = asrootpy(sum(pids1hists).ProjectionX())
        pids2_px = asrootpy(sum(pids2hists).ProjectionX())
        ratio1d = pids1_px / pids2_px
        ratio1d.SetTitle(est_dir.GetName())
        ratio1d.color = int(800 + i*100.0/nesti) +1
        if scale:
            ratio1d.Scale(scale)
        stack.append(ratio1d)

    if not title:
        title = "{} div {}".format(str(pids1), str(pids2))
    c = plot_list_of_plottables(stack, title)
    c.name = "_".join(pids1) + "_div_" + "_".join(pids2)

    c.Write()


def _plot_particle_ratios_vs_refmult(f, sums, results_post, pids1, pids2, scale=None, ytitle=''):
    """
    plot and write to file the ratio of the two pid-lists (pids1/pids2). Plot is vs refmult.
    This function depends on the correlation histograms to be present in f
    """
    ratios = []
    for i, est_dir in enumerate(sums):
        h3d = asrootpy(est_dir.FindObject("fNch_pT_pid"))
        chist = asrootpy(results_post.correlations.Get("corr_hist_EtaLt05_vs_{}".format(est_dir.GetName())))
        try:
            g1 = create_graph_pided_refest_vs_pidcount(h3d, chist, pids1)
            g2 = create_graph_pided_refest_vs_pidcount(h3d, chist, pids2)
        except IndexError:
            log.info("Could not process pid ratio for {} for estimator {}".format(
                "_".join(pids1) + "_div_" + "_".join(pids2),
                est_dir.GetName()))
            continue
        remove_zero_value_points(g1)
        remove_zero_value_points(g2)
        remove_points_with_x_err_gt_1NchRef(g1)
        remove_points_with_x_err_gt_1NchRef(g2)
        remove_points_with_equal_x(g1)
        remove_points_with_equal_x(g2)
        remove_non_mutual_points(g1, g2)
        try:
            ratio = g1 / g2
        except ZeroDivisionError:
            print "ZeroDivisionError in {}".format(est_dir.GetName())
            continue
        if not ytitle:
            title = g1.title + " / " + g2.title
        else:
            title = ''
        if scale:
            ratio.Scale(scale)
        ratio.yaxis.title = ytitle
        ratio.title = est_dir.GetName()
        ratio.xaxis.title = g1.xaxis.title
        ratio.SetColor(800 + int(100.0/3)*(i) + 1)
        ratios.append(ratio)
    c = plot_list_of_plottables(ratios, title)
    c.name = "_".join(pids1) + "_div_" + "_".join(pids2)
    try:
        os.makedirs("./figures/pid_ratios_vs_ref_mult/")
    except OSError:
        pass
    #c.SaveAs("./figures/pid_ratios_vs_ref_mult/"+c.name+".pdf")
    c.Write()


def _make_dNdeta_and_event_counter(f, sums, results_post):
    # Loop over all estimators in the Sums list:
    for est_dir in sums:
        if est_dir.GetName() == "Total":
            continue
        log.info("Computing histograms for {0}".format(est_dir.GetName()))
        # and do everything for weighted and unweighted:
        for postfix in postfixes:
            h3d = asrootpy(est_dir.FindObject('fNch_pT_pid' + postfix))
            h2d = asrootpy(est_dir.FindObject('feta_Nch' + postfix))
            nt = asrootpy(est_dir.FindObject("fEventTuple"))

            f.cd("results_post/" + est_dir.GetName() + postfix)

            # Write out event counters vs. multiplicity
            h_event_counter = Hist1D(400, 0, 400,
                                     name=("event_counter" + postfix),
                                     title="Event counter vs. multiplicity in est region")
            if postfix:  # empty string is falsy
                nt.Project(h_event_counter.name, "nch")
            else:
                nt.Project(h_event_counter.name, "nch", "ev_weight")
            h_event_counter.write()

            esti_title = "({0})".format(h3d.GetTitle()[31:])
            ###########################################################
            # Category 1 on TWiki
            # create dN/deta stack for the current estimator
            hs = create_dNdeta_stack(h2d, h_event_counter)
            hs.SetTitle("dN/d#eta vs. #eta " + esti_title)
            c = plot_histogram_stack(hs)
            c.name = "dNdeta_summary"
            path = "./figures/results_post/{}/".format(est_dir.GetName())
            try:
                os.makedirs(path)
            except OSError:
                pass
            #c.SaveAs(path+c.name+".pdf")
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
            dirname = 'results_post/{}/pid_ratios'.format(est_dir.GetName()+postfix)
            try:
                f.mkdir(dirname, recurse=True)
            except:
                pass
            f.cd(dirname)            
            h3d_orig = asrootpy(est_dir.FindObject('fNch_pT_pid' + postfix))
            h3d = asrootpy(h3d_orig.RebinX(rebin_mult, h3d_orig.name+"rebinned"))
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
    # Create P(Nch) plots
    log.info("Creating P(Nch) plots")
    for postfix in postfixes:
        pNch_summary = HistStack()  # One HistStack!
        pNch_summary.title = "P(N_{ch}^{est}) summary "

        step_size = 10  # binning in Nch^est
        # make ntuples:
        nt0 = sums[0].FindObject("fEventTuple")
        nt0.SetAlias(sums[0].GetName(), "fEventTuple")
        #previous_estimator_names = [sums[0].GetName(), ]
        for est_dir in sums[1:]:
            nt0.AddFriend(est_dir.FindObject("fEventTuple"), est_dir.GetName())
        for est_dir in sums:
            ############################################################
            # Summary plot P(Nch^{est}):
            name = "PNch_"+est_dir.GetName()+ postfix
            axis_title = ";N_{ch}^{est}"
            h_summary = Hist1D(400, 0, 400, name=name,
                               title=est_dir.GetName() + axis_title)
            nt0.Project(h_summary.name,
                        "{}.nch".format(est_dir.GetName()),
                        "{}.ev_weight".format(est_dir.GetName()))
            h_summary.Scale(1.0/h_summary.Integral())
            pNch_summary.Add(h_summary)

            #############################################################
            # Plot binned in mult_ref vs. mutl_est and vice versa
            for ref_est in ref_ests:
                axis_title_vs_ref = ";N_{{ch}}^{{{}}}".format(ref_est)
                axis_title_vs_est = ";N_{{ch}}^{{{}}}".format(est_dir.GetName())
                pNch_ref_binned_Nch_est = HistStack()
                pNch_est_binned_Nch_ref = HistStack()
                pNch_ref_binned_Nch_est.title = "P(N_{{ch}}^{{{}}}) binned in N_{{ch}}^{{{}}}"\
                                                .format(ref_est, est_dir.GetName())
                pNch_est_binned_Nch_ref.title = "P(N_{{ch}}^{{{}}}) binned in N_{{ch}}^{{{}}}"\
                                       .format(est_dir.GetName(), ref_est)

                for nch_interv_min in range(0, 401, step_size):
                    nch_interv_max = nch_interv_min + step_size
                    name = "PNch_{}_lt_mult_ref_lt_{}".format(nch_interv_min, nch_interv_max)
                    h_tmp_vs_ref = Hist1D(400, 0, 400, name=name)
                    name = "PNch_{}_lt_mult_est_lt_{}".format(nch_interv_min, nch_interv_max)
                    h_tmp_vs_est = Hist1D(400, 0, 400, name=name)
                    h_tmp_vs_ref.title = "{} < N_{{ch}}^{{{}}} < {}".\
                                         format(nch_interv_min, est_dir.GetName(), nch_interv_max) + axis_title_vs_ref
                    h_tmp_vs_est.title = "{} < N_{{ch}}^{{{}}} < {}".\
                                         format(nch_interv_min, ref_est, nch_interv_max) + axis_title_vs_est
                    nt0.Project(h_tmp_vs_ref.name,
                                "{}.nch".format(ref_est),
                                "{2}.ev_weight*({0} < {2}.nch && {2}.nch < {1})"\
                                .format(nch_interv_min, nch_interv_max, est_dir.GetName()))
                    nt0.Project(h_tmp_vs_est.name,
                                "{}.nch".format(est_dir.GetName()),
                                "{2}.ev_weight*({0} < {2}.nch && {2}.nch < {1})"\
                                .format(nch_interv_min, nch_interv_max, ref_est))
                    try:
                        h_tmp_vs_ref.Scale(1.0/h_tmp_vs_ref.Integral())
                    except ZeroDivisionError:
                        pass
                    else:
                        # only add to stack if there are values in this Nch interval
                        pNch_ref_binned_Nch_est.Add(h_tmp_vs_ref)
                    try:
                        h_tmp_vs_est.Scale(1.0/h_tmp_vs_est.Integral())
                    except ZeroDivisionError:
                        pass
                    else:
                        # only add to stack if there are values in this Nch interval
                        pNch_est_binned_Nch_ref.Add(h_tmp_vs_est)
                f.cd("results_post/" + est_dir.GetName() + postfix)
                # plot vs ref mult 
                c = plot_histogram_stack(pNch_ref_binned_Nch_est)
                c.name = "PNch{}_binned_in_Nch{}".format(ref_est, est_dir.GetName())
                c.FindObject("plot").SetLogy(1)
                c.write()

                # plot vs est mult 
                c = plot_histogram_stack(pNch_est_binned_Nch_ref)
                c.name = "PNch{}_binned_in_Nch{}".format(est_dir.GetName(), ref_est)
                c.FindObject("plot").SetLogy(1)
                c.write()


        # Write summary to disk
        c = plot_histogram_stack(pNch_summary)
        c.FindObject("plot").SetLogy(1)
        c.name = "PNch_summary" + postfix
        f.cd("results_post")
        c.write()


def _make_mult_vs_pt_plots(f, sums, results_post):
    log.info("Makeing 2D mult pt plots for each particle kind")
    for est_dir in sums:
        if est_dir.GetName() == "Total":
            continue
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


def _create_ratio_to_mb_stack(stack):
    """Given a stack of histograms, create the ratio of each multiplicity bin to the MB (all binns combined)"""
    mb_name = 'mb_dNdeta'
    mb_hist = next(h for h in stack.GetHists() if h.name==mb_name)
    ratio_stack = HistStack()
    for h in stack:
        if h.name == mb_name:
            break
        h_tmp = h / mb_hist
        ratio_stack.Add(h_tmp)
    return ratio_stack

        
def _make_PNch_ratio_plots(f, sums, results_post):
    # Create ratio plots; depends on the previously created histograms
    log.info("Creating ratios of dN/deta plots for each multiplicity bin")
    for postfix in postfixes:
        for est_dir in sums:
            if est_dir.GetName() == "Total":
                continue
            res_dir_str = "results_post/" + est_dir.GetName() + postfix
            stack = asrootpy(results_post\
                             .Get(est_dir.GetName() + postfix)\
                             .Get("dNdeta_summary")\
                             .FindObject('dNdeta_stack'))
            ratio_stack = _create_ratio_to_mb_stack(stack)
            ratio_stack.title = 'dN/d#eta|_{mult} / dN/d#eta|_{MB} '
            #ratio_stack.name = "dNdeta_ratio_to_mb"
            c = plot_histogram_stack(ratio_stack)
            c.name = "dNdeta_ratio_to_mb_canvas"
            f.cd(res_dir_str)
            c.Write()
            
        
        # # Get the  dN/deta stack for each estimator in order to calc 
        # # the cannonical average for all _other_ estimators later on
        # # P(N_ch): one stack containing all estimators!
        # dNdeta_stacks = []        # List of stacks!
        # #pNch_stack = HistStack()  # One HistStack!
        # #pNch_stack.name = "P(N_{ch}^{est}) summary "
        # for est_dir in sums:
        #     
        #     dNdeta_stacks.append(asrootpy(f.results_post\
        #                            .Get(est_dir.GetName() + postfix)\
        #                            .Get("dNdeta_summary")\
        #                            .FindObject('dNdeta_stack')))

        # # looping over file again in order to have the estimator name handy,
        # for i, est_dir in enumerate(sums):
        #     res_dir_str = "results_post/" + est_dir.GetName() + postfix + "/ratios_to_other_est"
        #     try:
        #         f.mkdir(res_dir_str, recurse=True)
        #     except:
        #         pass
        #     # calc cannonical avg withouth the current hist
        #     other_dNdeta_stacks = dNdeta_stacks[:i] + dNdeta_stacks[i+1:]
        #     avg_stack = create_canonnical_avg_from_stacks(other_dNdeta_stacks)
        #     ratio = divide_stacks(dNdeta_stacks[i], avg_stack)
        #     ratio.title = ('Ratio of '
        #                    + dNdeta_stacks[i].title
        #                    + ' to cannonical average')

        #     c = plot_histogram_stack(ratio)
        #     c.name = dNdeta_stacks[i].name + '_ratio_cannonical_avg'
        #     c.Update()
        #     f.cd(res_dir_str)
        #     c.Write()

        #     # create ratio between P(N_ch) and cannonical_avg
        #     # avg = pNch_stack[0].Clone()
        #     # avg.Clear()
        #     # other_PNch = pNch_stack[:i] + pNch_stack[i+1:]
        #     # [avg.Add(h) for h in other_PNch]
        #     # avg.Scale(1.0/(len(other_PNch)))

        #     # ratio = pNch_stack[i] / avg
        #     # ratio.name = "{0}_div_by_can_avg".format(pNch_stack[i].name)
        #     # ratio.title = "Ratio of {} over cannonical avg".format(pNch_stack[i].title)
        #     # f.cd(res_dir_str)
        #     # ratio.Write()

def _make_correlation_plots(f, sums, results_post):
    # Make correlations between estimators
    log.info("Correlating N_ch of each estimator")
    corr_dir = 'results_post/correlations'
    try:
        f.mkdir(corr_dir, recurse=True)
    except:
        pass
    # only for weighted case
    # Take ntuple from the first estimator and then add friends to this one
    nt0 = sums[0].FindObject("fEventTuple")
    nt0.SetAlias(sums[0].GetName(), "fEventTuple")
    #previous_estimator_names = [sums[0].GetName(), ]
    # build ntuple
    for est_dir in sums[1:]:
        nt0.AddFriend(est_dir.FindObject("fEventTuple"), est_dir.GetName())
    for ref_est in ref_ests:
        for est_dir in sums:
            corr_hist = Hist2D(400, 0, 400,
                               400, 0, 400,
                               name="corr_hist_{}_vs_{}".format(ref_est, est_dir.GetName()))
            # Lables are deliberatly swaped, see Projection below!
            corr_hist.title = ("Correlation N_{{ch}} in {0} and {1};N_{{ch}} {1};N_{{ch}} {0}"\
                               .format(ref_est, est_dir.GetName()))

            # this projects onto y:x, to make coding more adventurous 
            nt0.Project(corr_hist.name, "{0}.nch:{1}.nch".format(ref_est, est_dir.GetName()),
                        "ev_weight")
            corr_hist.drawstyle = 'colz'
            f.cd(corr_dir)
            corr_hist.write()

def _make_pid_ratio_plots(f, sums, results_post):
    log.info("Correlating identified particle ratios")
    ratio_vs_ref_dir = 'results_post/pid_ratios_vs_refmult'
    ratio_vs_est_dir = 'results_post/pid_ratios_vs_estmult'
    try:
        f.mkdir(ratio_vs_ref_dir, recurse=True)
    except:
        pass
    f.cd(ratio_vs_ref_dir)

    # Proton / pi_ch
    _plot_particle_ratios_vs_refmult(f, sums, results_post, ['-2212','2212'], ['-211', '211'],
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
    _plot_particle_ratios_vs_refmult(f, sums, results_post, ['-211', '211'], ['111',],
                                     ytitle="#pi^{+-}/#pi^{0}")
    # proton / pi0
    _plot_particle_ratios_vs_refmult(f, sums, results_post, ['-2212', '2212'], ['111',],
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
                                     scale=.5, title="(K^{+} + K^{-}) / (2*K_{S}^{0})")


def _reset_results_dir():
    with root_open(sys.argv[1], 'update') as f:
        try:
        # delete old result directory
            f.rmdir('results_post')
        except:
            pass
        f.mkdir('results_post')

if __name__ == "__main__":
    # go into batch mode
    ROOT.gROOT.SetBatch(True)

    log = log["/post"]  # set name of this script in logger
    log.info("IsBatch: {0}".format(ROOT.gROOT.IsBatch())) # Results in "DEBUG:myapp] Hello"

    # Rebin multiplicity with factor:
    rebin_mult = 10
    ref_ests = ['EtaLt05',]
    postfixes = ["",]  # "_unweighted"]

    _reset_results_dir()
    with root_open(sys.argv[1], 'update') as f:
        sums = f.Sums
        for est_dir in sums:
            try:
                resdir = f.results_post.mkdir(est_dir.GetName())
                resdir.Write()
            except:
                pass
        results_post = f.results_post
        _make_dNdeta_and_event_counter(f, sums, results_post)
        _make_mult_vs_pt_plots(f, sums, results_post)
        _make_PNch_ratio_plots(f, sums, results_post)
        _make_hists_vs_pt(f, sums, results_post)
        _make_PNch_plots(f, sums, results_post)
        _make_correlation_plots(f, sums, results_post)
        _make_pid_ratio_plots(f, sums, results_post)



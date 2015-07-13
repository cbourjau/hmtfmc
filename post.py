"""Run post analysis on the given file"""
import sys
if len(sys.argv) != 2:
    print "Usage: python ./post.py path_to_root_file.root"
    quit()

from rootpy.io import root_open
from rootpy import asrootpy, ROOT, log
from rootpy.plotting import HistStack, Hist1D, Hist2D
from post_utils import create_dNdeta_stack,\
    plot_histogram_stack, create_stack_pid_ratio_over_pt,\
    create_hist_pid_ratio_over_mult,\
    create_canonnical_avg_from_stacks,\
    divide_stacks

# go into batch mode
ROOT.gROOT.SetBatch(True)

log = log["/post"]  # set name of this script in logger
log.info("IsBatch: {0}".format(ROOT.gROOT.IsBatch())) # Results in "DEBUG:myapp] Hello"

# Rebin multiplicity with factor:
rebin_mult = 10

with root_open(sys.argv[1], 'update') as f_post:
    try:
        # delete old result directory
        f_post.rmdir('results_post')
    except:
        pass
    # Loop over all estimators in the Sums list:
    for est_dir in f_post.Sums:
        log.info("Computing histograms for {0}".format(est_dir.GetName()))
        # and do everything for weighted and unweighted:
        for postfix in ["", "_unweighted"]:
            h3d = est_dir.FindObject('festi_pT_pid' + postfix)
            h3d = asrootpy(h3d)
            h3d.RebinX(rebin_mult)
            #h3d.Scale(1.0/rebin_mult)
            
            h2d = est_dir.FindObject('fdNdeta' + postfix)
            h2d = asrootpy(h2d)
            h2d.RebinY(rebin_mult)

            nt = asrootpy(est_dir.FindObject("fevent_counter"))
            
            res_dir = f_post.mkdir("results_post/" + est_dir.GetName() + postfix, recurse=True)
            res_dir.write()

            # res_dir = f_post.results_post.FindObject(est_dir.GetName())
            f_post.results_post.cd(est_dir.GetName() + postfix)

            # Write out event counters vs. multiplicity
            h_event_counter = Hist1D(100, 0, 100,
                                     name=("event_counter" + postfix),
                                     title="Event counter vs. multiplicity in est region")
            if postfix:  # empty string is falsy
                nt.Project(h_event_counter.name, "nch")
            else:
                nt.Project(h_event_counter.name, "nch", "ev_weight")
            h_event_counter.write()

            # rescale for later operations
            h_event_counter.Rebin(rebin_mult)

            esti_title = "({0})".format(h3d.title[38:])
            ###########################################################
            # Category 1 on TWiki
            # create dN/deta stack for the current estimator
            hs = create_dNdeta_stack(h2d, h_event_counter)
            hs.SetTitle("dN/d#eta vs. #eta " + esti_title)
            c = plot_histogram_stack(hs)
            c.name = "dNdeta_summary"
            c.write()
            quit()

            # P(N_ch)
            h_PN_ch = asrootpy(h_event_counter.clone("PN_ch"))
            h_PN_ch.Scale(1.0/h_PN_ch.Integral())
            h_PN_ch.title = "P(N_{ch})" + esti_title
            h_PN_ch.yaxis.title = "P(N_{ch})"
            h_PN_ch.write()

            ###########################################################
            # Category 2 on TWiki
            # create particle ratio vs pT plots
            # Ratios to pich
            # hs = create_stack_pid_ratio_over_pt(h3d, [0], [5,6])
            # hs.title = "p/#pi^{+-} vs. p_{T} " + "{}".format(esti_title)
            # c = plot_histogram_stack(hs)
            # c.name = "proton_over_pich__vs__pt"
            # c.write()

            # hs = create_stack_pid_ratio_over_pt(h3d, [8], [5,6])
            # hs.title= "#Xi/#pi^{+-} vs. p_{T} " + "{}".format(esti_title)
            # c = plot_histogram_stack(hs)
            # c.name = "Xi_over_pich__vs__pt"
            # c.write()

            # hs = create_stack_pid_ratio_over_pt(h3d, [9,10], [5,6])
            # hs.title= "\Omega_{ch}/\pi^{+-} vs. p_{T} " + "{}".format(esti_title)
            # c = plot_histogram_stack(hs)
            # c.name = "OmegaCh_over_pich__vs__pt"
            # c.write()

            # # Ratios to pi0
            # # hs = create_stack_pid_ratio_over_pt(h3d, [5,6], [7])
            # # hs.title = "p^{+-}/\pi^{0} vs. p_{T} " + "{}".format(esti_title)
            # # c = plot_histogram_stack(hs)
            # # c.name = "pich_over_pi0__vs__pt"
            # # c.write()

            # # c = plot_histogram_stack(create_stack_pid_ratio_over_pt(h3d, [0], [7]))
            # # c.name = "proton_over_pi0__vs__pt"
            # # c.title= "p/#pi^{0} vs. p_{T} " + "{}".format(est_dir.GetName())
            # # c.write()

            # # c = plot_histogram_stack(create_stack_pid_ratio_over_pt(h3d, [2], [7]))
            # # c.name = "K0S_over_pi0__vs__pt"
            # # c.title= "K0S/#pi^{0} vs. p_{T} " + "{}".format(est_dir.GetName())
            # # c.write()

            # # c = plot_histogram_stack(create_stack_pid_ratio_over_pt(h3d, [1], [7]))
            # # c.name = "Lambda_over_pi0__vs__pt"
            # # c.title= "#Lambda/#pi^{0} vs. p_{T} " + "{}".format(est_dir.GetName())
            # # c.write()

            # # c = plot_histogram_stack(create_stack_pid_ratio_over_pt(h3d, [8], [7]))
            # # c.name = "Xi_over_pi0__vs__pt"
            # # c.title= "#Xi/#pi^{0} vs. p_{T} " + "{}".format(est_dir.GetName())
            # # c.write()

            # # c = plot_histogram_stack(create_stack_pid_ratio_over_pt(h3d, [9,10], [7]))
            # # c.name = "OmegaCh_over_pi0__vs__pt"
            # # c.title= "#Omega_{ch}/#pi^{0} vs. p_{T} " + "{}".format(est_dir.GetName())
            # # c.write()

            # # Ratios to K0S
            # hs = create_stack_pid_ratio_over_pt(h3d, [0], [2])
            # hs.title= "p/K^{0}_{S} vs. p_{T} " + "{}".format(esti_title)
            # c = plot_histogram_stack(hs)
            # c.name = "proton_over_K0S__vs__pt"
            # c.write()

            # hs = create_stack_pid_ratio_over_pt(h3d, [1], [2])
            # hs.title= "#Lambda/K^{0}_{S} vs. p_{T} " + "{}".format(esti_title)
            # c = plot_histogram_stack(hs)
            # c.name = "Lambda_over_K0S__vs__pt"
            # c.write()

            # hs = create_stack_pid_ratio_over_pt(h3d, [8], [2])
            # hs.title= "#Xi/K^{0}_{S} vs. p_{T} " + "{}".format(esti_title)
            # c = plot_histogram_stack(hs)
            # c.name = "Xi_over_K0S__vs__pt"
            # c.write()

            # hs = create_stack_pid_ratio_over_pt(h3d, [9,10], [2])
            # hs.title= "#Omega_{ch}/K^{0}_{S} vs. p_{T} " + "{}".format(esti_title)
            # c = plot_histogram_stack(hs)
            # c.name = "OmegaCh_over_K0S__vs__pt"
            # c.write()
            # #############################################################33
            # #Category 3 on Twiki
            # h = create_hist_pid_ratio_over_mult(h3d, [0], [5,6])
            # h.SetTitle('p/#pi^{+-} vs. N_{ch}|_{est} ' + "{}".format(esti_title))
            # h.name = "proton_over_pich__vs__mult"
            # h.write()

            # h = create_hist_pid_ratio_over_mult(h3d, [2], [5,6])
            # h.SetTitle('K^{0}_{S}/#pi^{+-} vs. N_{ch}|_{est} ' + "{}".format(esti_title))
            # h.name = "K0S_over_pich__vs__mult"
            # h.write()

            # h = create_hist_pid_ratio_over_mult(h3d, [1], [5,6])
            # h.title= "#Lambda/#pi^{+-} vs. N_{ch}|_{est} " + "{}".format(esti_title)
            # h.name = "Lambda_over_pich__vs__mult"
            # h.write()

            # h = create_hist_pid_ratio_over_mult(h3d, [1], [2])
            # h.name = "lambda_over_K0S__vs__mult"
            # h.write()

# Create P(Nch) plots
with root_open(sys.argv[1], 'update') as f:
    log.info("Creating P(Nch) plots")
    for postfix in ["", "_unweighted"]:
        pNch_summary = HistStack()  # One HistStack!
        pNch_summary.title = "P(N_{ch}^{est}) summary "

        step_size = 40  # binning in Nch^est
        ref_ests = ['Total', 'EtaLt05']

        # make ntuples:
        nt0 = f.Sums[0].FindObject("fevent_counter")
        nt0.SetAlias(f.Sums[0].GetName(), "fevent_counter")
        previous_estimator_names = [f.Sums[0].GetName(), ]
        for est_dir in f.Sums[1:]:
            nt0.AddFriend(est_dir.FindObject("fevent_counter"), est_dir.GetName())
        for est_dir in f.Sums:
            ############################################################
            # Summary plot:
            name = "PNch_"+est_dir.GetName()
            axis_title = ";N_{ch}^{est}"
            h_summary = Hist1D(200, 0, 400, name=name,
                               title=est_dir.GetName() + axis_title)
            nt0.Project(h_summary.name,
                        "{}.nch".format(est_dir.GetName()))
            h_summary.Scale(1.0/h_summary.Integral())
            pNch_summary.Add(h_summary)

            #############################################################
            # Plot binned in mult_est vs. mutl_total
            for ref_est in ref_ests:
                axis_title = ";N_{{ch}}^{{{}}}".format(ref_est)
                pNch_per_est_vs_total = HistStack()
                hists = []
                pNch_per_est_vs_total.title = "P(N_{{ch}}^{}) for {}".format(ref_est, est_dir.GetName())
                for nch_est_min in range(0, 401, step_size):
                    nch_est_max = nch_est_min + step_size
                    name = "PNch_{}_lt_mult_tot_lt_{}".format(nch_est_min, nch_est_max)
                    h_tmp = Hist1D(50, 0, 100, name=name)
                    h_tmp.title = "{} < N_{{ch}}^{{est}} < {}".format(nch_est_min, nch_est_max) + axis_title
                    nt0.Project(h_tmp.name,
                                "{}.nch".format(est_dir.GetName()),
                                "{0} < {2}.nch && {2}.nch < {1}".format(nch_est_min, nch_est_max, ref_est))
                    try:
                        h_tmp.Scale(1.0/h_tmp.Integral())
                    except ZeroDivisionError:
                        pass
                    pNch_per_est_vs_total.Add(h_tmp)
                    hists.append(h_tmp)
                c = plot_histogram_stack(pNch_per_est_vs_total)
                c.name = "PNch_vs_mult_{}_for_{}".format(ref_est, est_dir.GetName())
                c.FindObject("plot").SetLogy(1)
                f.results_post.cd(est_dir.GetName() + postfix)
                c.write()
                del pNch_per_est_vs_total
                [h.Delete() for h in hists]
            
        # Write summary to disk
        c = plot_histogram_stack(pNch_summary)
        c.FindObject("plot").SetLogy(1)
        c.name = "PNch_summary" + postfix
        f.cd("results_post")
        c.write()


# Create ratio plots; depends on the previously created histograms
with root_open(sys.argv[1], 'update') as f:
    log.info("Creating ratios of dN/deta plots for each multiplicity bin")
    for postfix in ["", "_unweighted"]:
        # Get the  dN/deta stack for each estimator in order to calc 
        # the cannonical average for all _other_ estimators later on
        # P(N_ch): one stack containing all estimators!
        dNdeta_stacks = []        # List of stacks!
        pNch_stack = HistStack()  # One HistStack!
        pNch_stack.name = "P(N_{ch}^{est}) summary "
        for est_dir in f.Sums:
            dNdeta_stacks.append(asrootpy(f.Get("results_post")\
                                   .Get(est_dir.GetName() + postfix)\
                                   .Get("dNdeta_summary")\
                                   .FindObject('dNdeta_stack')))

        # looping over file again in order to have the estimator name handy,
        for i, est_dir in enumerate(f.Sums):
            res_dir_str = "results_post/" + est_dir.GetName() + postfix + "/ratios_to_other_est"
            try:
                f.mkdir(res_dir_str, recurse=True)
            except:
                pass
            # calc cannonical avg withouth the current hist
            other_dNdeta_stacks = dNdeta_stacks[:i] + dNdeta_stacks[i+1:]
            avg_stack = create_canonnical_avg_from_stacks(other_dNdeta_stacks)
            ratio = divide_stacks(dNdeta_stacks[i], avg_stack)
            ratio.title = ('Ratio of '
                           + dNdeta_stacks[i].title
                           + ' to cannonical average')
            
            c = plot_histogram_stack(ratio)
            c.name = dNdeta_stacks[i].name + '_ratio_cannonical_avg'
            c.Update()
            f.cd(res_dir_str)
            c.Write()

            # create ratio between P(N_ch) and cannonical_avg
            # avg = pNch_stack[0].Clone()
            # avg.Clear()
            # other_PNch = pNch_stack[:i] + pNch_stack[i+1:]
            # [avg.Add(h) for h in other_PNch]
            # avg.Scale(1.0/(len(other_PNch)))
                                            
            # ratio = pNch_stack[i] / avg
            # ratio.name = "{0}_div_by_can_avg".format(pNch_stack[i].name)
            # ratio.title = "Ratio of {} over cannonical avg".format(pNch_stack[i].title)
            # f.cd(res_dir_str)
            # ratio.Write()

# Make correlations between estimators
with root_open(sys.argv[1], 'update') as f:
    log.info("Correlating N_ch of each estimator")
    corr_dir = 'results_post/correlations'
    try:
        f.mkdir(corr_dir, recurse=True)
    except:
        pass
    # only for weighted case
    # Take ntuple from the first estimator and then add friends to this one
    nt0 = f.Sums[0].FindObject("fevent_counter")
    nt0.SetAlias(f.Sums[0].GetName(), "fevent_counter")
    previous_estimator_names = [f.Sums[0].GetName(), ]
    for est_dir in f.Sums[1:]:
        nt0.AddFriend(est_dir.FindObject("fevent_counter"), est_dir.GetName())
        for est_name in previous_estimator_names:
            corr_hist = Hist2D(200, 0, 200,
                               200, 0, 200,
                               name="corr_hist_{}_vs_{}".format(est_name, est_dir.GetName()))
            # Lables are deliberatly swaped, see Projection below!
            corr_hist.title = ("Correlation N_{{ch}} in {0} and {1};N_{{ch}} {1};N_{{ch}} {0}"\
                               .format(est_name, est_dir.GetName()))

            # this projects onto y:x, to make coding more adventurous 
            nt0.Project(corr_hist.name, "{0}.nch:{1}.nch".format(est_name, est_dir.GetName()),
                        "ev_weight")
            corr_hist.drawstyle = 'colz'
            f.cd(corr_dir)
            corr_hist.write()
        previous_estimator_names.append(est_dir.GetName())

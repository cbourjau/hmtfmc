"""Run post analysis on the given file"""
import sys
if len(sys.argv) != 2:
    print "Usage: python ./post.py path_to_root_file.root"
    quit()

from rootpy.io import root_open
from rootpy import asrootpy, ROOT, log
from rootpy.plotting import HistStack
from post_utils import create_dNdeta_stack,\
    plot_histogram_stack, create_stack_pid_ratio_over_pt,\
    create_hist_pid_ratio_over_mult,\
    create_canonnical_avg_from_stacks,\
    divide_stacks

# go into batch mode
ROOT.gROOT.SetBatch(True)

log = log["/post"]  # set name of this script in logger
log.info("IsBatch: {0}".format(ROOT.gROOT.IsBatch())) # Results in "DEBUG:myapp] Hello"

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
            h3d = f_post.Sums.FindObject(est_dir.GetName()).FindObject('festi_pT_pid' + postfix)
            h3d = asrootpy(h3d)

            h2d = f_post.Sums.FindObject(est_dir.GetName()).FindObject('fdNdeta' + postfix)
            h2d = asrootpy(h2d)

            h_event_counter = asrootpy(f_post.Sums\
                                       .FindObject(est_dir.GetName())\
                                       .FindObject('fEventcounter' + postfix))
            res_dir = f_post.mkdir("results_post/" + est_dir.GetName() + postfix, recurse=True)
            res_dir.write()

            # res_dir = f_post.results_post.FindObject(est_dir.GetName())
            f_post.results_post.cd(est_dir.GetName() + postfix)

            esti_title = "({0})".format(h3d.title[38:])
            ###########################################################
            # Category 1 on TWiki
            # create dN/deta stack for the current estimator
            hs = create_dNdeta_stack(h2d, h_event_counter)
            hs.title = "dN/d\\eta \\ vs.\\ \\eta\\ " + esti_title
            c = plot_histogram_stack(hs)
            c.name = "dNdeta_summary"
            c.Update()
            c.write()

            # P(N_ch)
            h_PN_ch = asrootpy(h_event_counter.clone("PN_ch"))
            h_PN_ch.Scale(1.0/h_PN_ch.Integral())
            h_PN_ch.title = "P(N_{ch})\\ " + esti_title
            h_PN_ch.yaxis.title = "P(N_{ch})"
            h_PN_ch.write()

            ###########################################################
            # Category 2 on TWiki
            # create particle ratio vs pT plots
            # Ratios to pich
            hs = create_stack_pid_ratio_over_pt(h3d, h_event_counter, [0], [5,6])
            hs.title = "p/\\pi^{+-}\\ vs.\\ p_{T} " + "({})".format(esti_title)
            c = plot_histogram_stack(hs)
            c.name = "proton_over_pich__vs__pt"
            c.Update()
            c.write()
            
            hs = create_stack_pid_ratio_over_pt(h3d, h_event_counter, [2], [5,6])
            hs.SetTitle(r'K^{0}_{S}/\\pi^{+-}\\ vs.\\ p_{T}\\ ' + "({})".format(esti_title))
            c = plot_histogram_stack(hs)
            c.name = "K0S_over_pich__vs__pt"
            c.Update()
            c.write()

            hs = create_stack_pid_ratio_over_pt(h3d, h_event_counter,[1], [5,6])
            hs.title= r"\\Lambda/\\pi^{+-}\\ vs.\\ p_{T}\\ " + "({})".format(esti_title)
            c = plot_histogram_stack(hs)
            c.name = "Lambda_over_pich__vs__pt"
            c.Update()
            c.write()

            hs = create_stack_pid_ratio_over_pt(h3d, h_event_counter,[8], [5,6])
            hs.title= "\\Xi/\\pi^{+-}\\ vs.\\ p_{T}\\ " + "({})".format(esti_title)
            c = plot_histogram_stack(hs)
            c.name = "Xi_over_pich__vs__pt"
            c.Update()
            c.write()

            hs = create_stack_pid_ratio_over_pt(h3d, h_event_counter,[9,10], [5,6])
            hs.title= "\Omega_{ch}/\pi^{+-} vs. p_{T} " + "({})".format(esti_title)
            c = plot_histogram_stack(hs)
            c.name = "OmegaCh_over_pich__vs__pt"
            c.Update()
            c.write()

            # Ratios to pi0
            # hs = create_stack_pid_ratio_over_pt(h3d, h_event_counter,[5,6], [7])
            # hs.title = "p^{+-}/\pi^{0} vs. p_{T} " + "({})".format(esti_title)
            # c = plot_histogram_stack(hs)
            # c.name = "pich_over_pi0__vs__pt"
            # c.Update()
            # c.write()

            # c = plot_histogram_stack(create_stack_pid_ratio_over_pt(h3d, h_event_counter,[0], [7]))
            # c.name = "proton_over_pi0__vs__pt"
            # c.title= "p/#pi^{0} vs. p_{T} " + "({})".format(est_dir.GetName())
            # c.Update()
            # c.write()

            # c = plot_histogram_stack(create_stack_pid_ratio_over_pt(h3d, h_event_counter,[2], [7]))
            # c.name = "K0S_over_pi0__vs__pt"
            # c.title= "K0S/#pi^{0} vs. p_{T} " + "({})".format(est_dir.GetName())
            # c.Update()
            # c.write()

            # c = plot_histogram_stack(create_stack_pid_ratio_over_pt(h3d, h_event_counter,[1], [7]))
            # c.name = "Lambda_over_pi0__vs__pt"
            # c.title= "#Lambda/#pi^{0} vs. p_{T} " + "({})".format(est_dir.GetName())
            #c.Update()
            # c.write()

            # c = plot_histogram_stack(create_stack_pid_ratio_over_pt(h3d, h_event_counter,[8], [7]))
            # c.name = "Xi_over_pi0__vs__pt"
            # c.title= "#Xi/#pi^{0} vs. p_{T} " + "({})".format(est_dir.GetName())
            # c.Update()
            # c.write()

            # c = plot_histogram_stack(create_stack_pid_ratio_over_pt(h3d, h_event_counter,[9,10], [7]))
            # c.name = "OmegaCh_over_pi0__vs__pt"
            # c.title= "#Omega_{ch}/#pi^{0} vs. p_{T} " + "({})".format(est_dir.GetName())
            # c.Update()
            # c.write()

            # Ratios to K0S
            hs = create_stack_pid_ratio_over_pt(h3d, h_event_counter,[0], [2])
            hs.title= "p/K^{0}_{S} vs. p_{T} " + "({})".format(esti_title)
            c = plot_histogram_stack(hs)
            c.name = "proton_over_K0S__vs__pt"
            c.Update()
            c.write()

            hs = create_stack_pid_ratio_over_pt(h3d, h_event_counter,[1], [2])
            hs.title= "\Lambda/K^{0}_{S} vs. p_{T} " + "({})".format(esti_title)
            c = plot_histogram_stack(hs)
            c.name = "Lambda_over_K0S__vs__pt"
            c.Update()
            c.write()

            hs = create_stack_pid_ratio_over_pt(h3d, h_event_counter,[8], [2])
            hs.title= "\Xi/K^{0}_{S} vs. p_{T} " + "({})".format(esti_title)
            c = plot_histogram_stack(hs)
            c.name = "Xi_over_K0S__vs__pt"
            c.Update()
            c.write()

            hs = create_stack_pid_ratio_over_pt(h3d, h_event_counter,[9,10], [2])
            hs.title= "\Omega_{ch}/K^{0}_{S} vs. p_{T} " + "({})".format(esti_title)
            c = plot_histogram_stack(hs)
            c.name = "OmegaCh_over_K0S__vs__pt"
            c.Update()
            c.write()
            #############################################################33
            #Category 3 on Twiki
            h = create_hist_pid_ratio_over_mult(h3d, h_event_counter,[0], [5,6])
            h.name = "proton_over_pich__vs__mult"
            h.write()

            h = create_hist_pid_ratio_over_mult(h3d, h_event_counter,[1], [2])
            h.name = "lambda_over_K0S__vs__mult"
            h.write()



# Create ratio plots; depends on the previously created histograms
with root_open(sys.argv[1], 'update') as f:
    log.info("Creating ratios of dN/deta plots for each multiplicity bin")
    for postfix in ["", "_unweighted"]:
        # Get the  dN/deta stack for each estimator in order to calc 
        # the cannonical average for all _other_ estimators later on
        # P(N_ch): one stack containing all estimators!
        dNdeta_stacks = []        # List of stacks!
        pNch_stack = HistStack()  # One HistStack!
        pNch_stack.name = "pNch_stack"
        for est_dir in f.Sums:
            dNdeta_stacks.append(asrootpy(f.Get("results_post")\
                                   .Get(est_dir.GetName() + postfix)\
                                   .Get("dNdeta_summary")\
                                   .FindObject('dNdeta_stack')))
            pNch_stack.Add(asrootpy(f.Get("results_post")\
                                    .Get(est_dir.GetName() + postfix)\
                                    .Get("PN_ch")))

        # Write P(Nch) stack to disk before proceeding
        pNch_stack.title = "P(N_{ch}) of various estimators"
        c = plot_histogram_stack(pNch_stack)
        c.name = "PN_ch_summary" + postfix
        c.Update()
        f.cd("results_post")
        c.write()

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
            ratio.title = (r'Ratio\\ of\\ '
                           + dNdeta_stacks[i].title
                           + r'\\ to\\ cannonical\\ average')
            
            c = plot_histogram_stack(ratio)
            c.name = dNdeta_stacks[i].name + '_ratio_cannonical_avg'
            c.Update()
            f.cd(res_dir_str)
            c.Write()

            # create ratio between P(N_ch) and cannonical_avg
            avg = pNch_stack[0].Clone()
            avg.Clear()
            other_PNch = pNch_stack[:i] + pNch_stack[i+1:]
            [avg.Add(h) for h in other_PNch]
            avg.Scale(1.0/(len(other_PNch)))
                                            
            ratio = pNch_stack[i] / avg
            ratio.name = "{0}_div_by_can_avg".format(pNch_stack[i].name)
            ratio.title = "Ratio\\ of\\ {}\\ over\\ cannonical\\ avg".format(pNch_stack[i].title)
            f.cd(res_dir_str)
            ratio.Write()

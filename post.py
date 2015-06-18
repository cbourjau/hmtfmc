"""Run post analysis on the given file"""
import sys

from rootpy.io import root_open
from rootpy import asrootpy
from rootpy.plotting import HistStack, Canvas, Legend, Pad, Hist1D
from post_utils import create_dNdeta_stack, make_stack_of_mult_bins_for_pids,\
    plot_stack_of_estimators, create_stack_pid_ratio_over_pt,\
    create_hist_pid_ratio_over_mult

if len(sys.argv) != 2:
    print "Usage: python ./post.py path_to_root_file.root"
    exit

f = root_open(sys.argv[1], 'r')

f_post = root_open(sys.argv[1], 'update')
try:
    f_post.rmdir('Results_post')
except:
    pass

# Loop over all estimators in the Sums list:
for est_dir in f.Sums:
    # and do everything for weighted and unweighted:
    for postfix in ["", "_unweighted"]:
        h3d = f.Sums.FindObject(est_dir.GetName()).FindObject('festi_pT_pid' + postfix)
        h3d = asrootpy(h3d)

        h2d = f.Sums.FindObject(est_dir.GetName()).FindObject('fdNdeta' + postfix)
        h2d = asrootpy(h2d)
        res_dir = f_post.mkdir("Results_post/" + est_dir.GetName() + postfix, recurse=True)
        res_dir.write()

        # res_dir = f_post.Results_post.FindObject(est_dir.GetName())
        f_post.Results_post.cd(est_dir.GetName() + postfix)

        ###########################################################
        # Category 1 on TWiki
        # create dN/deta stack for the current estimator
        hs = create_dNdeta_stack(h2d)
        hs.title = "$dN/d\eta$ vs. $\eta$ " + "({})".format(h3d.title[30:])
        c = plot_stack_of_estimators(hs)
        c.name = "dN_deta_summary"
        c.Update()
        c.write()

        ###########################################################
        # Category 2 on TWiki
        # create particle ratio vs pT plots
        # Ratios to pich
        hs = create_stack_pid_ratio_over_pt(h3d, [0], [5,6])
        hs.title = "p/#pi^{+-} vs. multiplicity " + "({})".format(h3d.title[30:])
        c = plot_stack_of_estimators(hs)
        c.name = "proton_over_pich__vs__pt"
        c.Update()
        c.write()

        hs = create_stack_pid_ratio_over_pt(h3d, [2], [5,6])
        hs.title= "K0S/#pi^{+-} vs. multiplicity " + "({})".format(h3d.title[30:])
        c = plot_stack_of_estimators(hs)
        c.name = "K0S_over_pich__vs__pt"
        c.Update()
        c.write()

        hs = create_stack_pid_ratio_over_pt(h3d, [1], [5,6])
        hs.title= "#Lambda/#pi^{+-} vs. multiplicity " + "({})".format(h3d.title[30:])
        c = plot_stack_of_estimators(hs)
        c.name = "Lambda_over_pich__vs__pt"
        c.Update()
        c.write()

        hs = create_stack_pid_ratio_over_pt(h3d, [8], [5,6])
        hs.title= "#Xi/#pi^{+-} vs. multiplicity " + "({})".format(h3d.title[30:])
        c = plot_stack_of_estimators(hs)
        c.name = "Xi_over_pich__vs__pt"
        c.Update()
        c.write()

        hs = create_stack_pid_ratio_over_pt(h3d, [9,10], [5,6])
        hs.title= "#Omega_{ch}/#pi^{+-} vs. multiplicity " + "({})".format(h3d.title[30:])
        c = plot_stack_of_estimators(hs)
        c.name = "OmegaCh_over_pich__vs__pt"
        c.Update()
        c.write()

        # Ratios to pi0
        # c = plot_stack_of_estimators(create_stack_pid_ratio_over_pt(h3d, [5,6], [7]))
        # c.name = "pich_over_pi0__vs__pt"
        # c.title= "p^{+-}/#pi^{0} vs. multiplicity " + "({})".format(est_dir.GetName())
        # c.Update()
        # c.write()

        # c = plot_stack_of_estimators(create_stack_pid_ratio_over_pt(h3d, [0], [7]))
        # c.name = "proton_over_pi0__vs__pt"
        # c.title= "p/#pi^{0} vs. multiplicity " + "({})".format(est_dir.GetName())
        # c.Update()
        # c.write()

        # c = plot_stack_of_estimators(create_stack_pid_ratio_over_pt(h3d, [2], [7]))
        # c.name = "K0S_over_pi0__vs__pt"
        # c.title= "K0S/#pi^{0} vs. multiplicity " + "({})".format(est_dir.GetName())
        # c.Update()
        # c.write()

        # c = plot_stack_of_estimators(create_stack_pid_ratio_over_pt(h3d, [1], [7]))
        # c.name = "Lambda_over_pi0__vs__pt"
        # c.title= "#Lambda/#pi^{0} vs. multiplicity " + "({})".format(est_dir.GetName())
        #c.Update()
        # c.write()

        # c = plot_stack_of_estimators(create_stack_pid_ratio_over_pt(h3d, [8], [7]))
        # c.name = "Xi_over_pi0__vs__pt"
        # c.title= "#Xi/#pi^{0} vs. multiplicity " + "({})".format(est_dir.GetName())
        # c.Update()
        # c.write()

        # c = plot_stack_of_estimators(create_stack_pid_ratio_over_pt(h3d, [9,10], [7]))
        # c.name = "OmegaCh_over_pi0__vs__pt"
        # c.title= "#Omega_{ch}/#pi^{0} vs. multiplicity " + "({})".format(est_dir.GetName())
        # c.Update()
        # c.write()

        # Ratios to K0S
        hs = create_stack_pid_ratio_over_pt(h3d, [0], [2])
        hs.title= "p/#pi^{+} vs. multiplicity " + "({})".format(h3d.title[30:])
        c = plot_stack_of_estimators(hs)
        c.name = "proton_over_K0S__vs__pt"
        c.Update()
        c.write()

        hs = create_stack_pid_ratio_over_pt(h3d, [1], [2])
        hs.title= "#Lambda/#pi^{-} vs. multiplicity " + "({})".format(h3d.title[30:])
        c = plot_stack_of_estimators(hs)
        c.name = "Lambda_over_K0S__vs__pt"
        c.Update()
        c.write()

        hs = create_stack_pid_ratio_over_pt(h3d, [8], [2])
        hs.title= "#Xi/#pi^{-} vs. multiplicity " + "({})".format(h3d.title[30:])
        c = plot_stack_of_estimators(hs)
        c.name = "Xi_over_K0S__vs__pt"
        c.Update()
        c.write()

        hs = create_stack_pid_ratio_over_pt(h3d, [9,10], [2])
        hs.title= "#Omega_{ch}/#pi^{-} vs. multiplicity " + "({})".format(h3d.title[30:])
        c = plot_stack_of_estimators(hs)
        c.name = "OmegaCh_over_K0S__vs__pt"
        c.Update()
        c.write()
        #############################################################33
        #Category 3 on Twiki
        h = create_hist_pid_ratio_over_mult(h3d, [0], [5,6])
        h.name = "proton_over_pich__vs__mult"
        h.write()

        h = create_hist_pid_ratio_over_mult(h3d, [1], [2])
        h.name = "lambda_over_K0S__vs__mult"
        h.write()
f_post.close()

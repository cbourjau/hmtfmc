from rootpy.io import root_open
from rootpy import asrootpy
from rootpy.plotting import HistStack, Canvas, Legend, Pad, Hist1D
from post_utils import make_stack_of_mult_bins_for_pids,\
    plot_stack_of_estimators, create_stack_pid_ratio_over_pt,\
    create_hist_pid_ratio_over_mult

f = root_open("../caf/MC_estimators.root", 'r')
h3d = f.Sums.FindObject('EtaLt05').FindObject('festi_pT_pid')
h3d = asrootpy(h3d)

f_post = root_open("../caf/MC_estimators.root", 'update')

for est_dir in f.Sums:
    try:
        f_post.rmdir('Results_post')
    except:
        pass
    res_dir = f_post.mkdir("Results_post/" + est_dir.GetName(), recurse=True)
    res_dir.write()

    # res_dir = f_post.Results_post.FindObject(est_dir.GetName())
    f_post.Results_post.cd(est_dir.GetName())
    ###########################################################
    # Category 2 on TWiki
    # create particle ratio vs pT plots
    # Ratios to pich
    h = plot_stack_of_estimators(create_stack_pid_ratio_over_pt(h3d, [0], [5,6]))
    h.name = "proton_over_pich__vs__pt"
    h.title= "p/#pi^{+} vs. multiplicity " + "({})".format(est_dir.GetName())
    h.write()
    
    h = plot_stack_of_estimators(create_stack_pid_ratio_over_pt(h3d, [2], [5,6]))
    h.name = "K0S_over_pich__vs__pt"
    h.title= "K0S/#pi^{-} vs. multiplicity " + "({})".format(est_dir.GetName())
    h.write()
    
    h = plot_stack_of_estimators(create_stack_pid_ratio_over_pt(h3d, [1], [5,6]))
    h.name = "Lambda_over_pich__vs__pt"
    h.title= "#Lambda/#pi^{-} vs. multiplicity " + "({})".format(est_dir.GetName())
    h.write()
    
    h = plot_stack_of_estimators(create_stack_pid_ratio_over_pt(h3d, [8], [5,6]))
    h.name = "Xi_over_pich__vs__pt"
    h.title= "#Xi/#pi^{-} vs. multiplicity " + "({})".format(est_dir.GetName())
    h.write()
    
    h = plot_stack_of_estimators(create_stack_pid_ratio_over_pt(h3d, [9,10], [5,6]))
    h.name = "OmegaCh_over_pich__vs__pt"
    h.title= "#Omega_{ch}/#pi^{-} vs. multiplicity " + "({})".format(est_dir.GetName())
    h.write()

    # Ratios to pi0
    # h = plot_stack_of_estimators(create_stack_pid_ratio_over_pt(h3d, [5,6], [7]))
    # h.name = "pich_over_pi0__vs__pt"
    # h.title= "p^{+-}/#pi^{0} vs. multiplicity " + "({})".format(est_dir.GetName())
    # h.write()
    
    # h = plot_stack_of_estimators(create_stack_pid_ratio_over_pt(h3d, [0], [7]))
    # h.name = "proton_over_pi0__vs__pt"
    # h.title= "p/#pi^{0} vs. multiplicity " + "({})".format(est_dir.GetName())
    # h.write()
    
    # h = plot_stack_of_estimators(create_stack_pid_ratio_over_pt(h3d, [2], [7]))
    # h.name = "K0S_over_pi0__vs__pt"
    # h.title= "K0S/#pi^{0} vs. multiplicity " + "({})".format(est_dir.GetName())
    # h.write()
    
    # h = plot_stack_of_estimators(create_stack_pid_ratio_over_pt(h3d, [1], [7]))
    # h.name = "Lambda_over_pi0__vs__pt"
    # h.title= "#Lambda/#pi^{0} vs. multiplicity " + "({})".format(est_dir.GetName())
    # h.write()
    
    # h = plot_stack_of_estimators(create_stack_pid_ratio_over_pt(h3d, [8], [7]))
    # h.name = "Xi_over_pi0__vs__pt"
    # h.title= "#Xi/#pi^{0} vs. multiplicity " + "({})".format(est_dir.GetName())
    # h.write()
    
    # h = plot_stack_of_estimators(create_stack_pid_ratio_over_pt(h3d, [9,10], [7]))
    # h.name = "OmegaCh_over_pi0__vs__pt"
    # h.title= "#Omega_{ch}/#pi^{0} vs. multiplicity " + "({})".format(est_dir.GetName())
    # h.write()
    
    # Ratios to K0S
    h = plot_stack_of_estimators(create_stack_pid_ratio_over_pt(h3d, [0], [2]))
    h.name = "proton_over_K0S__vs__pt"
    h.title= "p/#pi^{+} vs. multiplicity " + "({})".format(est_dir.GetName())
    h.write()
    
    h = plot_stack_of_estimators(create_stack_pid_ratio_over_pt(h3d, [1], [2]))
    h.name = "Lambda_over_K0S__vs__pt"
    h.title= "#Lambda/#pi^{-} vs. multiplicity " + "({})".format(est_dir.GetName())
    h.write()
    
    h = plot_stack_of_estimators(create_stack_pid_ratio_over_pt(h3d, [8], [2]))
    h.name = "Xi_over_K0S__vs__pt"
    h.title= "#Xi/#pi^{-} vs. multiplicity " + "({})".format(est_dir.GetName())
    h.write()
    
    h = plot_stack_of_estimators(create_stack_pid_ratio_over_pt(h3d, [9,10], [2]))
    h.name = "OmegaCh_over_K0S__vs__pt"
    h.title= "#Omega_{ch}/#pi^{-} vs. multiplicity " + "({})".format(est_dir.GetName())
    h.write()
    #############################################################33
    #Category 3 on Twiki
    c = create_hist_pid_ratio_over_mult(h3d, [0], [5,6])
    c.name = "proton_over_pich__vs__mult"
    c.write()

    c = create_hist_pid_ratio_over_mult(h3d, [1], [2])
    c.name = "lambda_over_K0S__vs__mult"
    c.write()
f_post.close()

"""
This is the entry point to run the post analysis. The file post_plotting contains the logic for the visual
presentation of the plots. The logic which extracts the data from the primary generated "Sums" but is not concerned with any visual representation can be found in post_data_extractor.py. Lastly, the file post_utils.py contains small helper function for cleaning up data or generating names and titles.
"""

import sys

from rootpy import log, ROOT

from post_plotting import Plotting

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage: python post.py file.root {Inel, InelGt0, V0AND}"
        quit()
    # go into batch mode
    ROOT.gROOT.SetBatch(True)

    log = log["/post"]  # set name of this script in logger
    log.info("IsBatch: {0}".format(ROOT.gROOT.IsBatch()))

    try:
        global_trigger = sys.argv[2]
    except IndexError:
        global_trigger = ""
    sums_dir_name = "Sums" + global_trigger
    results_dir_name = "results_post" + global_trigger

    plotting = Plotting(f_name=sys.argv[1], sums_dir_name=sums_dir_name, results_dir_name=results_dir_name)

    # run the actual plots:
    plotting.plot_dNdetas()
    plotting.plot_PNch()
    plotting.plot_mult_vs_pt()
    plotting.plot_meanpt_vs_ref_mult_for_pids()
    plotting.plot_pt_distribution_ratios()
    plotting.plot_pid_ratio_vs_refmult()
    plotting.plot_dNdpT()
    plotting.plot_pT_HM_div_pt_MB(scale_nMPI=False)
    plotting.plot_pT_HM_div_pt_MB(scale_nMPI=True)

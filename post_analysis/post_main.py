"""
This is the entry point to run the post analysis. The file post_plotting contains the logic for the visual
presentation of the plots. The logic which extracts the data from the primary generated "Sums" but is not concerned with any visual representation can be found in post_data_extractor.py. Lastly, the file post_utils.py contains small helper function for cleaning up data or generating names and titles.
"""

import sys

from rootpy import log, ROOT

from post_plotting import Plotting
from roofie.beamify import Beamerdoc

############
# Settings #
############

# A list of the estimators which will be included in the plots (if available in the input file)
considered_ests = ['EtaLt05', 'EtaLt08', 'EtaLt15', 'Eta08_15', 'V0M',  # 'V0A', 'V0C',
                   'ZDC', 'nMPI', 'Q2', 'spherocity', 'sphericity']

# Ranges of percentiles which should be considered in the plots These
# ranges are then translated into ranges of bins in multiplicity, nMPI
# or whatever is applicable
std_perc_bins = [(1, 0.7), (.5, .4), (.1, .05), (0.001, 0.0)]
percentile_bins = {
    'EtaLt05': std_perc_bins,
    'EtaLt08': std_perc_bins,
    'EtaLt15': std_perc_bins,
    'Eta08_15': std_perc_bins,
    'V0M': std_perc_bins,
    'V0A': std_perc_bins,
    'V0C': std_perc_bins,
    'ZDC': [(1, 0.7), (.7, .3), (.3, .05), (0.001, 0.0)],
    'nMPI': [(1, 0.7), (.7, .4), (.3, .05), (0.001, 0.0)],
    'Q2': [(1, 0.7), (.7, .4), (.3, .05), (0.001, 0.0)],
    'spherocity': [(1, 0.7), (.7, .4), (.3, .05), (0.001, 0.0)],
    'sphericity': [(1, 0.7), (.7, .4), (.3, .05), (0.001, 0.0)],
}

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print """Usage: python post.py file.root {Inel, InelGt0, V0AND} "<summary name>" """
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

    plotting = Plotting(f_name=sys.argv[1], sums_dir_name=sums_dir_name, results_dir_name=results_dir_name,
                        percentile_bins=percentile_bins, considered_ests=considered_ests)
    latexdoc = Beamerdoc(author="Christian Bourjau", title=sys.argv[3])

    # run the actual plots:
    sec = latexdoc.add_section(r"$dN/d\eta$")
    [sec.add_figure(fig) for fig in plotting.plot_dNdetas(ratio_to_mb=False)]

    sec = latexdoc.add_section(r"$dN/d\eta$ over MB result")
    [sec.add_figure(fig) for fig in plotting.plot_dNdetas(ratio_to_mb=True)]

    sec = latexdoc.add_section(r"$P(N_{ch})$ summary")
    [sec.add_figure(fig) for fig in plotting.plot_PNch_summary()]

    sec = latexdoc.add_section(r"$P(N_{ch})$")
    [sec.add_figure(fig) for fig in plotting.plot_PNch()]

    plotting.plot_mult_vs_pt()

    sec = latexdoc.add_section(r"$\left< p_T \right>$ vs. ref multiplicity")
    [sec.add_figure(fig) for fig in plotting.plot_meanpt_vs_ref_mult_for_pids()]

    sec = latexdoc.add_section(r"Ratios for various species vs $p_T$")
    [sec.add_figure(fig) for fig in plotting.plot_pt_distribution_ratios()]

    sec = latexdoc.add_section(r"Ratios for various species vs ref. multiplicity")
    [sec.add_figure(fig) for fig in plotting.plot_pid_ratio_vs_refmult()]

    sec = latexdoc.add_section(r"$dN/dp_T$")
    [sec.add_figure(fig) for fig in plotting.plot_dNdpT()]

    sec = latexdoc.add_section(r"$\left[ dN_{HM}/dp_T\right] / \left[ dN_{MB}/dp_T\right]$")
    [sec.add_figure(fig) for fig in plotting.plot_pT_HM_div_pt_MB(scale_nMPI=False)]

    sec = latexdoc.add_section(r"$\left[ dN_{HM}/dp_T\right] / \left[ dN_{MB}/dp_T\right] \times \left[ \left<N_{MPI}^{MB}\right> / \left<N_{MPI}^{HM}\right>\right]$")
    [sec.add_figure(fig) for fig in plotting.plot_pT_HM_div_pt_MB(scale_nMPI=True)]

    sec = latexdoc.add_section(r"$nMPI(N_{ch})$")
    [sec.add_figure(fig) for fig in plotting.plot_nMPI_vs_Nch()]

    latexdoc.finalize_document()

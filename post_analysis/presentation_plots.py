from rootpy.io import root_open
from rootpy.plotting import Graph, Hist1D

from roofie.figure import Figure, get_color_generator

from post_utils import make_estimator_title

triggers = ['Inel', 'InelGt0', 'V0AND']

# !!! The two generator lists have to have the same order!
generators = ['Pythia6_MB_7TeV_Perugia0',
              'Pythia6_MB_7TeV_Perugia2011',
              'Pythia6_MB_noCR',
              'Pythia8_MB_7TeV_Monash_noCR',
              'Pythia8_MB_7TeV_Monash_CR',
              'Pythia6_MB_CR',
              'Dipsy']

generators_short = ['Py6 Per0',
                    'Py6 Per2011',
                    'Py6 MB noCR',
                    'Py8 Mon. noCR',
                    'Py8 Mon. CR',
                    'Py6 MB CR',
                    'Dipsy']

estimators = ['EtaLt05', 'Eta08_15', 'V0M']  # 'sphericity']


def get_generator(name):
    """
    Return TDirectory representing this generator
    """
    fs = {'Pythia6_MB_7TeV_Perugia0':
          '368_Pythia6_MB_7TeV_Perugia0/AnalysisResults.root',
          'Pythia6_MB_7TeV_Perugia2011':
          '369_Pythia6_MB_7TeV_Perugia2011/AnalysisResults.root',
          'Pythia8_MB_7TeV_Monash_CR':
          '371_Pythia8_MB_7TeV_Monash_CR/AnalysisResults.root',
          'Pythia8_MB_7TeV_Monash_noCR':
          '372_Pythia8_MB_7TeV_Monash_noCR/AnalysisResults.root',
          'Dipsy':
          '374_Dipsy/AnalysisResults.root',
          'Pythia6_MB_CR':
          '375_Pythia6_MB_CR/AnalysisResults.root',
          'Pythia6_MB_noCR':
          '376_Pythia6_MB_noCR/AnalysisResults.root'}

    if name in fs.keys():
        f = root_open(fs.get(name), 'read')
        tdir = f.MultEstimators
        return tdir


def get_trigger(name, gen):
    """Return the TList/TDirectory with for a given trigger from the given
    generator file"""
    dir_name = 'results_post' + name
    return gen.__getattr__(dir_name)


def get_estimator(name, trigger):
    """
    Return the TList/TDirectory for an estimator from a given trigger TList.
    """
    return trigger.__getattr__(name)


summary_fig = Figure()
# summary_fig.plot.ymin = 4
summary_fig.ytitle = r"""M_{0 - 0.1%}/<M> avg'd over |#eta| < 0.8"""
summary_fig.plot.gridx = True
summary_fig.plot.xmin = 0
summary_fig.plot.ymin = 2

summary_fig.legend.position = 'seperate'
colors = list(get_color_generator(palette='colorblind'))
eta_range = (-0.8, 0.8)  # to integrate over

# hack around to name the axis bins
# leave some 1 bins room for the legend
hist_labels = Hist1D(len(generators) + 1, 0, len(generators) + 1)
[hist_labels.xaxis.SetBinLabel(n + 1, gen_name)
 for n, gen_name in enumerate(generators_short)]
summary_fig.add_plottable(hist_labels, use_as_frame=True)

# Add the lables for the estimator colors
[summary_fig.add_plottable(hist_labels, legend_title=make_estimator_title(est),
                           markerstyle=34, color=colors[nest])
 for nest, est in enumerate(estimators)]

# Add the lables for the trigger markers
[summary_fig.add_plottable(hist_labels, legend_title=trig, markerstyle=20 + ntrig, color=1)
 for ntrig, trig in enumerate(triggers)]

for ngen, gen_name in enumerate(generators):
    g = get_generator(gen_name)
    for ntrig, trig in enumerate(triggers):
        for nest, estname in enumerate(estimators):
            est = get_estimator(estname, get_trigger(trig, g))
            tmp_fig = Figure()
            tmp_fig.import_plottables_from_canvas(est.dNdeta_MB_ratio_summary)
            # Find the HM histogram
            for obj in tmp_fig._plottables:
                if obj['legend_title'] == '0.0%-0.1%':
                    break
            obj['p'].xaxis.SetRangeUser(*eta_range)
            nbins_in_range = (obj['p'].xaxis.FindBin(eta_range[1])
                              - obj['p'].xaxis.FindBin(eta_range[0]))
            tmp_graph = Graph()  # one bin per generator
            tmp_graph.SetPoint(0, ngen + 0.25 * (nest + 1),
                               obj['p'].Integral() / nbins_in_range)
            if ngen == 0:
                summary_fig.add_plottable(tmp_graph,
                                          # legend_title=trig,
                                          color=colors[nest],
                                          markerstyle=20 + ntrig)
            else:
                summary_fig.add_plottable(tmp_graph, color=colors[nest],
                                          markerstyle=20 + ntrig)
#summary_fig.save_to_file("./", "m_over_mean_m.pdf")

###############
# Omega / pi0 #
###############
generators = ['Pythia6_MB_7TeV_Perugia2011',
              'Pythia8_MB_7TeV_Monash_noCR',
              'Pythia8_MB_7TeV_Monash_CR',
              'Dipsy']

generators_short = ['Py6 Per2011',
                    'Py8 Mon. noCR',
                    'Py8 Mon. CR',
                    'Dipsy']
omega_pi0_fig = Figure()

for ngen, gen_name in enumerate(generators):
    g = get_generator(gen_name)
    tmp_fig = Figure()
    tmp_fig.import_plottables_from_canvas(
        get_trigger('V0AND', g).pid_ratios_vs_refmult.__getattr__('3122_div_111')
    )
    # Find the V0M histogram
    for obj in tmp_fig._plottables:
        if obj['legend_title'] == 'V0M':
            break
    omega_pi0_fig.add_plottable(obj['p'], legend_title=generators_short[ngen])

omega_pi0_fig.plot.xmax = 50
omega_pi0_fig.plot.ymin = 0.005
omega_pi0_fig.plot.ymax = 0.06
omega_pi0_fig.plot.palette = 'colorblind'
omega_pi0_fig.legend.title = 'Selected in V0M'
omega_pi0_fig.xtitle = r'N_{ch}^{|#eta|<0.5}'
omega_pi0_fig.ytitle = r'#Lambda / #pi^{0}'

#omega_pi0_fig.save_to_file("./", "lambda_pi0_refmult_V0M.pdf")

##############################
# HM dN/dpT normalized to MB #
##############################
generators = ['Pythia6_MB_noCR',
              'Pythia6_MB_CR',
              'Pythia8_MB_7TeV_Monash_noCR',
              'Pythia8_MB_7TeV_Monash_CR',
              'Dipsy']

generators_short = ['Py6 MB noCR',
                    'Py6 MB CR',
                    'Py8 Mon. noCR',
                    'Py8 Mon. CR',
                    'Dipsy']
estimator = 'V0M'
trigger = 'Inel'
dndpt_to_mb = Figure()
for ngen, gen_name in enumerate(generators):
    g = get_generator(gen_name)
    tmp_fig = Figure()
    tmp_fig.import_plottables_from_canvas(
        get_estimator(estimator, get_trigger(trigger, g)).pt_hm_div_pt_mb
    )
    # Find the HM histogram
    for obj in tmp_fig._plottables:
        if obj['legend_title'] == '0.0%-0.1%':
            break
    dndpt_to_mb.add_plottable(obj['p'], legend_title=generators_short[ngen])

dndpt_to_mb.plot.ymin = 2
dndpt_to_mb.plot.ymax = 80 if estimator == 'EtaLt05' else 20
dndpt_to_mb.plot.palette = 'colorblind'
dndpt_to_mb.xtitle = "p_{T} (GeV)"
dndpt_to_mb.ytitle = "(dN^{HM}/dp_{T}) / (dN^{MB}/dp_{T})"
dndpt_to_mb.legend.title = estimator + " and " + trigger
#dndpt_to_mb.save_to_file("./", "dndpt_to_mb_{}_{}.pdf".format(estimator, trigger))

###################
# dN/deta summary #
###################
generators = [
    'Pythia8_MB_7TeV_Monash_noCR',
    'Pythia6_MB_noCR',
    'Pythia8_MB_7TeV_Monash_CR',
    'Pythia6_MB_CR',
    'Pythia6_MB_7TeV_Perugia0',
    'Dipsy'
]

generators_short = [
    'Py8 Mon. noCR',
    'Py6 MB noCR',
    'Py8 Mon. CR',
    'Py6 MB CR',
    'Py6 Per0',
    'Dipsy'
]
estimator = 'V0M'
trigger = 'V0AND'
dndeta = Figure()
for ngen, gen_name in enumerate(generators):
    g = get_generator(gen_name)
    tmp_fig = Figure()
    tmp_fig.import_plottables_from_canvas(
        get_estimator(estimator, get_trigger(trigger, g)).dNdeta_summary
    )
    # Find the HM histogram
    for obj in tmp_fig._plottables:
        if obj['legend_title'] == '0.0%-0.1%':
            break
    dndeta.add_plottable(obj['p'], legend_title=generators_short[ngen])

dndeta.plot.palette = 'colorblind'
dndeta.legend.title = "0-0.1% in V0M; V0AND"
dndeta.xtitle = "#eta"
dndeta.ytitle = r'dN_{ch}/d#eta'
dndeta.plot.ymax = 75
#dndeta.save_to_file("./", "hm_dNdeta_gen_comp.pdf")

##############################
# P(Nch)_eta<0.5 in V0M bins #
##############################
generators = [
    'Pythia8_MB_7TeV_Monash_noCR',
    'Pythia6_MB_noCR',
    'Pythia8_MB_7TeV_Monash_CR',
    'Pythia6_MB_CR',
    'Pythia6_MB_7TeV_Perugia0',
    'Dipsy'
]

generators_short = [
    'Py8 Mon. noCR',
    'Py6 MB noCR',
    'Py8 Mon. CR',
    'Py6 MB CR',
    'Py6 Per0',
    'Dipsy'
]
estimator = 'V0M'
trigger = 'V0AND'
pnch = Figure()
for ngen, gen_name in enumerate(generators):
    g = get_generator(gen_name)
    tmp_fig = Figure()
    tmp_fig.import_plottables_from_canvas(
        get_estimator(estimator, get_trigger(trigger, g)).PNchEtaLt05_binned_in_NchEst
    )
    # Find the HM histogram
    for obj in tmp_fig._plottables:
        if obj['legend_title'] == '0.0%-0.1%':
            break
    pnch.add_plottable(obj['p'], legend_title=generators_short[ngen])

pnch.plot.palette = 'colorblind'
pnch.xtitle = r"N_{{ch}}^{{{}}}".format(make_estimator_title("EtaLt05"))
pnch.ytitle = r"P(N_{{ch}}), scaled to N_{{ch}}^{{{}}}".format(make_estimator_title("EtaLt05"))
pnch.legend.position = 'tr'
pnch.plot.xmax = 90
pnch.legend.title = "V0M, V0AND"
pnch.save_to_file("./", "hm_pnch_gen_comp.pdf")

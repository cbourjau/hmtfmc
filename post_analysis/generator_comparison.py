from roofi.figure import Figure

from rootpy.io import root_open

colors = ([(0.0, 0.4470588235294118, 0.6980392156862745),
           (0.0, 0.6196078431372549, 0.45098039215686275),
           (0.8352941176470589, 0.3686274509803922, 0.0),
           (0.8, 0.4745098039215686, 0.6549019607843137)])

fs = [(root_open('train_out/259_Pythia6_MB_7TeV_Perugia2011/AnalysisResults.root', 'read'), 'Py.6 Perugia2011'),
      (root_open('train_out/262_Pythia8_MB_7TeV_Monash_CR/AnalysisResults.root', 'read'), 'Py.8 Monash CR'),
      (root_open('train_out/263_Pythia8_MB_7TeV_Monash_noCR/AnalysisResults.root', 'read'), 'Py.8 Monash no CR'),
      # (root_open('train_out/222_Pythia8_MB_7TeV_CR/AnalysisResults.root', 'read'), 'Py.8 MB CR'),
      # Virtually identical to pythia6 perugia 2011:
      # (root_open('train_out/223_Pythia8_MB_7TeV_noCR/AnalysisResults.root', 'read'), 'Py.8 MB no CR')
]
fout = root_open('hmtf_pres_plots.root', 'recreate')


###########
# dN/deta #
###########

fig = Figure()
fig.plot.palette = 'colorblind'
fig.xtitle = r'#eta'
fig.ytitle = r'1/N_{evts} dN/d#eta'
fig.plot.ymin = 0
fig.plot.ymax = 69
hms = []
mbs = []
for f in fs:
    c = f[0].MultEstimators.results_post.V0M.dNdeta_summary
    tmp_fig = Figure()
    tmp_fig.import_plottables_from_canvas(c)
    # tmp_fig._plottables = tmp_fig._plottables[-2:]
    hms.append(tmp_fig._plottables[-2])
    mbs.append(tmp_fig._plottables[-1])
    fig.legend.title = '{} and {} (V0M)'.format(hms[-1]['legend_title'], mbs[-1]['legend_title'])
    hms[-1]['legend_title'] = f[1]
for hm, mb, color in zip(hms, mbs, colors):
    fig.add_plottable(hm['p'], legend_title=hm['legend_title'], color=color)
    fig.add_plottable(mb['p'], color=color, markerstyle=25)
fig.save_to_root_file(fout, 'dNdeta_gen_comp_V0M')

fig = Figure()
fig.plot.palette = 'colorblind'
fig.xtitle = r'#eta'
fig.ytitle = r'1/N_{evts} dN/d#eta'
fig.plot.ymin = 0
fig.plot.ymax = 69
hms = []
mbs = []
for f in fs:
    c = f[0].MultEstimators.results_post.EtaLt15.dNdeta_summary
    tmp_fig = Figure()
    tmp_fig.import_plottables_from_canvas(c)
    # tmp_fig._plottables = tmp_fig._plottables[-2:]
    hms.append(tmp_fig._plottables[-2])
    mbs.append(tmp_fig._plottables[-1])
    fig.legend.title = '{} and {} (|#eta|#leq1.5)'.format(hms[-1]['legend_title'], mbs[-1]['legend_title'])
    hms[-1]['legend_title'] = f[1]
for hm, mb, color in zip(hms, mbs, colors):
    fig.add_plottable(hm['p'], legend_title=hm['legend_title'], color=color)
    fig.add_plottable(mb['p'], color=color, markerstyle=25)
fig.save_to_root_file(fout, 'dNdeta_gen_comp_EtaLt15')


fig = Figure()
fig.plot.palette = 'colorblind'
fig.xtitle = r'#eta'
fig.ytitle = r'1/N_{evts} dN/d#eta (1/MB)'
for f in fs:
    c = f[0].MultEstimators.results_post.V0M.dNdeta_MB_ratio_summary
    tmp_fig = Figure()
    tmp_fig.import_plottables_from_canvas(c)
    tmp_fig._plottables = [tmp_fig._plottables[-1
    ]]
    fig.legend.title = tmp_fig._plottables[0]['legend_title'] + ' in V0M (1/MB)'
    tmp_fig._plottables[0]['legend_title'] = f[1]
    fig._plottables += tmp_fig._plottables
fig.save_to_root_file(fout, 'dNdeta_MB_ratio_gen_comp_V0M')

fig = Figure()
fig.plot.palette = 'colorblind'
fig.xtitle = r'#eta'
fig.ytitle = r'1/N_{evts} dN/d#eta (1/MB)'
for f in fs:
    c = f[0].MultEstimators.results_post.EtaLt15.dNdeta_MB_ratio_summary
    tmp_fig = Figure()
    tmp_fig.import_plottables_from_canvas(c)
    tmp_fig._plottables = [tmp_fig._plottables[-1
    ]]
    fig.legend.title = tmp_fig._plottables[0]['legend_title'] + ' in |#eta|#leq1.5 (1/MB)'
    tmp_fig._plottables[0]['legend_title'] = f[1]
    fig._plottables += tmp_fig._plottables
fig.save_to_root_file(fout, 'dNdeta_MB_ratio_gen_comp_EtaLt15')


##########
# P(Nch) #
##########

fig = Figure()
fig.legend.title = 'V0M, Py.8, Monash'
fig.plot.palette = 'colorblind'
fig.xtitle = r'N_{ch}^{|#eta|#leq0.5}'
fig.ytitle = r'P(N_{ch}^{|#eta|#leq0.5})'
fig.legend.position = 'tr'
fig.plot.xmin = 0
fig.plot.xmax = 150
fig.plot.ymin = 0.00001
fig.plot.ymax = 0.5
fig.plot.logy = True
fig.import_plottables_from_canvas(fs[1][0].MultEstimators.results_post.V0M.PNchEtaLt05_binned_in_NchEst)
fig.import_plottables_from_canvas(fs[2][0].MultEstimators.results_post.V0M.PNchEtaLt05_binned_in_NchEst)

for p_CR, p_noCR, color in zip(fig._plottables[:4], fig._plottables[4:], colors):
    # p_CR['markerstyle'] = 'square'
    p_noCR['markerstyle'] = 25
    p_noCR['color'] = p_CR['color'] = color
    p_noCR['legend_title'] = ''
fig.save_to_root_file(fout, 'PNchEst_binned_in_NchEtaLt05_Py8_Monash_noCR')

fig = Figure()
fig.legend.title = '0.1%-0.0% in V0M'
fig.plot.palette = 'colorblind'
fig.xtitle = r'N_{ch}^{|#eta|#leq0.5}'
fig.ytitle = r'P(N_{ch}^{|#eta|#leq0.5})'
fig.legend.position = 'tr'
fig.plot.xmin = 0
fig.plot.xmax = 150
fig.plot.ymin = 0.00001
fig.plot.ymax = 0.5
fig.plot.logy = True

for f in fs:
    c = f[0].MultEstimators.results_post.V0M.PNchEtaLt05_binned_in_NchEst
    tmp_fig = Figure()
    tmp_fig.import_plottables_from_canvas(c)
    tmp_fig._plottables = [tmp_fig._plottables[-1]]
    tmp_fig._plottables[0]['legend_title'] = f[1]
    fig._plottables += tmp_fig._plottables

fig.save_to_root_file(fout, 'PNchEst_binned_in_NchEtaLt05_gen_comp')


# p / pi in different estimators; two plots, one for V0M and one for EtaLt15
fig = Figure()
fig.legend.title = '|#eta|#leq1.5'
fig.plot.palette = 'colorblind'
fig.xtitle = r'N_{ch}^{|#eta|#leq0.5}'
fig.ytitle = r'p / #pi^{#pm}'
fig.plot.xmin = 0
fig.plot.xmax = 60
fig.plot.ymin = 0.045
fig.plot.ymax = 0.09
for f in fs:
    c = getattr(f[0].MultEstimators.results_post.pid_ratios_vs_refmult, "-2212_2212_div_-211_211")
    tmp_fig = Figure()
    tmp_fig.import_plottables_from_canvas(c)
    tmp_fig._plottables = [tmp_fig._plottables[2]]  # EtaLt15 should be the 3 in the list...
    print tmp_fig._plottables
    tmp_fig._plottables[0]['legend_title'] = f[1]
    fig._plottables += tmp_fig._plottables
fig.save_to_root_file(fout, 'p_div_pi_vs_mult_EtaLt15')

fig = Figure()
fig.legend.title = 'V0M'
fig.plot.palette = 'colorblind'
fig.xtitle = r'N_{ch}^{|#eta|#leq0.5}'
fig.ytitle = r'p / #pi^{#pm}'
fig.plot.xmin = 0
fig.plot.xmax = 60
fig.plot.ymin = 0.045
fig.plot.ymax = 0.09
for f in fs:
    c = getattr(f[0].MultEstimators.results_post.pid_ratios_vs_refmult, "-2212_2212_div_-211_211")
    tmp_fig = Figure()
    tmp_fig.import_plottables_from_canvas(c)
    tmp_fig._plottables = [tmp_fig._plottables[-1]]
    tmp_fig._plottables[0]['legend_title'] = f[1]
    fig._plottables += tmp_fig._plottables
fig.save_to_root_file(fout, 'p_div_pi_vs_mult_V0M')


# Xi / pi in different generators; once for V0M and once for EtaLt15
fig = Figure()
fig.legend.title = '|#eta|#leq1.5'
fig.plot.palette = 'colorblind'
fig.xtitle = r'N_{ch}^{|#eta|#leq0.5}'
fig.ytitle = r'#Xi / #pi^{#pm}'
fig.plot.xmin = 0
fig.plot.xmax = 60
fig.plot.ymin = 0.0005
fig.plot.ymax = 0.0018
for f in fs:
    c = getattr(f[0].MultEstimators.results_post.pid_ratios_vs_refmult, "3312_div_-211_211")
    tmp_fig = Figure()
    tmp_fig.import_plottables_from_canvas(c)
    tmp_fig._plottables = [tmp_fig._plottables[2]]
    tmp_fig._plottables[0]['legend_title'] = f[1]
    fig._plottables += tmp_fig._plottables
fig.save_to_root_file(fout, 'Xi_div_pi_vs_mult_EtaLt15')

fig = Figure()
fig.legend.title = 'V0M'
fig.plot.palette = 'colorblind'
fig.xtitle = r'N_{ch}^{|#eta|#leq0.5}'
fig.ytitle = r'#Xi / #pi^{#pm}'
fig.plot.xmin = 0
fig.plot.xmax = 60
fig.plot.ymin = 0.0005
fig.plot.ymax = 0.0018
for f in fs:
    c = getattr(f[0].MultEstimators.results_post.pid_ratios_vs_refmult, "3312_div_-211_211")
    tmp_fig = Figure()
    tmp_fig.import_plottables_from_canvas(c)
    tmp_fig._plottables = [tmp_fig._plottables[-1]]
    tmp_fig._plottables[0]['legend_title'] = f[1]
    fig._plottables += tmp_fig._plottables

fig.save_to_root_file(fout, 'Xi_div_pi_vs_mult_V0M')



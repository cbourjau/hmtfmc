"""
Given a folder structure with "estimator directories", find all
plot files and group them together. Also, name get a meaningfull name
for that group
"""

import subprocess
import sys
import os

from dump_plots_to_files import dump_plots_to_files
if len(sys.argv) != 3:
    print (
        "Usage:\n"
        'python make_summary_pres.py file.root "title of presentation"\n'
    )
    quit()

# create the plot files
dump_plots_to_files(sys.argv[1])

considered_ests = ['EtaLt05', 'EtaLt08', 'EtaLt15', 'Eta08_15', 'V0M', 'V0A', 'V0C', 'ZDC']
base_path = './' + sys.argv[1].split('.')[0] + '/'
process = subprocess.Popen("find . -type f | grep pdf",
                           shell=True,
                           stdout=subprocess.PIPE,
                           cwd=base_path
                           )
stdout_list = process.communicate()[0].split('\n')
stdout_list.remove('')
est0_substring = "/{}/".format(considered_ests[0])
est0_plots = [path for path in stdout_list if est0_substring in path]

# plots which summarize several ests and are thus not in an "est folder"
summary_plots = [path for path in stdout_list if est0_substring not in path]


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        yield l[i:i + n]


def make_frames(title, fig_paths):
    """Return tex code with up to 4 figures per frame"""
    tex = ''
    for chunk in chunks(fig_paths, 4):
        frame = (
            """\\begin{{frame}}[plain]
            \\frametitle{{{}}}
            \\begin{{columns}}
            """).format(title)
        ig_lines = ['', '', '', '']
        for i, path in enumerate(chunk):
            ig_lines[i] = r'\includegraphics[width=\textwidth]{{{}}}'.format(path)
        frame += (
            """\\begin{{column}}{{.45\\textwidth}}
            {}\\\\
            {}
            \\end{{column}}
            \\begin{{column}}{{.45\\textwidth}}
            {}\\\\
            {}
            \\end{{column}}
            \\end{{columns}}
            \\end{{frame}}
            """).format(*ig_lines)
        tex += frame
    return tex

groups = []
# cycle through the plots in the first estimator folder and find the equivalent plots
# in the other estimator folders
for plot_path in est0_plots:
    # get same plot of other estimators
    group = []
    # use the file name as the group name
    group_name = plot_path.split('/')[-1].split('.')[0].replace('_', ' ')
    for est in considered_ests:
        group.append(plot_path.replace(est0_substring, "/{}/".format(est)))
    # kick out plots which are not presents as files in the 'find' command
    group = [path for path in group if path in stdout_list]
    groups.append((group_name, group))

# go through all the plots which lie in the "top directory"
summary_group = []
for plot_path in summary_plots:
    # if the plots are in a folder use folder name as group name, else, use file name
    if plot_path.split('/')[-2] == 'results_post':  # the plot is in the top directory
        summary_group.append(plot_path)
groups.append(("Summary plots", summary_group))

# PID ratios vs refmult
pid_ratio_ref_group = []
for plot_path in summary_plots:
    if 'pid_ratios_vs_refmult' in plot_path:
        pid_ratio_ref_group.append(plot_path)
groups.append(("PID vs ref. multiplicity", pid_ratio_ref_group))

# PID ratios vs est
pid_ratio_est_group = []
for plot_path in summary_plots:
    if 'pid_ratios_vs_estmult' in plot_path:
        pid_ratio_est_group.append(plot_path)
groups.append(("PID vs est. multiplicity", pid_ratio_est_group))


tex_body = "\n".join([make_frames(name, paths) for (name, paths) in groups])

tex_doc = (
    r"""
    \documentclass[xcolor=dvipsnames]{{beamer}}

    \usepackage{{graphicx,subfigure,url, tikz}}

    % example themes
    \usetheme[nat,style=simple]{{Frederiksberg}}

    % put page numbers
    \setbeamertemplate{{footline}}[frame number]{{}}
    % remove navigation symbols
    \setbeamertemplate{{navigation symbols}}{{}}
    \author{{Christian Bourjau}}
    \institute[NBI, Copenhagen]{{Niels Bohr Institute, Copenhagen\\HMTF MC benchmark studies\\Supervisor: Michele Floris}}

    \title{{{}}}
    \begin{{document}}

    \frame[plain]{{\titlepage}}
    """.format(sys.argv[2])
)
tex_doc += tex_body
tex_doc += r'\end{document}'
tex_file_name = 'summary.tex'
with open(base_path + '/' + tex_file_name, 'w') as f:
    f.write(tex_doc)

# latex_out = subprocess.Popen('pdflatex -file-line-error -interaction=nonstopmode {}'.format(tex_file_name),
#                              cwd=base_path,
#                              shell=True)
latex_out = subprocess.check_output(
    ['pdflatex', '-file-line-error', '-interaction=nonstopmode', format(tex_file_name)],
    cwd=os.path.abspath(base_path))
latex_out = subprocess.check_output(
    ['pdflatex', '-file-line-error', '-interaction=nonstopmode', format(tex_file_name)],
    cwd=os.path.abspath(base_path))

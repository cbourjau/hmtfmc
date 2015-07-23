from rootpy import asrootpy
from rootpy.plotting import HistStack, Canvas, Legend, Pad, Hist1D, Graph
import ROOT


def create_dNdeta_stack(h2d, event_counter):
    """
    Create dN/deta stack for various multiplicity bins from given 2D histogram.
    xaxis: eta
    yaxis: multiplicity bins
    """
    stack = HistStack(name="dNdeta_stack")
    nbins = h2d.yaxis.GetNbins()
    for mult_bin in range(1, nbins):
        h2d.yaxis.set_range(mult_bin, mult_bin)
        stack.Add(asrootpy(h2d.projection_x()))
        # named colors of the ROOT TColor colorwheel are between 800 and 900, +1 to make them look better
        stack[-1].color = 800 + int(100.0/nbins)*(mult_bin-1) + 1
        stack[-1].name = str(mult_bin)
        stack[-1].title = (str(h2d.yaxis.get_bin_low_edge(mult_bin))
                           +' #leq N_{ch} #leq ' +
                           str(h2d.yaxis.get_bin_up_edge(mult_bin)))
        # scale by the number of events in this mult_bin
        # import ipdb; ipdb.set_trace()
        stack[-1].Scale(1.0/float(event_counter.Integral(mult_bin, mult_bin)))
    stack.Draw('nostack')
    stack.xaxis.SetTitle("#eta")
    stack.yaxis.SetTitle('1/N dN_{ch}/d#eta')
    return stack

def make_stack_of_mult_bins_for_pids(h3d, pids):
    """
    Make a histstack for the different mult bins. Takes the 3D hist and an itrerable of pids.
    The pids are added befor the stacks are created (eg. charged pions)
    Return a HistStack
    """
    h3d = asrootpy(h3d)
    stack = HistStack()
    pid_hists = []
    for pid in pids:
        h3d.zaxis.SetRange(pid+1, pid+1)  # pid + 1 = histogram bin :P 
        pid_h2d = asrootpy(h3d.Project3D("xy"))
        pid_h2d.name = pid_h2d.name[:-2] + "pid" + str(pid)
        pid_hists.append(pid_h2d)
    pid_sum_hist = sum(pid_hists)

    n_mult_classes = pid_h2d.yaxis.get_nbins()
    for ibin in range(1, n_mult_classes):
        pid_sum_hist.yaxis.set_range(ibin, ibin)
        tmp = asrootpy(pid_sum_hist.projection_x())
        if tmp.GetEntries() < 1:
            continue
        tmp.name = tmp.name + str(ibin)
        tmp.set_title(str(pid_sum_hist.yaxis.get_bin_low_edge(ibin))
                      + '#leq N_{ch}^{est} #leq '
                      + str(pid_sum_hist.yaxis.get_bin_up_edge(ibin)))
        stack.Add(tmp)
    return stack


def plot_histogram_stack(stack):
    """
    Plot a stack of histograms. The legend is generated from the titles of the histograms.
    """
    stack = asrootpy(stack)
    c = Canvas()
    pad1 = Pad(0., 0., .7, 1., name="plot")
    pad2 = Pad(.7, 0, 1., 1., name="legend")
    pad1.Draw()
    pad2.Draw()
    nesti = len(stack.GetHists())

    c.cd()
    pad1.cd(0)
    leg = Legend(entries=nesti, leftmargin=0, rightmargin=0, entrysep=0, entryheight=.04)
    leg.SetBorderSize(0)
    maximum = 0
    for mult_bin, h in enumerate(stack):
        h.color = 800 + int(100.0/len(stack))*(mult_bin) + 1
        h.markerstyle = 'circle'
        leg.AddEntry(h)
        if h.GetMaximum() > maximum:
            maximum = h.GetMaximum()
    stack.SetMaximum(maximum + maximum * 0.1)
    stack.Draw('nostack')
    stack.xaxis.SetTitle(stack.GetHists()[0].GetXaxis().GetTitle())
    stack.yaxis.SetTitle(stack.GetHists()[0].GetYaxis().GetTitle())
    c.cd()
    pad2.cd()
    leg.Draw()
    return c

def plot_list_of_plottables(l, title=''):
    """
    Plot the plottable objects, given as a list in a nice graph with legend. Legend is the plottables' titles.
    `title` is an optional title for the entire plot.
    """
    c = Canvas(width=1620, height=1000)
    pad1 = Pad(0., 0., .62, 1., name="plot")
    pad2 = Pad(.62, 0, 1., 1., name="legend")
    pad1.Draw()
    pad2.Draw()
    ntot = len(l)

    c.cd()
    pad1.cd(0)
    leg = Legend(entries=ntot, leftmargin=0, rightmargin=0, entrysep=0, entryheight=.04)
    leg.SetBorderSize(0)
    #maximum, minimum = 0, 0  # extreme values of all plottables
    isfirst = True           # switch to turn on drawingoption 'same'
    for i, p in enumerate(l):
        p.color = 800 + int(100.0/len(l))*(i) + 1
        p.markerstyle = 'circle'
        leg.AddEntry(p)
        # if p.GetMaximum() > maximum:
        #    maximum = p.GetMaximum()
        if isfirst:
            p.title = title
            p.Draw()
        else:
            p.Draw('same')
    # stack.SetMaximum(maximum + maximum * 0.1)
    # stack.xaxis.SetTitle(stack.GetHists()[0].GetXaxis().GetTitle())
    # stack.yaxis.SetTitle(stack.GetHists()[0].GetYaxis().GetTitle())
    c.cd()
    pad2.cd()
    leg.Draw()
    return c


def create_stack_pid_ratio_over_pt(h3d, pid1, pid2):
    """
    Create a hist stack, where hists are the ratio of particles of species 1 over 2 vs. pt, stack is
    binned in multiplicity bins.
    pidx must be a list, these particles are added together befor dividing (eg. pi charged)
    """
    h3d = asrootpy(h3d)   
    stack1 = make_stack_of_mult_bins_for_pids(h3d, pid1)
    stack2 = make_stack_of_mult_bins_for_pids(h3d, pid2)
    outstack = HistStack()
    for ibin, (hpid1, hpid2) in enumerate(zip(stack1.GetHists(), stack2.GetHists())):
        tmp = hpid1/hpid2
        tmp.title = hpid1.title
        outstack.Add(tmp)
    outstack.Draw('nostack')
    outstack.title = "pid{0} / pid{1}".format(pid1, pid2)
    outstack.xaxis.SetTitle("p_{T}")
    return outstack


def create_hist_pid_ratio_over_mult(h3d, pid1, pid2):
    """
    Create a histogram whith the ratio of particles of species 1 over 2 vs. multiplicity
    pidx must be a list, these particles are added together befor dividing (eg. pi charged)
    """
    h3d = asrootpy(h3d)
    edges = ([h3d.xaxis.get_bin_low_edge(i) for i in range(1, h3d.xaxis.get_nbins() +1)]
             + [h3d.xaxis.get_bin_up_edge(h3d.xaxis.get_nbins())])
    h_1 = Hist1D(edges)
    h_2 = Hist1D(edges)
    
    stack1 = make_stack_of_mult_bins_for_pids(h3d, pid1)
    stack2 = make_stack_of_mult_bins_for_pids(h3d, pid2)
    
    for ibin, (hpid1, hpid2) in enumerate(zip(stack1.GetHists(), stack2.GetHists())):
        err = ROOT.Double()
        cont = hpid1.IntegralAndError(1, hpid1.xaxis.GetNbins(), err)
        h_1.set_bin_content(ibin + 1, cont)
        h_1.set_bin_error(ibin + 1, err)

        cont = hpid2.IntegralAndError(1, hpid2.xaxis.GetNbins(), err)
        h_2.set_bin_content(ibin + 1, cont)
        h_2.set_bin_error(ibin + 1, err)
    ratio = h_1/h_2
    ratio.title = "pid{0} / pid{1}, {2}".format(pid1, pid2, h3d.title[30:])
    ratio.xaxis.title = "Multiplicity in estimator region"
    return ratio


def create_canonnical_avg_from_stacks(stacks):
    """
    stacks: list of HistStack
    Each stack has histograms for each mult_bin
    return: stack with the same hist binning in mult_bins 
    """
    n_estimators = len(stacks)
    avg_stack = stacks[0].Clone()
    avg_stack.name = "avg_stack"
    for s in stacks[1:]:
        if len(s) != len(avg_stack):
            raise ValueError("Given list of HistStacks do not contain equal number of hists")
        for h_sum, h_this_est in zip(avg_stack, s):
            h_sum.Add(h_this_est)
    # scale histograms in avg_stack by the number of estimators
    for h in avg_stack:
        h.Scale(1.0/n_estimators)
    return avg_stack

def divide_stacks(stack1, stack2):
    """Divide two given stacks of histograms if they contain the same number of hists"""
    if len(stack1) != len(stack2):
        raise ValueError('Given HistStacks do not contain same number of hists')
    outstack = HistStack()
    outstack.name = "{0}_div_by_{1}".format(stack1.name, stack2.name)
    for hpid1, hpid2 in zip(stack1.GetHists(), stack2.GetHists()):
        tmp = hpid1/hpid2
        tmp.title = hpid1.title
        outstack.Add(tmp)
    outstack.Draw('nostack')
    return outstack


def create_graph_pided_refest_vs_pidcount(h3d, corr_hist, pids):
    """Creates a graph with the count of the given pids (iterable) vs Nch of the REF mult.
    The mapping from Nch_est to Nch_ref is done using the correlation hist. 
    The ref est needs to be on the y axis of the correlation histogram.
    pids are given as strings.
    """
    profx = asrootpy(corr_hist.ProfileX())
    graphs =[]
    for pid in pids:
        # set pid
        graphs.append(Graph())
        pid_bin = h3d.zaxis.find_bin(pid)
        if pid_bin==0:
            raise ValueError("given pid does not exist")
        h3d.GetZaxis().SetRange(pid_bin, pid_bin)

        h2d = h3d.Project3D("yx", )
        count = asrootpy(h2d.ProjectionX())
        # start counting at 1 since the first bin is empty with INEL>0
        prof_bins = [b for b in profx.bins()]
        count_bins= [b for b in count.bins()]
        for i, (nch_ref_bin, counter_bin) in enumerate(zip(prof_bins[1:], count_bins[1:])): 
            if (counter_bin.value == 0.0):
                break
            graphs[-1].SetPoint(i, nch_ref_bin.value, counter_bin.value)
            xerr, yerr = nch_ref_bin.error/2,counter_bin.error/2
            graphs[-1].SetPointError(i,xerr, xerr,yerr, yerr)
    sgraphs = sum(graphs)
    sgraphs.title = ", ".join(pids)
    sgraphs.xaxis.title = "N_{{ch}}^{{{}}}".format(corr_hist.yaxis.title[7:])
    return sgraphs

from array import array
import sys, os, string, random

from rootpy import asrootpy
from rootpy.plotting import HistStack, Canvas, Legend, Pad, Hist1D, Graph

import ROOT


def gen_random_name():
    """Generate a random name for temp hists"""
    return ''.join(random.choice(string.ascii_letters + string.digits) for _ in range(25))


def create_dNdeta_stack(h2d, event_counter, with_mb=True):
    """
    Create dN/deta stack for various multiplicity bins from given 2D histogram. If `with_mb` is `True`,
    also add dNdeta for minimum bias to the stack.
    xaxis: eta
    yaxis: multiplicity bins
    """
    stack = HistStack(name="dNdeta_stack")
    nbins = h2d.yaxis.GetNbins()
    nch_step = 10
    last_mult_bin = False
    nevent_threshold = 5000
    for mult_bin in range(1, nbins, nch_step):
        mult_bin_upper = mult_bin + nch_step - 1
        h2d.yaxis.set_range(mult_bin, mult_bin_upper)
        h = asrootpy(h2d.projection_x())
        # Combine left over event classes?
        if h.Integral(mult_bin, mult_bin+nch_step) < nevent_threshold:
            last_mult_bin = True
            mult_bin_upper = nbins
            h2d.yaxis.set_range(mult_bin, mult_bin_upper)
            h = asrootpy(h2d.projection_x())

        # named colors of the ROOT TColor colorwheel are between 800 and 900, +1 to make them look better
        h.color = 800 + int(100.0/(nbins/nch_step)*(mult_bin/nch_step-1)) + 1
        h.name = str(mult_bin)
        h.title = (str(h2d.yaxis.get_bin_low_edge(mult_bin))
                   +' #leq N_{ch}^{est} #leq ' +
                   str(h2d.yaxis.get_bin_up_edge(mult_bin_upper)))
        h.yaxis.title = '1/N dN_{ch}^{est}/d#eta'
        # scale by the number of events in this mult_bin
        try:
            h.Scale(1.0/float(event_counter.Integral(mult_bin, mult_bin_upper)))
        except ZeroDivisionError:
            # No events in multiplicity range; break loop here
            break
        stack.Add(h)
        if last_mult_bin:
            break
    if with_mb:
        h2d.yaxis.set_range(0,-1)
        h = asrootpy(h2d.projection_x())
        h.title = 'MB'
        h.name  = 'mb_dNdeta'
        h.yaxis.title = '1/N dN_{ch}^{est}/d#eta'
        h.Scale(1.0/float(event_counter.Integral(0,-1)))
        stack.Add(h)
    stack.Draw('nostack')
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
        pid_bin = h3d.zaxis.find_bin(str(pid))
        if pid_bin==0:
            raise ValueError("given pid ({}) does not exist in histogram".format(pid))
        h3d.zaxis.SetRange(pid_bin, pid_bin)  # pid + 1 = histogram bin :P 
        pid_h2d = asrootpy(h3d.Project3D("xy"))
        pid_h2d.name = pid_h2d.name[:-2] + "pid" + str(pid)
        pid_hists.append(pid_h2d)
    pid_sum_hist = sum(pid_hists)

    n_mult_classes = pid_h2d.yaxis.get_nbins()
    for ibin in range(1, n_mult_classes):
        pid_sum_hist.yaxis.set_range(ibin, ibin)
        tmp = asrootpy(pid_sum_hist.projection_x())
        # skip if there are too few entries
        if tmp.GetEntries() < 1000:
            continue
        # rebin with davids binning
        bin_edges = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8, 0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,
                     2.2,2.4,2.6,2.8,3.0,3.3,3.6,3.9,4.2,4.6,5,5.4,5.9,
                     6.5,7,7.5,8,8.5,9.2,10,11,12,13.5,15,17,20]
        tmp_rebinned = asrootpy(tmp.rebin(len(bin_edges)-1, tmp.name+"rebinned", array("d", bin_edges)))
        tmp_rebinned.name = tmp_rebinned.name + str(ibin)
        tmp_rebinned.set_title(str(pid_sum_hist.yaxis.get_bin_low_edge(ibin))
                               + '#leq N_{ch}^{est} #leq '
                               + str(pid_sum_hist.yaxis.get_bin_up_edge(ibin)))
        stack.Add(tmp_rebinned)
    return stack


def plot_histogram_stack(stack):
    """
    Plot a stack of histograms. The legend is generated from the titles of the histograms.
    """
    stack = asrootpy(stack)
    c = Canvas()
    pad1 = Pad(0., 0., .8, 1., name="plot")
    pad1.SetRightMargin(.015)
    pad2 = Pad(.8, 0, 1., 1., name="legend")
    pad2.SetLeftMargin(.0)
    pad1.Draw()
    pad2.Draw()

    nesti = len(stack.GetHists())

    c.cd()
    pad1.cd(0)
    leg = Legend(entries=nesti, leftmargin=0, rightmargin=0, entrysep=0, entryheight=.04, textsize=.1)
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
    c = Canvas()
    pad1 = Pad(0., 0., .8, 1., name="plot")
    pad2 = Pad(.8, 0, 1., 1., name="legend")
    pad2.SetLeftMargin(.03)
    pad1.Draw()
    pad2.Draw()
    ntot = len(l)

    c.cd()
    pad1.cd(0)
    leg = Legend(entries=ntot, leftmargin=0, rightmargin=0, entrysep=0, entryheight=.04, textsize=.15)
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
            p.Draw('ALP')
            isfirst = False
        else:
            p.Draw('LPsame')

    c.cd()
    pad2.cd()
    leg.Draw()
    return c


def create_stack_pid_ratio_over_pt(mult_pt_dir, pids1, pids2, nch_max):
    """
    Create a hist stack, where hists are the ratio of particles of species 1 over 2 vs. pt, stack is
    binned in multiplicity bins.
    pidx must be a list, these particles are added together befor dividing (eg. pi charged)
    """
    mult_bin_size = 10
    nmult_bins = mult_pt_dir.FindObject(pids1[0]).GetNbinsX()

    def make_list_of_pt_plots_binned_in_multclasses(pids, nch_max):
        pt_hists_in_multclasses_in_pids = []
        for ipid, pid in enumerate(pids):
            h2d = asrootpy(mult_pt_dir.FindObject(pid))
            pt_hists_in_multclasses = []
            is_last_bin = False
            for mult_bin_lower in range(0, nmult_bins, 10):
                mult_bin_upper = mult_bin_lower + mult_bin_size
                if mult_bin_upper > nch_max:
                    mult_bin_upper = h2d.xaxis.GetNbins()
                    is_last_bin = True
                h2d.xaxis.SetRange(mult_bin_lower, mult_bin_lower + mult_bin_size)
                pt_hists_in_multclasses.append(asrootpy(h2d.ProjectionY(gen_random_name())))
                pt_hists_in_multclasses[-1].title = "{} #leq N_{{ch}}^{{est}} #leq {}".\
                                                    format(mult_bin_lower, mult_bin_upper)
                if is_last_bin:
                    pt_hists_in_multclasses[-1].title = "{} #leq N_{{ch}}^{{est}}".\
                                                        format(mult_bin_lower)
                    break
            pt_hists_in_multclasses_in_pids.append(pt_hists_in_multclasses)
        # sum over pids -> list of length multclass
        # for that, I make convert the npid lists with length nmultclass to
        # nmultclass tuples of lenghtn npid
        return map(sum, zip(*pt_hists_in_multclasses_in_pids))
    pt_hists_pids1 = make_list_of_pt_plots_binned_in_multclasses(pids1, nch_max)
    pt_hists_pids2 = make_list_of_pt_plots_binned_in_multclasses(pids2, nch_max)

    # divide the sum of the two pids lists; results in list of ratios for each mutlclass
    outstack = HistStack()    
    for h1,h2 in zip(pt_hists_pids1, pt_hists_pids2):
        #tmp = hpid1/hpid2
        #tmp.title = hpid1.title
        #outstack.Add(tmp)
        if h1.Integral()==0.0 or h2.Integral()==0.0:
            # break if one of the histograms is empty
            break
        try:
            outstack.Add(h1/h2)
        except ZeroDivisionError:
            # break at first zero division
            break

    outstack.Draw('nostack')
    outstack.title = "pid{0} / pid{1}".format(pids1, pids2)
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


def remap_x_values(hist, map_hist):
    """Map the x values of hist to the y values of map_hist.
    In order to do so, it is necessary that the x values of hist are also present as x-values in map_hist.

    Parameters
    ----------
    hist : Hist1D
    map_hist : Hist1D

    Returns
    -------
    Graph
            Graph of the remapped hist. Errors are ??? TODO
    """
    hist = asrootpy(hist)
    map_hist = asrootpy(map_hist)
    rt_graph = Graph()
    for i, (nch_ref_bin, counter_bin) in enumerate(zip(map_hist.bins(), hist.bins())):
        rt_graph.SetPoint(i, nch_ref_bin.value, counter_bin.value)
        xerr, yerr = nch_ref_bin.error/2.0, counter_bin.error/2.0
        rt_graph.SetPointError(i,xerr, xerr,yerr, yerr)
    return rt_graph

def remove_zero_value_points(g):
    # Remove the points backwards, since the index would change if we do it forwards
    # The first point has index 0!
    points_to_remove = []
    for i, (x, y) in enumerate(g):
        if not y > 0.0:
            points_to_remove.append(i)
    for p in points_to_remove[::-1]:
        g.RemovePoint(p)


def remove_points_with_equal_x(g):
    """Remove all points which are on already occupied x values. Ie. The first point is kept, all later ones removed"""
    points_to_remove = []
    seen_x = []
    for i, (x, y) in enumerate(g):
        if x in seen_x:
            points_to_remove.append(i)
        else:
            seen_x.append(x)
    for p in points_to_remove[::-1]:
        g.RemovePoint(p)

        
def remove_points_with_x_err_gt_1NchRef(g):
    npoints = g.GetN()
    points_to_remove = []
    for idx in xrange(0,npoints):
        if g.GetErrorX(idx) > 1:
            points_to_remove.append(idx)
    for p in points_to_remove[::-1]:
        g.RemovePoint(p)
    
            
        
def remove_non_mutual_points(g1, g2):
    """Remove all points with do no have a corresponding point at the same x-value in the other hist"""
    points_to_remove1 = []
    points_to_remove2 = []
    xs1 = [p[0] for p in g1]
    xs2 = [p[0] for p in g2]
    for i, x in enumerate(xs1):
        if x not in xs2:
            points_to_remove1.append(i)
    for i, x in enumerate(xs2):
        if x not in xs1:
            points_to_remove2.append(i)
    for p in points_to_remove1[::-1]:
        g1.RemovePoint(p)
    for p in points_to_remove2[::-1]:
        g2.RemovePoint(p)


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
        pid_bin = h3d.zaxis.find_bin(pid)
        if pid_bin==0:
            raise ValueError("given pid ({}) does not exist in histogram".format(pid))

        h3d.GetZaxis().SetRange(pid_bin, pid_bin)
        h2d = h3d.Project3D("yx", )
        count = asrootpy(h2d.ProjectionX())
        
        graphs.append(remap_x_values(count, profx))
    sgraphs = sum(graphs)
    sgraphs.title = ", ".join(pids)
    sgraphs.xaxis.title = "N_{{ch}}^{{{}}}".format(corr_hist.yaxis.title[7:])
    return sgraphs

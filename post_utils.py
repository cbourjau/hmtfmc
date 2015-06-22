from rootpy.io import root_open
from rootpy import asrootpy
from rootpy.plotting import HistStack, Canvas, Legend, Pad, Hist1D
import ROOT


def create_dNdeta_stack(h2d, event_counter):
    """
    Create dN/deta stack for various multiplicity bins from given 2D histogram.
    xaxis: eta
    yaxis: multiplicity bins
    """
    stack = HistStack()
    nbins = h2d.yaxis.GetNbins()
    for mult_bin in range(1, nbins):
        h2d.yaxis.set_range(mult_bin, mult_bin)
        stack.Add(asrootpy(h2d.projection_x()))
        # named colors of the ROOT TColor colorwheel are between 800 and 900, +1 to make them look better
        stack[-1].color = 800 + int(100.0/nbins)*(mult_bin-1) + 1
        stack[-1].name = str(mult_bin)
        stack[-1].title = (str(h2d.yaxis.get_bin_low_edge(mult_bin))
                           +'$ \le N_{ch} < $' +
                           str(h2d.yaxis.get_bin_up_edge(mult_bin)))
        # scale by the number of events in this mult_bin
        stack[-1].Scale(1.0/float(event_counter.Integral(mult_bin, mult_bin)))
    stack.Draw('nostack')
    stack.xaxis.SetTitle("$\eta$")
    stack.yaxis.SetTitle('$1/N dN_{ch}/d\eta$')
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
                      + '$ \le N_{ch} \le $'
                      + str(pid_sum_hist.yaxis.get_bin_up_edge(ibin)))
        stack.Add(tmp)
    return stack


def plot_stack_of_estimators(stack):
    """
    Plot a stack of histograms where each represents a different estimator bin.
    Returns a canvas
    """
    stack = asrootpy(stack)
    c = Canvas()
    pad1 = Pad(0., 0., .7, 1.,)
    pad2 = Pad(.7, 0, 1., 1.)
    pad1.Draw()
    pad2.Draw()
    nesti = len(stack.GetHists())

    c.cd()
    pad1.cd(0)
    leg = Legend(entries=nesti, leftmargin=0, rightmargin=0, header="Estimator bins", entrysep=0, entryheight=.04)
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
    c.Update()
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
    outstack.title = "pid{0} / pid{1}, estimator: |#eta | < 0.5".format(pid1, pid2)
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

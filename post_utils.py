import string
import random

from rootpy import asrootpy
from rootpy.plotting import HistStack, Canvas, Legend, Pad, Graph
from rootpy.plotting.utils import get_limits

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
        if h.Integral(mult_bin, mult_bin + nch_step) < nevent_threshold:
            last_mult_bin = True
            mult_bin_upper = nbins
            h2d.yaxis.set_range(mult_bin, mult_bin_upper)
            h = asrootpy(h2d.projection_x())

        # named colors of the ROOT TColor colorwheel are between 800 and 900, +1 to make them look better
        h.color = 800 + int(100.0 / (nbins / nch_step) * (mult_bin / nch_step - 1)) + 1
        h.name = str(mult_bin)
        h.title = (str(h2d.yaxis.get_bin_low_edge(mult_bin))
                   + ' #leq N_{ch}^{est} #leq ' +
                   str(h2d.yaxis.get_bin_up_edge(mult_bin_upper)))
        h.yaxis.title = '1/N dN_{ch}^{est}/d#eta'
        # scale by the number of events in this mult_bin
        try:
            h.Scale(1.0 / float(event_counter.Integral(mult_bin, mult_bin_upper)))
        except ZeroDivisionError:
            # No events in multiplicity range; break loop here
            break
        stack.Add(h)
        if last_mult_bin:
            break
    if with_mb:
        h2d.yaxis.set_range(0, -1)
        h = asrootpy(h2d.projection_x())
        h.title = 'MB'
        h.name = 'mb_dNdeta'
        h.yaxis.title = '1/N dN_{ch}^{est}/d#eta'
        h.Scale(1.0 / float(event_counter.Integral(0, -1)))
        stack.Add(h)
    stack.Draw('nostack')
    return stack


def remap_x_values(hist, corr_hist):
    """
    Map the x values of hist to the y values of map_hist.
    In order to do so, it is necessary that the x values of hist are also present as x-values in map_hist.

    Parameters
    ----------
    hist : Hist1D
    corr_hist : Hist2D
                Correlations between the quantity on hist's x-axis (also corr_hist's xaxis) and the new
                quantity to plot agains (on corr_hist's y-axis.

    Returns
    -------
    Graph
            Graph of the remapped hist. Errors are ??? TODO
    """
    hist = asrootpy(hist)
    corr_hist = asrootpy(corr_hist)
    profx = asrootpy(corr_hist.ProfileX(gen_random_name()))

    rt_graph = Graph()
    for i, (nch_ref_bin, counter_bin) in enumerate(zip(profx.bins(), hist.bins())):
        rt_graph.SetPoint(i, nch_ref_bin.value, counter_bin.value)
        xerr, yerr = nch_ref_bin.error / 2.0, counter_bin.error / 2.0
        rt_graph.SetPointError(i, xerr, xerr, yerr, yerr)
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
    for idx in xrange(0, npoints):
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
    graphs = []
    for pid in pids:
        # set pid
        pid_bin = h3d.zaxis.find_bin(pid)
        if pid_bin == 0:
            raise ValueError("given pid ({}) does not exist in histogram".format(pid))

        h3d.GetZaxis().SetRange(pid_bin, pid_bin)
        h2d = h3d.Project3D("yx", )
        count = asrootpy(h2d.ProjectionX())
        graphs.append(remap_x_values(count, profx))
    sgraphs = sum(graphs)
    sgraphs.title = ", ".join(pids)
    sgraphs.xaxis.title = "N_{{ch}}^{{{}}}".format(corr_hist.yaxis.title[7:])
    return sgraphs

import unittest

from rootpy.plotting import Hist1D, Graph
from ROOT import TCanvas, TLegend
from roofi import Figure

import ROOT
ROOT.gROOT.SetBatch(False)


class Test_Figure(unittest.TestCase):
    def test_initiate_figure(self):
        f = Figure()
        self.assertIsInstance(f, Figure)

    def test_add_plottable(self):
        f = Figure()
        h = Hist1D(10, 0, 10)
        f.add_plottable(h)
        # check if a deep copy was performed
        self.assertNotEqual(h, f._plottables[0])

        # add histogram with legend
        self.assertEqual(len(f._legend_labels), 0)
        f.add_plottable(h, legend_title="cool hist")
        self.assertEqual(len(f._legend_labels), 1)

        # no old plottables if I make a new one:
        f = Figure()
        self.assertEqual(len(f._plottables), 0)


class Test_draw_to_canvas(unittest.TestCase):
    def test_draw_without_plottables(self):
        f = Figure()
        self.assertRaises(IndexError, f.draw_to_canvas)

    def test_draw_to_canvas_no_legend(self):
        f = Figure()
        h = Hist1D(10, 0, 10)
        h.Fill(5)
        f.add_plottable(h)
        c = f.draw_to_canvas()
        self.assertIsInstance(c, TCanvas)
        c.SaveAs("test_no_legend.pdf")

    def test_draw_to_canvas_legend(self):
        f = Figure()
        h = Hist1D(10, 0, 10)
        h.Fill(5)
        f.add_plottable(h, legend_title="some title")
        c = f.draw_to_canvas()
        c.SaveAs("test_legend.pdf")

    def test_draw_to_canvas_several_hists(self):
        f = Figure()
        h1 = Hist1D(10, 0, 10)
        h1.Fill(5)
        h2 = Hist1D(10, 0, 10)
        h2.Fill(2)
        f.add_plottable(h1, legend_title="hist 1")
        f.add_plottable(h2, legend_title="hist 1")
        c = f.draw_to_canvas()
        c.SaveAs("test_several.pdf")

    def test_draw_to_canvas_hists_and_graphs(self):
        f = Figure()
        h1 = Hist1D(10, 0, 10)
        h1.Fill(5)
        h2 = Hist1D(10, 0, 10)
        h2.Fill(2)
        gr = Graph()
        gr.SetPoint(0, 0, 0)
        gr.SetPoint(1, 1, 1)
        gr.SetPoint(2, 7, 2)
        f.add_plottable(h1, "hist")
        f.add_plottable(gr, "graph")
        c = f.draw_to_canvas()
        c.SaveAs("test_h_and_graph.pdf")


class Test_plot_options(unittest.TestCase):
    def test_log_scales(self):
        f = Figure()
        h1 = Hist1D(10, 0, 10)
        h1.Fill(5)
        f.add_plottable(h1)
        f.plot.logx = True
        f.plot.logy = True
        c = f.draw_to_canvas()
        self.assertEqual(c.FindObject("plot").GetLogy(), 1)


class Test_legend_options(unittest.TestCase):
    def setUp(self):
        self.f = Figure()
        h1 = Hist1D(10, 0, 10)
        h1.Fill(5)
        h2 = Hist1D(10, 0, 10)
        h2.Fill(2)
        self.f.add_plottable(h1, legend_title="hist 1")
        self.f.add_plottable(h2, legend_title="hist 1")

    def test_legend_postion(self):
        self.f.legend.position = 'tl'
        c = self.f.draw_to_canvas()
        leg = next(obj for obj in c.FindObject("plot").GetListOfPrimitives() if isinstance(obj, TLegend))
        self.assertGreater(leg.y1, .5)

        self.f.legend.position = 'bl'
        c = self.f.draw_to_canvas()
        leg = next(obj for obj in c.FindObject("plot").GetListOfPrimitives() if isinstance(obj, TLegend))
        self.assertLess(leg.y1, .5)

        self.f.legend.position = 'tr'
        c = self.f.draw_to_canvas()
        leg = next(obj for obj in c.FindObject("plot").GetListOfPrimitives() if isinstance(obj, TLegend))
        self.assertGreater(leg.x1, .5)
        self.assertGreater(leg.y1, .5)

        self.f.legend.position = 'br'
        c = self.f.draw_to_canvas()
        leg = next(obj for obj in c.FindObject("plot").GetListOfPrimitives() if isinstance(obj, TLegend))
        self.assertGreater(leg.x1, .5)
        self.assertLess(leg.y1, .5)
        c.SaveAs("leg_br.pdf")

        self.f.legend.position = 'seperate'
        c = self.f.draw_to_canvas()
        pad = c.FindObject("legend")
        self.assertIsInstance(pad, ROOT.TPad)
        import ipdb; ipdb.set_trace()


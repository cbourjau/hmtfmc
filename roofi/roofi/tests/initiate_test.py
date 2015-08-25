import unittest

from rootpy.plotting import Hist1D, Graph
from rootpy.interactive import wait
from rootpy.io import File
from ROOT import TCanvas, TLegend, TFile, TDirectoryFile
from roofi import Figure

import ROOT


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

    def test_reset_figure(self):
        f = Figure()
        h = Hist1D(10, 0, 10)
        f.add_plottable(h)
        f.delete_plottables()
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

    def test_axis_labels_dont_overlap(self):
        ROOT.gROOT.SetBatch(False)
        f = Figure()
        f.xtitle = 'N_{ch}#times#eta / #phi'
        f.ytitle = '1/N_{ch}^{supscr}#pi^{#pm}/#Xi'
        h1 = Hist1D(10, 0, 10,)
        h1.Fill(5)
        f.add_plottable(h1)
        c = f.draw_to_canvas()
        wait()

    def test_draw_half_width(self):
        ROOT.gROOT.SetBatch(False)
        f = Figure()
        f.style = 'presentation_half'
        f.xtitle = 'N_{ch}#times#eta / #phi'
        f.ytitle = '1/N_{ch}^{supscr}#pi^{#pm}/#Xi'
        h1 = Hist1D(10, 0, 10,)
        h1.Fill(5)
        f.add_plottable(h1)
        c = f.draw_to_canvas()

    def test_SetRangeUser(self):
        f = Figure()
        h1 = Hist1D(10, 0, 10)
        h1.Fill(5)
        f.add_plottable(h1)
        f.plot.xmin = -5
        f.plot.xmax = 5
        f.plot.ymin = -5
        f.plot.ymax = 5
        c = f.draw_to_canvas()
        plot_pad = c.FindObject("plot")
        self.assertEqual(plot_pad.GetListOfPrimitives()[-1].GetXaxis().GetXmin(), -5)
        self.assertEqual(plot_pad.GetListOfPrimitives()[-1].GetXaxis().GetXmax(), 5)
        self.assertEqual(plot_pad.GetListOfPrimitives()[-1].GetMinimum(), -5)
        self.assertEqual(plot_pad.GetListOfPrimitives()[-1].GetMaximum(), 5)


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


class Test_write_to_root_file(unittest.TestCase):
    def setUp(self):
        self.fig = Figure()
        h1 = Hist1D(10, 0, 10)
        h1.Fill(5)
        self.fig.add_plottable(h1, legend_title="hist 1")

    @unittest.skip(True)
    def test_write_to_TFile(self):
        f = TFile("test.root", "recreate")
        f = self.fig.save_to_root_file(f, 'myname')
        # Close might raise an error here!
        f.Close()

    def test_write_to_root_file(self):
        f = File("test.root", "recreate")
        f = self.fig.save_to_root_file(f, 'myname')
        f.close()
        f = TFile("test.root", "read")
        self.assertIsInstance(f.Get("myname"), TCanvas)

    def test_write_to_root_file_with_path(self):
        f = File("test.root", "recreate")
        f = self.fig.save_to_root_file(f, 'myname', path='folder/')
        f.Close()
        f = TFile("test.root", "read")
        self.assertIsInstance(f.Get("folder"), TDirectoryFile)
        self.assertIsInstance(f.Get("folder").Get("myname"), TCanvas)

import string
import random

from rootpy import asrootpy
from rootpy.plotting import Legend, Canvas, Pad
from rootpy.plotting.utils import get_limits

import ROOT


def gen_random_name():
    """Generate a random name for temp hists"""
    return ''.join(random.choice(string.ascii_letters + string.digits) for _ in range(25))


def get_color_generator(palett='set2'):
    """Returns a generator for n colors"""
    # generated with sns.palplot(sns.color_palette("colorblind", 10))
    if 'colorblind':
        return iter([(0.0, 0.4470588235294118, 0.6980392156862745),
                     (0.0, 0.6196078431372549, 0.45098039215686275),
                     (0.8352941176470589, 0.3686274509803922, 0.0),
                     (0.8, 0.4745098039215686, 0.6549019607843137),
                     (0.9411764705882353, 0.8941176470588236, 0.25882352941176473),
                     (0.33725490196078434, 0.7058823529411765, 0.9137254901960784),
                     (0.0, 0.4470588235294118, 0.6980392156862745),
                     (0.0, 0.6196078431372549, 0.45098039215686275),
                     (0.8352941176470589, 0.3686274509803922, 0.0),
                     (0.8, 0.4745098039215686, 0.6549019607843137)])
    # return iter(colorblind_colors)
    if 'set2':
        return iter([(0.40000000596046448, 0.7607843279838562, 0.64705884456634521),
                     (0.98131487965583808, 0.55538641635109398, 0.38740485135246722),
                     (0.55432528607985565, 0.62711267120697922, 0.79595541393055635),
                     (0.90311419262605563, 0.54185316071790801, 0.76495195557089413),
                     (0.65371782148585622, 0.84708959004458262, 0.32827375098770734),
                     (0.9986312957370983, 0.85096502233954041, 0.18488274134841617),
                     (0.89573241682613591, 0.76784315109252932, 0.58182240093455595),
                     (0.70196080207824707, 0.70196080207824707, 0.70196080207824707),
                     (0.40000000596046448, 0.7607843279838562, 0.64705884456634521),
                     (0.98131487965583808, 0.55538641635109398, 0.38740485135246722)])
    else:
        raise ValueError("Unknonw palette")


class Figure(object):
    def __init__(self):
        # User settable parameters:
        self.title = ''
        self.xtitle = ''
        self.ytitle = ''
        self.logx, self.logy = False, False
        # ncolors = 10
        # Private:
        self._plottables = []
        self._legend_labels = []
        self.style = 'presentation'
        self.plot = self.Plot()
        self.legend = self.Legend()

    class Plot(object):
        logx = False
        logy = False
        # xmin, xmax, ymin, ymax = None, None, None, None

    class Legend(object):
        title = ''
        position = 'tl'

    def _create_legend(self):
        nentries = len(self._legend_labels)
        leg = Legend(nentries, leftmargin=0, rightmargin=0, entrysep=0.02, textsize=14, textfont=63, margin=0.1, )
        leg.SetBorderSize(0)  # no box
        leg.SetFillStyle(0)   # transparent background of legend TPave(!)
        return leg

    def _theme_plottable(self, obj):
        axes = obj.GetXaxis(), obj.GetYaxis()
        for axis in axes:
            axis.SetLabelSize(14)
            axis.SetLabelFont(63)

            axis.SetTitleFont(63)
            axis.SetTitleSize(14)
            axis.SetTitleFont(63)
            axis.SetTitleSize(14)

    def add_plottable(self, obj, legend_title=''):
        """
        Add a plottable objet to this figure. This function performs a
        copy of the passed object and assigns it a random name. Once
        commited, these should not be touched any more by the user!!!

        Parameters
        ----------
        obj : Hist1D, Graph
              A root plottable object

        """
        self._plottables.append(asrootpy(obj.Clone(gen_random_name())))
        if legend_title:
            # check if we already have a legend:
            self._legend_labels.append((self._plottables[-1], legend_title))

    def draw_to_canvas(self):
        """
        Draw this figure to a canvas, which is then returned.
        """
        if len(self._plottables) == 0:
            raise IndexError("No plottables defined")
        c = Canvas(width=340, height=300, size_includes_decorations=True)
        if self.legend.position == 'seperate':
            legend_width = .2
            pad_legend = Pad(1 - legend_width, 0, 1., 1., name="legend")
            pad_legend.SetLeftMargin(0.0)
            pad_legend.SetFillStyle(0)  # make this pad transparent
            pad_legend.Draw()
        else:
            legend_width = 0
        pad_plot = Pad(0., 0., 1 - legend_width, 1., name="plot", )
        pad_plot.SetRightMargin(.04)
        pad_plot.Draw()
        pad_plot.cd()

        xmin, xmax, ymin, ymax = get_limits(self._plottables, logx=self.plot.logx, logy=self.plot.logy)

        colors = get_color_generator()

        is_first = True
        for i, obj in enumerate(self._plottables):
            obj.markerstyle = 'circle'
            color = next(colors)
            obj.color = color

            obj.SetMinimum(ymin)
            obj.SetMaximum(ymax)
            obj.GetYaxis().SetLimits(ymin, ymax)
            obj.GetYaxis().SetRangeUser(ymin, ymax)

            obj.GetXaxis().SetTitle(self.xtitle)
            obj.GetYaxis().SetTitle(self.ytitle)

            # Set the title to the given title:
            obj.title = self.title

            self._theme_plottable(obj)

            if isinstance(obj, ROOT.TGraph):
                if is_first:
                    drawoption = 'APL'  # root think axis are over rated for graphs...
                else:
                    drawoption = 'PLsame'
            elif isinstance(obj, ROOT.TH1):
                obj.SetStats(0)
                if is_first:
                    drawoption = ''
                else:
                    drawoption = 'same'
            else:
                raise TypeError("Un-plottable type given.")
            is_first = False
            obj.Draw(drawoption)
        pad_plot.SetTicks()
        pad_plot.SetLogx(self.logx)
        pad_plot.SetLogy(self.logy)

        if len(self._legend_labels) > 0:
            leg = self._create_legend()
            longest_label = 0
            for obj, lable in self._legend_labels:
                leg.AddEntry(obj, lable)
                if len(lable) > longest_label:
                    longest_label = len(lable)

            # Set the legend position
            # vertical:
            if self.legend.position.startswith('t'):
                leg_hight = leg.y2 - leg.y1
                leg.y2 = 1 - pad_plot.GetTopMargin()
                leg.y1 = leg.y2 - leg_hight
            elif self.legend.position.startswith('b'):
                leg_hight = leg.y2 - leg.y1
                leg.y1 = pad_plot.GetBottomMargin()
                leg.y2 = leg.y1 + leg_hight
            # horizontal:
            if self.legend.position[1:].startswith('l'):
                pass
            elif self.legend.position[1:].startswith('r'):
                leg_width = 0.3
                leg.x2 = 1 - pad_plot.GetRightMargin()
                leg.x1 = leg.x2 - leg_width
            if self.legend.position == 'seperate':
                with pad_legend:
                    leg.Draw()
            else:
                leg.Draw()
        if self.plot.logx:
            pad_plot.SetLogx(True)
        if self.plot.logy:
            pad_plot.SetLogy(True)
        return c

    def delete_plottables(self):
        """
        Delete all plottables in this figure so that it can be filled with
        new ones while keeping the lables.
        """
        self._plottables = []
        self._legend_labels = []

    def save_to_root_file(self, in_f, name, path=''):
        """
        Save the current figure to the given root file under the given path
        Parameters
        ----------
        f : TFile
            Root file object open in writable mode
        name : str
            Name for the canvas in the root file
        path : str
            The path where the figure should be saved within the root file
        Returns
        -------
        TFile :
            The file where the object was written to
        """
        f = asrootpy(in_f)
        c = self.draw_to_canvas()
        c.name = name
        try:
            f.mkdir(path, recurse=True)
        except ValueError:
            pass
        f.cd(path)
        c.Write()
        return f
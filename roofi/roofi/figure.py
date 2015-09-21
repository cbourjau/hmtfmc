import string
import random
import logging
import os

from rootpy import asrootpy, log
from rootpy.plotting import Legend, Canvas, Pad
from rootpy.plotting.utils import get_limits

import ROOT

# from external import husl


class Styles(object):
    # Define names of plot layouts:
    class _Default_Style(object):
        pt_per_cm = 28.4527625
        titlefont = 43
        labelfont = 43
        markerSizepx = 4  # number of pixels of the marker

    class Presentation_full(_Default_Style):
        axisTitleSize = 14
        axisLabelSize = 14
        legendSize = 14
        canvasWidth = 340
        canvasHeight = 300
        plot_margins = (.13, .05, .13, .04)   # left, right, bottom, top
        plot_ytitle_offset = 1.15  # factor of the normal offset :P, may lay outside of the canvas

    class Presentation_half(_Default_Style):
        axisTitleSize = 10
        axisLabelSize = 10
        legendSize = 10
        canvasWidth = 170
        canvasHeight = 150
        plot_margins = (.3, .08, .2, .04)
        plot_ytitle_offset = 1

    class Public_full(_Default_Style):
        axisTitleSize = 10
        axisLabelSize = 8
        legendSize = 8
        canvasWidth = 340
        canvasHeight = 300
        plot_margins = (.13, .05, .13, .04)
        plot_ytitle_offset = 1.15


logging.basicConfig(level=logging.DEBUG)
log = log["/roofi"]


def gen_random_name():
    """Generate a random name for temp hists"""
    return ''.join(random.choice(string.ascii_letters + string.digits) for _ in range(25))


def get_color_generator(palette='root', ncolors=10):
    """
    Returns a generator for n colors.
    Parameters
    ----------
    ncolors : int
              number of colors this palette should have, it might be ignored by some palettes!
    Returns
    -------
    generator :
               colors which can be digested by _rootpy_
    """
    # generated with sns.palplot(sns.color_palette("colorblind", 10))
    if palette == 'colorblind':
        colors = ([(0.0, 0.4470588235294118, 0.6980392156862745),
                   (0.0, 0.6196078431372549, 0.45098039215686275),
                   (0.8352941176470589, 0.3686274509803922, 0.0),
                   (0.8, 0.4745098039215686, 0.6549019607843137),
                   (0.9411764705882353, 0.8941176470588236, 0.25882352941176473),
                   (0.33725490196078434, 0.7058823529411765, 0.9137254901960784)])
    if palette == 'set2':
        colors = ([(0.40000000596046448, 0.7607843279838562, 0.64705884456634521),
                   (0.98131487965583808, 0.55538641635109398, 0.38740485135246722),
                   (0.55432528607985565, 0.62711267120697922, 0.79595541393055635),
                   (0.90311419262605563, 0.54185316071790801, 0.76495195557089413),
                   (0.65371782148585622, 0.84708959004458262, 0.32827375098770734),
                   (0.9986312957370983, 0.85096502233954041, 0.18488274134841617),
                   (0.89573241682613591, 0.76784315109252932, 0.58182240093455595),
                   (0.70196080207824707, 0.70196080207824707, 0.70196080207824707)])
    if palette == 'husl':
        colors = [(0.9677975592919913, 0.44127456009157356, 0.5358103155058701),
                  (0.8616090647292522, 0.536495730113334, 0.19548899031476086),
                  (0.6804189127793346, 0.6151497514677574, 0.19405452111445337),
                  (0.46810256823426105, 0.6699492535792404, 0.1928958739904499),
                  (0.20125317221201128, 0.6907920815379025, 0.47966761189275336),
                  (0.21044753832183283, 0.6773105080456748, 0.6433941168468681),
                  (0.2197995660828324, 0.6625157876850336, 0.7732093159317209),
                  (0.433280341176423, 0.6065273407962815, 0.9585467098271748),
                  (0.8004936186423958, 0.47703363533737203, 0.9579547196007522),
                  (0.962272393509669, 0.3976451968965351, 0.8008274363432775)]
    if palette == 'root':
        # named colors of the ROOT TColor colorwheel are between 800 and 900, +1 to make them look better
        colors = []
        for i in range(0, ncolors):
            colors.append((800 + int(100.0 / ncolors) * i) + 1)
    if colors:
        for color in colors:
            yield color
    else:
        raise ValueError("Unknonw palette")


class Figure(object):
    def __init__(self):
        # User settable parameters:
        self.title = ''
        self.xtitle = ''
        self.ytitle = ''
        self.plot = self.Plot()
        self.legend = self.Legend()

        # Private:
        self._plottables = []
        self._legend_labels = []
        self.style = Styles.Presentation_full

    class Plot(object):
        logx = False
        logy = False
        palette = 'husl'
        palette_ncolors = 10
        xmin, xmax, ymin, ymax = None, None, None, None

    class Legend(object):
        title = None
        position = 'tl'

    def _create_legend(self):
        nentries = len(self._legend_labels)
        leg = Legend(nentries, leftmargin=0, rightmargin=0, entrysep=0.01,
                     textsize=self.style.legendSize, textfont=43, margin=0.1, )
        if self.legend.title:
            leg.SetHeader(self.legend.title)
        leg.SetBorderSize(0)  # no box
        leg.SetFillStyle(0)   # transparent background of legend TPave(!)
        return leg

    def _theme_plottable(self, obj):
        axes = obj.GetXaxis(), obj.GetYaxis()
        for axis in axes:
            axis.SetLabelSize(self.style.axisLabelSize)
            axis.SetLabelFont(self.style.labelfont)
            axis.SetTitleFont(self.style.titlefont)
            axis.SetTitleSize(self.style.axisTitleSize)
        # yaxis only settings:
        axes[1].SetTitleOffset(self.style.plot_ytitle_offset)
        # apply styles, this might need to get more fine grained
        # markers are avilable in children of TAttMarker
        if isinstance(obj, ROOT.TAttMarker):
            # marker size 1 == 8 px, and never scales with canvas...
            obj.SetMarkerSize(self.style.markerSizepx / 8.0)

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
        c = Canvas(width=self.style.canvasWidth,
                   height=self.style.canvasHeight,
                   size_includes_decorations=True)
        if self.legend.position == 'seperate':
            legend_width = .2
            pad_legend = Pad(1 - legend_width, 0, 1., 1., name="legend")
            pad_legend.SetLeftMargin(0.0)
            pad_legend.SetFillStyle(0)  # make this pad transparent
            pad_legend.Draw()
        else:
            legend_width = 0
        pad_plot = Pad(0., 0., 1 - legend_width, 1., name="plot", )
        pad_plot.SetMargin(*self.style.plot_margins)
        pad_plot.Draw()
        pad_plot.cd()

        xmin, xmax, ymin, ymax = get_limits(self._plottables, logx=self.plot.logx, logy=self.plot.logy)
        # overwrite these ranges if defaults are given
        if self.plot.xmin is not None:
            xmin = self.plot.xmin
        if self.plot.xmax is not None:
            xmax = self.plot.xmax
        if self.plot.ymax is not None:
            ymax = self.plot.ymax
        if self.plot.ymin is not None:
            ymin = self.plot.ymin

        colors = get_color_generator(self.plot.palette, self.plot.palette_ncolors)

        is_first = True
        for i, obj in enumerate(self._plottables):
            obj.markerstyle = 'circle'
            try:
                color = next(colors)
            except StopIteration:
                log.warning("Ran out of colors; defaulting to black")
                color = 1
            obj.color = color

            xaxis = obj.GetXaxis()
            yaxis = obj.GetYaxis()

            obj.SetMinimum(ymin)
            obj.SetMaximum(ymax)
            yaxis.SetLimits(ymin, ymax)  # for unbinned data
            yaxis.SetRangeUser(ymin, ymax)

            xaxis.SetLimits(xmin, xmax)

            xaxis.SetTitle(self.xtitle)
            yaxis.SetTitle(self.ytitle)

            xtick_length = xaxis.GetTickLength()
            ytick_length = yaxis.GetTickLength()

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
        pad_plot.SetLogx(self.plot.logx)
        pad_plot.SetLogy(self.plot.logy)

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
                leg.y2 = 1 - pad_plot.GetTopMargin() - ytick_length
                leg.y1 = leg.y2 - leg_hight
            elif self.legend.position.startswith('b'):
                leg_hight = leg.y2 - leg.y1
                leg.y1 = pad_plot.GetBottomMargin() + ytick_length
                leg.y2 = leg.y1 + leg_hight
            # horizontal:
            if self.legend.position[1:].startswith('l'):
                leg_width = 0.3
                leg.x1 = pad_plot.GetLeftMargin() + xtick_length
                leg.x2 = leg.x1 + leg_width
            elif self.legend.position[1:].startswith('r'):
                leg_width = 0.3
                leg.x2 = 1 - pad_plot.GetRightMargin() - xtick_length
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
        success = c.Write()
        if success == 0:
            raise ValueError("Could not write to file!")
        return f

    def save_to_file(self, path, name):
        """
        Save the current figure to the given root file under the given path
        Parameters
        ----------
        name : string
            Name of the file including its extension
        path : string
            Path excluding the file name, relative files are interpreted relative to the working dir
        """
        # check if the name has the right extension
        if len(name.split('.')) != 2:
            raise ValueError("Filename must be given with extension")
        if name.split('.')[1] != 'pdf':
            raise NotImplementedError("Only PDF export is implemented at the moment")
        # strip of tailing / if any
        # this is not compatible with windows, I guess!
        if path.endswith('/'):
            path = path[:-1]
        c = self.draw_to_canvas()
        try:
            os.makedirs(path)
        except OSError:
            pass

        # The order of the following is important! First, set paper size, then draw the canvas and then create the pdf
        # Doin pdf.Range(10, 10) is not sufficient. it just does random shit
        # Be careful to reset the global gStyle when we are finished. Yeah! Globals!
        # Ok, Root does not like that either...
        # paper_width, paper_height = ROOT.Double(), ROOT.Double()
        # ROOT.gStyle.GetPaperSize(paper_width, paper_height)
        ROOT.gStyle.SetPaperSize(self.style.canvasWidth / self.style.pt_per_cm,
                                 self.style.canvasHeight / self.style.pt_per_cm,)
        c = self.draw_to_canvas()
        pdf = ROOT.TPDF("{}/{}".format(path, name))
        c.Draw()
        pdf.Close()

        # reset the page size
        # ROOT.gStyle.SetPaperSize(paper_width, paper_height)

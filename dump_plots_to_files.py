"""
Dump all plots in a root file to individual (pdf) files, while preserving the folder structure
"""
import sys
import os
import shutil

from wand.image import Image

from ROOT import TDirectory, TCanvas
from rootpy import asrootpy
from rootpy.io import root_open


def objects(in_folder):
    """
    Parameters
    ----------
    in_folder : TDirectory
        File or directory of which the objects should be returned
    Yields
    -------
    Tuple :
        (Canvas, root_file_path)
    """
    folder = asrootpy(in_folder)
    for obj in folder.objects():
        if isinstance(obj, TDirectory):
            # recursivly yield objects in sub-folders
            for obj in objects(obj):
                yield obj
        if isinstance(obj, TCanvas):
            yield (obj, folder.GetPath().split(':')[1])

if __name__ == '__main__':
    f_name = sys.argv[1]
    base_path = './' + f_name.split('.')[0] + '/'

    shutil.rmtree(base_path, ignore_errors=True)
    in_folder = root_open(f_name, 'read').MultEstimators.results_post
    for c, rfile_path in objects(in_folder):
        path = base_path + rfile_path
        path = path.replace('//', '/')
        try:
            os.makedirs(path)
        except OSError:
            pass
        filename = path + '/' + c.GetName()
        c.SaveAs(filename + ".pdf")
        with Image(filename=filename + ".pdf") as img:
            img.save(filename=filename + ".png")

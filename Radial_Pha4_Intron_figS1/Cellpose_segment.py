# ─── SCRIPT PARAMETERS ──────────────────────────────────────────────────────────

#@ File(label="Path to your image", description="Path where your image is") input_file_path
#@ File(label="Path to cellpose environment", style="directory", description="Path to conda env") cellpose_env
#@ Double(label="Object diameter", value=1) chr_diameter_calibrated
#@ String(choices={"flow_threshold","cellproba_threshold"}, style="listBox", label="Value to try") variable_to_change
#@ Double(label="Starting value", value=0.0) start_value
#@ Double(label="Stopping value", value=0.0) stop_value
#@ Double(label="Step", value=0.5) step_value
#@ CommandService command

# ─── REQUIREMENTS ───────────────────────────────────────────────────────────────

# List of update sites needed for the code
# * BIOP

# ─── Imports ──────────────────────────────────────────────────────────────────

import os
import sys
from ij import IJ, Prefs, ImagePlus
from ij.gui import WaitForUserDialog
from ij.plugin import Duplicator, ZProjector, ImageCalculator

# Bioformats imports
from loci.plugins import BF, LociExporter
from loci.plugins.in import ImporterOptions
from loci.plugins.out import Exporter

from itertools import count, takewhile

from ch.epfl.biop.wrappers.cellpose.ij2commands import Cellpose_SegmentImgPlusAdvanced, CellposePrefsSet


# ─── Functions ────────────────────────────────────────────────────────────────

def BFImport(indivFile):
    """Import using BioFormats

    Parameters
    ----------
    indivFile : str
        Path of the file to open

    Returns
    -------
    imps : ImagePlus
        Image opened via BF
    """
    options = ImporterOptions()
    options.setId(str(indivFile))
    options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE)
    return BF.openImagePlus(options)

def frange(start, stop, step):
    return takewhile(lambda x: x< stop, count(start, step))

def BFExport(imp, savepath):
    """Export using BioFormats

    Parameters
    ----------
    imp : ImagePlus
        ImagePlus of the file to save
    savepath : str
        Path where to save the image

    """

    # print('Savepath: ', savepath)
    paramstring = "outfile=[" + savepath + "] windowless=true compression=Uncompressed saveROI=true"
    plugin     = LociExporter()
    plugin.arg = paramstring

    exporter   = Exporter(plugin, imp)
    exporter.run()

def drange(start, end, increment, round_decimal_places=None):
    result = []
    if start == end:
        return [start]
    if start < end:
        # Counting up, e.g. 0 to 0.4 in 0.1 increments.
        if increment < 0:
            raise Exception("Error: When counting up, increment must be positive.")
        while start <= end:
            result.append(start)
            start += increment
            if round_decimal_places is not None:
                start = round(start, round_decimal_places)
    else:
        # Counting down, e.g. 0 to -0.4 in -0.1 increments.
        if increment > 0:
            raise Exception("Error: When counting down, increment must be negative.")
        while start >= end:
            result.append(start)
            start += increment
            if round_decimal_places is not None:
                start = round(start, round_decimal_places)
    return result

# ─── Main Code ────────────────────────────────────────────────────────

IJ.log("\\Clear")
IJ.log("STARTING")

command.run(CellposePrefsSet, False,
            "envType", "conda",
            "cellposeEnvDirectory", cellpose_env,
            "version", "2.0")

# Retrieve list of files
IJ.log("Opening the file")
input_file_path   = str(input_file_path)
imp = BFImport(input_file_path)[0]

cal = imp.getCalibration()

folder   = os.path.dirname(input_file_path)
basename = os.path.basename(imp.getTitle())
basename = os.path.splitext(basename)[0]

pxl_diameter = int(round(chr_diameter_calibrated / cal.pixelWidth))

for i in drange(start_value, stop_value, step_value, round_decimal_places=4):
    IJ.log("Now running cellpose with " + variable_to_change + " set to " + str(i))
    res = command.run(Cellpose_SegmentImgPlusAdvanced, False,
                        "imp", imp, "diameter", pxl_diameter,
                        "model", "cyto2", "nuclei_channel", 0,
                        "cyto_channel", 1,
                        variable_to_change, i,
                        "dimensionMode", "3D").get()

    imp_label_cells = res.getOutput("cellpose_imp")
    out_imp = os.path.join(folder, basename + "_" + variable_to_change + "_" + str(i) + ".tif")
    BFExport(imp_label_cells, out_imp)
    imp_label_cells.close()

imp.close()
IJ.log("Script finished")

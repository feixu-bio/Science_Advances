#@ File (label="Select the directory containing the images", style="directory") srcFolder
#@ String (label="Enter the image file extension (e.g. '.vsi' or '.tif')", value=".vsi") imageExt
#@ Integer (label="Select the channel to process") srcChannel
#@ Double (label="Diameter of the object to detect", min=0, value=1.000) blobDiameter
#@ Double (label="Quality threshold the object to detect", value=20000) qualityThresh
#@ String (label="Enter the spot file suffix (e.g. '_nucleus' or '_spot')", value="_nucleus") objectBioType
#@ Boolean (label="Enable batch process (no visual control)") batchMode

import sys
import os
 
from ij import IJ
from ij import WindowManager
from datetime import datetime as dt

# Bio-formats imports
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
from loci.formats import ImageReader, MetadataTools, Memoizer

from java.io import File 

from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import Settings
from fiji.plugin.trackmate import TrackMate
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate import Logger
from fiji.plugin.trackmate.io import TmXmlWriter
from fiji.plugin.trackmate.util import LogRecorder
from fiji.plugin.trackmate.detection import LogDetectorFactory
from fiji.plugin.trackmate.tracking import LAPUtils
from fiji.plugin.trackmate.tracking.sparselap import SparseLAPTrackerFactory
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettingsIO
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettings
from fiji.plugin.trackmate.visualization.table import AllSpotsTableView 
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer as HyperStackDisplayer
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter


############################### Local functions ###############################
def initializeFiji4TM():
	# We have to do the following to avoid errors with UTF8 chars generated in 
	# TrackMate that will mess with our Fiji Jython.
	reload(sys)
	sys.setdefaultencoding('utf-8')
	IJ.run("Close All", "")	

def findAllImages( path, imageExt ):
	"""
	Find all relevant images in the provided folder
	"""
	# We collect all files to open in a list.
	path_to_images = []
	# Replacing some abbreviations (e.g. $HOME on Linux).
	path = os.path.expanduser(path)
	path = os.path.expandvars(path)
	
	# Assumes a non recursive search
	for file_name in os.listdir(path):
	    full_path = os.path.join(path, file_name)
	    if os.path.isfile(full_path):
	        if file_name.endswith(imageExt):
				path_to_images.append(full_path)
				

	return path_to_images

"""
def loadImage( path ):
	# Get currently selected image
	# imp = WindowManager.getCurrentImage()
	imp = IJ.openImage('https://fiji.sc/samples/FakeTracks.tif')
	imp.show()
	print(type(imp))
	return imp
"""

def open_single_series_ch_with_BF(path_to_file, series_number, channel_number):
    options = ImporterOptions()
    options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE)
    options.setSeriesOn(series_number-1, True) # python starts at 0
    options.setSpecifyRanges(True)
    options.setCBegin(series_number-1, channel_number-1) # python starts at 0
    options.setCEnd(series_number-1, channel_number-1)
    options.setCStep(series_number-1, 1)
    options.setId(path_to_file)
    imps = BF.openImagePlus(options) # is an array of imp with one entry
    return imps[0]

def createModel():
	# Some of the parameters we configure below need to have
	# a reference to the model at creation. So we create an
	# empty model now.
	 
	model = Model()
	 
	# Send all messages to ImageJ log window.
	model.setLogger(Logger.IJ_LOGGER)
	return model

def prepSettings( imp, blobDiameter, qualityThresh ):
	settings = Settings( imp )
 
	# Configure detector - We use the Strings for the keys
	settings.detectorFactory = LogDetectorFactory()
	settings.detectorSettings = {
	    'DO_SUBPIXEL_LOCALIZATION' : True,
	    'RADIUS' : blobDiameter/2,
	    'TARGET_CHANNEL' : 1,
	    'THRESHOLD' : qualityThresh,
	    'DO_MEDIAN_FILTERING' : True,
	}  
	 
	## Configure spot filters - Classical filter on quality
	#filter1 = FeatureFilter('QUALITY', qualityThresh, True)
	#settings.addSpotFilter(filter1)
	 
	# Configure tracker => Configure the tracker just enough that it doesn't crash
	settings.trackerFactory = SparseLAPTrackerFactory()
	settings.trackerSettings = LAPUtils.getDefaultLAPSettingsMap() # almost good enough
	settings.trackerSettings['ALLOW_TRACK_SPLITTING'] = True
	settings.trackerSettings['ALLOW_TRACK_MERGING'] = True

	# Add the analyzers for some spot features.
	# Here we decide brutally to add all of them.
	settings.addAllAnalyzers()
	 
	"""
	# Add ALL the feature analyzers known to TrackMate. They will 
	# yield numerical features for the results, such as speed, mean intensity etc.
	settings.addAllAnalyzers()
	 
	# Configure track filters - We want to get rid of the two immobile spots at
	# the bottom right of the image. Track displacement must be above 10 pixels.
	 
	filter2 = FeatureFilter('TRACK_DISPLACEMENT', 10, True)
	settings.addTrackFilter(filter2)
	"""
	
	return settings

def sanityCheck(trackmate):

	ok = trackmate.checkInput()
	if not ok:
	    sys.exit(str(trackmate.getErrorMessage()))
	ok = trackmate.process()
	#trackmate.computeSpotFeatures( True )
	if not ok:
		print(str(trackmate.getErrorMessage()))
	    # sys.exit(str(trackmate.getErrorMessage()))
	return ok
 
def displayResults( model, imp ):
	# A selection.
	selectionModel = SelectionModel( model )
	 
	# Read the default display settings.
	ds = DisplaySettingsIO.readUserDefault()
	 
	displayer =  HyperStackDisplayer( model, selectionModel, imp, ds )
	displayer.render()
	displayer.refresh()
	 
	# Echo results with the logger we set at start:
	model.getLogger().log( str( model ) )
	return

def save2csv( model, sm, ds, savePath ):
	allSTV = AllSpotsTableView(model, sm, ds)
	allSTV.exportToCsv(savePath)
	print( "Results saved to: " + savePath )
	return

def save2xml( settings, logger, trackmate, savePath ):	
	writer = TmXmlWriter( savePath, logger )
	writer.appendLog( logger.toString() )
	writer.appendModel( trackmate.getModel() )
	writer.appendSettings( trackmate.getSettings() )
	writer.writeToFile()
	print( "Results saved to: " + savePath.toString() )
	return
	
############################### Main ###############################

# System initialization
initializeFiji4TM()

# Append all image paths to process
path_to_images = findAllImages( srcFolder.getAbsolutePath(), imageExt )

# Process each image independantly
for path in path_to_images:

	print( "Processing image: " + File(path).getName() )
	
	# Logger -> content will be saved in the XML file.
	logger = LogRecorder( Logger.VOID_LOGGER )
	logger.log( 'TrackMate analysis script\n' )
	dt_string = dt.now().strftime("%d/%m/%Y %H:%M:%S")
	logger.log( dt_string + '\n\n' )
	
	# Load image 
	#imp = loadImage(path)
	imp = open_single_series_ch_with_BF(path, 1, srcChannel)
	
	#----------------------------
	# Create the empty model object now
	#----------------------------
	model = createModel()
	
	#------------------------
	# Prepare settings object
	#------------------------
	settings = prepSettings( imp, blobDiameter, qualityThresh )
	
	#-------------------
	# Instantiate plugin
	#------------------- 
	trackmate = TrackMate( model, settings )
	# get the spot statistics

	 
	#--------
	# Process
	#--------
	ok = sanityCheck( trackmate ) # Make sure input and process are ok
	if not ok:
		print("Skipping file!")
		continue

	#----------------
	# Display results
	#----------------
	if batchMode==False:
		displayResults( model, imp )
			
	# Save output as csv
	savePathCsv = (objectBioType+'.csv').join(path.rsplit(imageExt, 1))
	save2csv( trackmate.getModel(), SelectionModel(model), DisplaySettings(), savePathCsv )

	# Save output as XML
	savePathXml = File((objectBioType+'.xml').join(path.rsplit(imageExt, 1)))
	save2xml(settings, logger, trackmate, savePathXml )	
	
	if batchMode==True:
		IJ.run("Close All", "")
		
	IJ.run("Collect Garbage", "")
			
print("Done")
IJ.log("Done")

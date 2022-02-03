
#--------------------------------------------------------------

# Global timestep output options
timeStepToStartOutputAt=0
forceOutputAtFirstCall=False

# Global screenshot output options
imageFileNamePadding=0
rescale_lookuptable=False

# Whether or not to request specific arrays from the adaptor.
requestSpecificArrays=False

# a root directory under which all Catalyst output goes
rootDirectory='data'

# makes a cinema D index table
make_cinema_table=False

#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# paraview version 5.8.1
#--------------------------------------------------------------

from paraview.simple import *
from paraview import coprocessing

# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 5.8.1

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      # trace generated using paraview version 5.8.1
      #
      # To ensure correct image size when batch processing, please search 
      # for and uncomment the line `# renderView*.ViewSize = [*,*]`

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # Create a new 'Render View'
      renderView1 = CreateView('RenderView')
      renderView1.ViewSize = [1144, 793]
      renderView1.InteractionMode = '2D'
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.CenterOfRotation = [9.99999, 5.0000055, 0.0]
      renderView1.StereoType = 'Crystal Eyes'
      renderView1.CameraPosition = [9.99999, 5.0000055, 10000.0]
      renderView1.CameraFocalPoint = [9.99999, 5.0000055, 0.0]
      renderView1.CameraFocalDisk = 1.0
      renderView1.CameraParallelScale = 9.239944961079138

      # init the 'GridAxes3DActor' selected for 'AxesGrid'
      renderView1.AxesGrid.Visibility = 1

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='B_%t.png', freq=100, fittoscreen=0, magnification=1, width=1144, height=793, cinema={}, compression=5)
      renderView1.ViewTime = datadescription.GetTime()

      SetActiveView(None)

      # ----------------------------------------------------------------
      # setup view layouts
      # ----------------------------------------------------------------

      # create new layout object 'Layout #1'
      layout1 = CreateLayout(name='Layout #1')
      layout1.AssignView(0, renderView1)

      # ----------------------------------------------------------------
      # restore active view
      SetActiveView(renderView1)
      # ----------------------------------------------------------------

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'Legacy VTK Reader'
      # create a producer from a simulation input
      particles = coprocessor.CreateProducer(datadescription, 'particles')

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from particles
      particlesDisplay = Show(particles, renderView1, 'UniformGridRepresentation')

      # get color transfer function/color map for 'B'
      bLUT = GetColorTransferFunction('B')
      bLUT.EnableOpacityMapping = 1
      bLUT.RGBPoints = [0.019445301344837917, 0.278431372549, 0.278431372549, 0.858823529412, 0.020610226289046678, 0.0, 0.0, 0.360784313725, 0.021767004904974258, 0.0, 1.0, 1.0, 0.0229400761774642, 0.0, 0.501960784314, 0.0, 0.024096854793391777, 1.0, 1.0, 0.0, 0.025261779737600538, 1.0, 0.380392156863, 0.0, 0.0264267046818093, 0.419607843137, 0.0, 0.0, 0.02759162962601806, 0.878431372549, 0.301960784314, 0.301960784314]
      bLUT.ColorSpace = 'RGB'
      bLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 'B'
      bPWF = GetOpacityTransferFunction('B')
      bPWF.Points = [0.019445301344837917, 1.0, 0.5, 0.0, 0.02759162962601806, 0.0, 0.5, 0.0]
      bPWF.ScalarRangeInitialized = 1

      # trace defaults for the display properties.
      particlesDisplay.Representation = 'Surface'
      particlesDisplay.ColorArrayName = ['POINTS', 'B']
      particlesDisplay.LookupTable = bLUT
      particlesDisplay.OSPRayScaleArray = 'B'
      particlesDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
      particlesDisplay.SelectOrientationVectors = 'B'
      particlesDisplay.ScaleFactor = 1.9999980000000002
      particlesDisplay.SelectScaleArray = 'B'
      particlesDisplay.GlyphType = 'Arrow'
      particlesDisplay.GlyphTableIndexArray = 'B'
      particlesDisplay.GaussianRadius = 0.0999999
      particlesDisplay.SetScaleArray = ['POINTS', 'B']
      particlesDisplay.ScaleTransferFunction = 'PiecewiseFunction'
      particlesDisplay.OpacityArray = ['POINTS', 'B']
      particlesDisplay.OpacityTransferFunction = 'PiecewiseFunction'
      particlesDisplay.DataAxesGrid = 'GridAxesRepresentation'
      particlesDisplay.PolarAxes = 'PolarAxesRepresentation'
      particlesDisplay.ScalarOpacityUnitDistance = 1.7888915082979384
      particlesDisplay.ScalarOpacityFunction = bPWF
      particlesDisplay.SliceFunction = 'Plane'

      # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
      particlesDisplay.ScaleTransferFunction.Points = [-0.01951385848224163, 0.0, 0.5, 0.0, 0.019511915743350983, 1.0, 0.5, 0.0]

      # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
      particlesDisplay.OpacityTransferFunction.Points = [-0.01951385848224163, 0.0, 0.5, 0.0, 0.019511915743350983, 1.0, 0.5, 0.0]

      # init the 'Plane' selected for 'SliceFunction'
      particlesDisplay.SliceFunction.Origin = [9.99999, 5.0000055, 0.0]

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for bLUT in view renderView1
      bLUTColorBar = GetScalarBar(bLUT, renderView1)
      bLUTColorBar.Title = 'B'
      bLUTColorBar.ComponentTitle = 'Magnitude'

      # set color bar visibility
      bLUTColorBar.Visibility = 1

      # show color legend
      particlesDisplay.SetScalarBarVisibility(renderView1, True)

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(particles)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'particles': [1]}
  coprocessor.SetUpdateFrequencies(freqs)
  if requestSpecificArrays:
    arrays = [['B', 0]]
    coprocessor.SetRequestedArrays('particles', arrays)
  coprocessor.SetInitialOutputOptions(timeStepToStartOutputAt,forceOutputAtFirstCall)

  if rootDirectory:
      coprocessor.SetRootDirectory(rootDirectory)

  if make_cinema_table:
      coprocessor.EnableCinemaDTable()

  return coprocessor


#--------------------------------------------------------------
# Global variable that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
coprocessor = CreateCoProcessor()

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView and the update frequency
coprocessor.EnableLiveVisualization(True, 1)

# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    coprocessor.LoadRequestedData(datadescription)

# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global coprocessor

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=rescale_lookuptable,
        image_quality=0, padding_amount=imageFileNamePadding)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)

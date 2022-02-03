
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
make_cinema_table=True

#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# paraview version 5.8.0
#--------------------------------------------------------------

from paraview.simple import *
from paraview import coprocessing

# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 5.8.0

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      # trace generated using paraview version 5.8.0
      #
      # To ensure correct image size when batch processing, please search 
      # for and uncomment the line `# renderView*.ViewSize = [*,*]`

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # Create a new 'Render View'
      renderView1 = CreateView('RenderView')
      renderView1.ViewSize = [2282, 1084]
      renderView1.InteractionMode = '2D'
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.CenterOfRotation = [9.99999, 5.0000055, 0.0]
      renderView1.StereoType = 'Crystal Eyes'
      renderView1.CameraPosition = [9.99999, 5.0000055, 10000.0]
      renderView1.CameraFocalPoint = [9.99999, 5.0000055, 0.0]
      renderView1.CameraFocalDisk = 1.0
      renderView1.CameraParallelScale = 6.311006735249734

      # init the 'GridAxes3DActor' selected for 'AxesGrid'
      renderView1.AxesGrid.Visibility = 1

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='RhoE0_%t.png', freq=100, fittoscreen=0, magnification=1, width=2282, height=1084, cinema={}, compression=5)
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

      # get color transfer function/color map for 'rhoe'
      rhoeLUT = GetColorTransferFunction('rhons0')
      rhoeLUT.EnableOpacityMapping = 1
      rhoeLUT.RGBPoints = [-1.0138441324234009, 0.278431372549, 0.278431372549, 0.858823529412, -0.8688644235045055, 0.0, 0.0, 0.360784313725, -0.7248985587039243, 0.0, 1.0, 1.0, -0.5789050056667151, 0.0, 0.501960784314, 0.0, -0.4349391408661337, 1.0, 1.0, 0.0, -0.2899594319472384, 1.0, 0.380392156863, 0.0, -0.14497972302834317, 0.419607843137, 0.0, 0.0, -1.4109447832311162e-08, 0.878431372549, 0.301960784314, 0.301960784314]
      rhoeLUT.ColorSpace = 'RGB'
      rhoeLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 'rhoe'
      rhoePWF = GetOpacityTransferFunction('rhons0')
      rhoePWF.Points = [-1.0138441324234009, 1.0, 0.5, 0.0, -1.4109447832311162e-08, 0.3482142984867096, 0.5, 0.0]
      rhoePWF.ScalarRangeInitialized = 1

      # trace defaults for the display properties.
      particlesDisplay.Representation = 'Surface'
      particlesDisplay.ColorArrayName = ['POINTS', 'rhons0']
      particlesDisplay.LookupTable = rhoeLUT
      particlesDisplay.OSPRayScaleArray = 'rhons0'
      particlesDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
      particlesDisplay.SelectOrientationVectors = 'None'
      particlesDisplay.ScaleFactor = 1.9999980000000002
      particlesDisplay.SelectScaleArray = 'rhons0'
      particlesDisplay.GlyphType = 'Arrow'
      particlesDisplay.GlyphTableIndexArray = 'rhons0'
      particlesDisplay.GaussianRadius = 0.0999999
      particlesDisplay.SetScaleArray = ['POINTS', 'rhons0']
      particlesDisplay.ScaleTransferFunction = 'PiecewiseFunction'
      particlesDisplay.OpacityArray = ['POINTS', 'rhons0']
      particlesDisplay.OpacityTransferFunction = 'PiecewiseFunction'
      particlesDisplay.DataAxesGrid = 'GridAxesRepresentation'
      particlesDisplay.PolarAxes = 'PolarAxesRepresentation'
      particlesDisplay.ScalarOpacityUnitDistance = 1.7888915082979384
      particlesDisplay.ScalarOpacityFunction = rhoePWF
      particlesDisplay.IsosurfaceValues = [-0.5069220732664244]
      particlesDisplay.SliceFunction = 'Plane'

      # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
      particlesDisplay.ScaleTransferFunction.Points = [-1.0138441324234009, 0.0, 0.5, 0.0, -1.4109447832311162e-08, 1.0, 0.5, 0.0]

      # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
      particlesDisplay.OpacityTransferFunction.Points = [-1.0138441324234009, 0.0, 0.5, 0.0, -1.4109447832311162e-08, 1.0, 0.5, 0.0]

      # init the 'Plane' selected for 'SliceFunction'
      particlesDisplay.SliceFunction.Origin = [9.99999, 5.0000055, 0.0]

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for rhoeLUT in view renderView1
      rhoeLUTColorBar = GetScalarBar(rhoeLUT, renderView1)
      rhoeLUTColorBar.Title = 'rhoe'
      rhoeLUTColorBar.ComponentTitle = ''

      # set color bar visibility
      rhoeLUTColorBar.Visibility = 1

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
    arrays = [['rhons0', 0]]
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

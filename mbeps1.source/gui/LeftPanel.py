import wx
import wx.stc as stc
import numpy as NP
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from subprocess import Popen, PIPE
import Image
import os

from Events import *

#class for the frame handling 
class RecordPanel(wx.Frame):

	currentlyWritingFiles = [] #Class variable, containing names of files currently recording

	def __init__(self,parent):
		wx.Frame.__init__(self, parent, -1, 'Movie Options',style=wx.FRAME_FLOAT_ON_PARENT|wx.DEFAULT_FRAME_STYLE)
		self.Bind(wx.EVT_CLOSE, self.OnClose)
		self.parent = parent
		self.SetupControls()
		self.Show()

	def SetupControls(self):
		self.sizer = wx.BoxSizer(orient=wx.VERTICAL)
		hs1 = wx.BoxSizer(orient=wx.HORIZONTAL)

		rec, text, overwrite = self.parent.getRecordStatus()

		self.cbox2 = wx.CheckBox(self, -1, "Overwrite Existing Movie File")
		self.cbox2.SetValue(overwrite)
		self.sizer.Add(item=self.cbox2)

		lbl1 = wx.StaticText(self, label="Movie File Name:")
		self.saveName = wx.TextCtrl(self, size=(140, -1),value=text)
		hs1.Add(item=lbl1)
		hs1.Add(item=self.saveName)

		self.cbox1 = wx.CheckBox(self, -1, "Check to begin recording")
		self.cbox1.SetValue(rec)
		self.sizer.Add(item=self.cbox1)
		self.sizer.Add(item=hs1)

		self.SetSizerAndFit(self.sizer)

	def OnClose(self, event):
		svname = self.saveName.GetValue()
		recval = self.cbox1.GetValue()
		overwrite = self.cbox2.GetValue()
		"""if not self.addActiveWritingFile(svname) and recval:
			wx.MessageBox('Another video is currently recording to file '+svname + ".  Not recording.", 'Cannot Record', wx.OK | wx.ICON_INFORMATION)
		elif recval and self.fileExists(svname) and not overwrite:
			True"""
		if recval and self.fileExists(svname) and not overwrite: #check if we can overwrite
			wx.MessageBox('File '+svname + " already exists.  Enable overwrite.", 'Cannot Record', wx.OK | wx.ICON_INFORMATION)
		else:  #No overwriting flags
			if not self.addActiveWritingFile(svname) and recval: #currently writing to stream
				wx.MessageBox('Another video is currently recording to file '+svname + ".  Not recording.", 'Cannot Record', wx.OK | wx.ICON_INFORMATION)
			else:
				self.parent.setRecordStatus(recval , svname, overwrite )
		if not recval:  #Stopping recording
			try:
				self.currentlyWritingFiles.remove(svname)
			except ValueError:
				True
		self.Hide()
		self.Destroy()

	def addActiveWritingFile(self, fname):
		if fname in self.currentlyWritingFiles:
			return False
		else:
			self.currentlyWritingFiles.append(fname)
			return True

	def fileExists(self, fname):
		return os.path.isfile(fname)


class LeftPanel(wx.Panel):
	def __init__(self, parent):
		self.mainframe = parent
		wx.Panel.__init__(self, parent, -1, wx.DefaultPosition, wx.DefaultSize)
		vsizer1 = wx.BoxSizer(orient=wx.VERTICAL)
		self.createGraph()
		vsizer1.Add(item=self.mycanvas, proportion=1, flag=wx.EXPAND | wx.ALL, border=0)
		vsizer1.Add(item=self.toolbar, flag=wx.EXPAND | wx.ALL, border=0)
		self.SetSizerAndFit(vsizer1)
		self.slopeStack = []
		self.measuring = False
		self.currentEvent = None
		self.arbGraphParameters = {"axesType":"Linear-Linear"}
		self.newCP = None
		#movie recording defaults
		self.recordVideo = False
		self.movieFileName = "Movie.avi"
		self.overwrite = False
		self.dpi = 100 #movie dots per inch
		self.moviePipe = None
		self.newframe = parent
		self.persistentVars = dict()

		EVT_CLOSEOP(self, self.OnCloseCP)
		EVT_REFRESHGRAPH(self, self.OnRefreshGraph)

	def setRecordStatus(self, record, fname, overwrite):
		self.movieFileName = fname
		self.overwrite = overwrite
		if (not self.recordVideo) and record:  #If record was off and changed to on
			#self.movieWriter.setup(self.figure, self.movieFileName, self.dpi)
			fps = 15
			#'mjpeg'
			self.moviePipe = Popen(['ffmpeg','-f', 'image2pipe','-framerate','24', '-vcodec', 'png', '-i','-', '-c:v', 'libx264', '-r','24','-y',self.movieFileName], stdin=PIPE)
			self.recordingSize = self.mycanvas.GetSize()
		if self.recordVideo and not record:  #Turned off record
			#self.movieWriter.finish()
			self.moviePipe.stdin.flush()
			self.moviePipe.stdin.close()
			self.moviePipe.wait()

		self.recordVideo = record
		if record:
			self.recButton.SetBitmapLabel(self.RecOnBmp)
		else:
			self.recButton.SetBitmapLabel(self.RecOffBmp)

	def getRecordStatus(self):
		return self.recordVideo, self.movieFileName, self.overwrite

	def resetGraph(self):
		self.figure.delaxes(self.axes)
		self.figure.clf()
		self.axes = self.figure.add_subplot(111)
		

	def createGraph(self):
		self.figure = matplotlib.figure.Figure()
		self.axes = self.figure.add_subplot(111)
		t = NP.arange(0.0,10,1.0)
		s = [0,1,0,1,0,2,1,2,1,0]
		self.y_max = 10
		self.mycanvas = FigureCanvas(self,-1,self.figure)
		self.mycanvas.SetSize((100,100))
		#self.axes.plot(t,s)
		#self.axes.imshow(mpimg.imread('gui/default.png'), interpolation='nearest', aspect='auto')
		self.mycanvas.mpl_connect('button_press_event',self.onclick)
		self.toolbar = wx.BoxSizer(orient=wx.HORIZONTAL)
		self.navMenu = NavigationToolbar2Wx(self.mycanvas)
		self.toolbar.Add(item=self.navMenu )
		self.slopeButton = wx.ToggleButton(self,-1,"Measure Slope")  #create slope measurement button
		self.toolbar.Add(item=self.slopeButton )

		self.graphOptionsButton = wx.Button(self,-1,"Graph Options")
		self.toolbar.Add(item=self.graphOptionsButton)

		self.RecOnBmp = wx.Bitmap("./gui/rec.png", wx.BITMAP_TYPE_ANY)
		self.RecOffBmp = wx.Bitmap("./gui/recoff.png", wx.BITMAP_TYPE_ANY)
		self.recButton = wx.BitmapButton(self, id=wx.ID_ANY, bitmap=self.RecOffBmp, size=(self.RecOffBmp.GetWidth()+10, self.RecOnBmp.GetHeight()+10))
		self.toolbar.Add(item=self.recButton)
		self.recButton.Bind(wx.EVT_BUTTON, self.OnRecord)

		self.graphOptionsButton.Bind(wx.EVT_BUTTON, self.OnOptions)
		self.slopeButton.Bind(wx.EVT_TOGGLEBUTTON, self.OnMeasureButton)

		#movie writer stuff
		"""self.FFMpegWriter = manimation.writers['ffmpeg']
		self.metadata = dict(title='Movie Test', artist='Matplotlib', comment='Movie support!')
		self.movieWriter = self.FFMpegWriter(fps=1, metadata=self.metadata)	"""

	def OnResult(self, event):
		"""Show Result status."""
		if event.data is None:
			# Thread aborted (using our convention of None return)
			self.status.SetLabel('Computation aborted')
		else:
			#self.figure = matplotlib.figure.Figure()
			#self.axes = self.figure.add_subplot(111)
			self.currentEvent = event
			self.currentEvent.data._PV = self.persistentVars  #Gives the plot a simple dictionary to save persistent vars
			self.DrawPlot()
			#Save movie frame
			if self.recordVideo:
				print "Writing " + str(self.currentEvent.data.simTime)
				#self.movieWriter.grab_frame()
				self.movieDim = self.mycanvas.get_width_height()
				im = Image.fromstring('RGB', self.movieDim, self.mycanvas.tostring_rgb())
				im = im.resize(self.recordingSize)
				im.save(self.moviePipe.stdin, 'PNG')  #Was JPEG

	def OnChangeGraph(self, event):
		try:
			self.currentEvent.data
			try:
				self.currentEvent._PV
				self.persistentVars = dict()
				self.currentEvent.data._PV = self.persistentVars
			except:
				True
		except:
			True

	def DrawPlot(self):
		try:
			self.currentEvent.data._PV
		except AttributeError:
			self.currentEvent.data._PV = self.persistentVars
		self.axes.cla()
		try:
			self.currentEvent.data.setParams(self.arbGraphParameters)  #pass paramters to plot
		except AttributeError:
			True
		#self.currentEvent.data.setAxesType(self.arbGraphParameters["axesType"])

		#self.figure.delaxes(self.axes)
		#self.figure.clf()
		#self.axes = self.figure.add_subplot(111)

		rax = self.currentEvent.data.drawPlot(self.figure, self.axes)
		if rax != None: #This allows the graph to reconfigure its own axes, instead of allowing the panel to handle it
			self.axes = rax
		self.axes.set_title(self.currentEvent.data.plottype, loc="right")
		
		self.mycanvas.draw()
		self.slopeStack = []
		self.measuring = False

	def PlotLine(self,x1,x2):
		self.axes.plot([x1[0],x2[0]],[x1[1],x2[1]])
		self.mycanvas.draw()

	def OnMeasureButton(self,event):
		if(len(self.slopeStack) == 0 and self.measuring == True):
			self.measuring = False
			self.slopeButton.SetValue(False)
			self.mainframe.status.SetStatusText("I await your command")
			return
		self.DrawPlot()
		self.slopeStack = []
		self.measuring = True
		self.slopeButton.SetValue(True)
		self.mainframe.status.SetStatusText("Click the two points that define your slope")

	def OnRefreshGraph(self,event):
		self.DrawPlot()

	def OnOptions(self, event): #Open options menu
		if self.newCP == None:  #Only allow one options window at a time
			#try: DEBUG:  Comment back in
				self.newCP = self.currentEvent.data.makeControlPanel(self)
			#except AttributeError:
			#	print "No bound control panel"

	def OnCloseCP(self, event):
		self.newCP.Hide()
		self.newCP.Destroy()
		self.newCP = None
		self.DrawPlot()  #refresh changes

	def onclick(self,event):
		if event.inaxes == None or self.measuring == False:
			return
		if len(self.slopeStack) < 2:
			self.slopeStack.append( NP.array([event.xdata,event.ydata]) )
			self.mainframe.status.SetStatusText("First Point selected is "+str((event.xdata,event.ydata) )+ ". Select next point to find slope")
		if len(self.slopeStack) == 2:
			x1 = self.slopeStack[0]
			x2 = self.slopeStack[1]
			mv = x2 - x1
			m = mv[1]/mv[0]
			self.PlotLine(x1,x2)
			self.mainframe.status.SetStatusText("The slope is "+str(m) ) 
			self.slopeStack = []
			self.measuring = False
			self.slopeButton.SetValue(False)

	def OnRecord(self,event):
		ne = RecordPanel(self)

	#Returns true if the graph is frozen (cannot switch to a new graph type)
	#For example, when recording a Velocity Plot, cannot change to a Phase plot
	def isFrozen(self):
		return self.recordVideo

import sys, os
from lib import *
import numpy as NP
import matplotlib.pyplot as plt

import wx
import wx.stc as stc
import threading
import copy
import matplotlib
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from collections import deque

from defaults import *
from NewFrame import *
from RunSim import *
from RightPanel import *
from GraphStack import *
from PipeSim import *
from Signals import *

#TODO:  Use nplots as default window size
#print pyfft1mod.fft1d.pycal()

# Button definitions
ID_START = wx.NewId()
ID_STOP = wx.NewId()		

programDefaults = DefaultLoader("foo.default")

class Hack:
	def __init__(self):
		self.nplot=4
		self.ntv = 1
		self.nts = 1
		self.ntp = 1
		self.ntw = 1

in1 = Hack()


# GUI Frame class that spins off the worker thread
class MainFrame(wx.Frame, Dispatcher, DefaultsCommLink):
	"""Class MainFrame."""
	def __init__(self, parent, id, loader, pipemode=None, que=None, timedir=None, events=None, outq=None):
		"""Create the MainFrame."""
		wx.Frame.__init__(self, parent, id, 'Thread Test')
		self.loader = loader
		self.loader.loadFromFile()
		self.initStack()
		
		self.status = self.CreateStatusBar()	#create status bar
		self.status.SetStatusText("Hola")

		self.rpanel = RightPanel(self,self)	#create interface
		self.sizer = wx.BoxSizer(wx.HORIZONTAL)
		self.sizer.Add(item=self.rpanel, proportion=1, flag=wx.ALL | wx.EXPAND, border=2)
		self.SetSizerAndFit(self.sizer)

		#self.Bind(wx.EVT_CLOSE, self.OnQuit)

		self.pEvents = events  #Communicate with parent process, various GUI events
		self.outq = outq
		self.pipemode = pipemode
		if pipemode is not None:
			self.worker = PipeSimulation(self, pipemode, que, timedir, outq)
		else:
			self.worker = RunSimulation(self)	#Initialize simulation
		EVT_RESULT(self,self.OnResultPre) #Link up events
		EVT_CONTROL(self, self.OnControl)
		EVT_NEWTIME(self, self.OnNewTime) #link up sim time
		EVT_CLEARGRAPHSTACK(self, self.OnClearGraphStack)

		nf = NewFrame(self, self.loader, self)  #create default frame
		self.windowList = [] #list of frames
		self.windowList.append(nf) #add default
		#Now create frame based on in1.nplot (NPLOT in input1 file)
		if in1.nplot == 1:
			nf.OnLayout1(None)
		elif in1.nplot == 2:
			nf.OnLayout2h(None)
		elif in1.nplot == 3:
			nf.OnLayout3(None)
		elif in1.nplot == 4:
			nf.OnLayout4(None)
		elif in1.nplot > 4:
			nf.OnLayout4(None)
			print "There is not currently a layout larger than 4.  Perhaps you should create one"
		#Now figure out which graphs to display based on defaults
		defaultGraphs = []
		if in1.ntv > 0 or pipemode is not None:
			defaultGraphs.append("DRAWVELOCITY")
		if in1.nts > 0 or pipemode is not None:
			defaultGraphs.append("DRAWPHASE")
		if in1.ntp > 0 or pipemode is not None:
			defaultGraphs.append("DRAWPOT")
		if in1.ntw > 0 or pipemode is not None:
			defaultGraphs.append("ENERGY")
		#Now populate the frame with the graphs
		if len(nf.displayAreas) > len(defaultGraphs):
			enumList = defaultGraphs
		else:
			enumList = nf.displayAreas
		for (i,d) in enumerate(enumList):
			print defaultGraphs[i]
			nf.displayAreas[i].setGraphByName(defaultGraphs[i])

		if pipemode is not None:  #Get pipe mode started
			#wx.CallAfter(self.rpanel.OnStartLong, None)  #Press the run continuously button
			#wx.CallAfter(self.Hide) #Hide the control window
			True

	def OnStart(self, event):
		"""Start Computation."""
		# Trigger the worker thread unless it's already busy
		self.worker.run()
		if not self.worker:
			self.status.SetLabel('Starting computation')

	def OnNewTime(self,event):
		self.rpanel.timerText.SetLabel("Time is " + str(event.time))

	def OnReset(self,event):
		self.pEvents.put( ResetSignal() )

	def OnControl(self,event):
		True  #An event to handle GUI control commands

	def OnResultPre(self,event):
		#Find a home for the event
		self.OnResult(event)  #Call the result handler in the Dispatcher mixin

	def OnQuit(self,event):
		if self.pipemode != None:
			self.pipemode.close()
			self.pEvents.close()
			self.que.close()

	def OnClearGraphStack(self,event):
		if hasattr(event,"codename"):
			print "Codename"
		else:
			print "Clearing GS"
			del self.dispatchers[:]  #Delete all objects in dispatchers, defined in GraphStack.py

def restart_program():
    """Restarts the current program.
    Note: this function does not return. Any cleanup action (like
    saving data) must be done before calling this function."""
    python = sys.executable
    os.execl(python, python, * sys.argv)

class MainApp(wx.App):
	"""Class Main App."""
	def __init__(self,arg, pipemode=None,que=None, timedir=None, events=None, outq = None):
		self.pipemode = pipemode
		self.que = que
		self.timeDir = timedir
		self.pEvents = events
		self.outq = outq
		wx.App.__init__(self,arg)

	def OnInit(self):
		"""Init Main App."""
		self.frame = MainFrame(None, -1, programDefaults, pipemode=self.pipemode, que=self.que, timedir=self.timeDir, events=self.pEvents, outq = self.outq)
		self.frame.Show(True)
		self.SetTopWindow(self.frame)
		return True

	def OnExit(self):
		self.pipemode.close()
		self.que.close()
		self.pEvents.close()
		#self.pEvents.put(ExitSignal())
		print "OnExit"

if __name__ == "__main__":
	#main()
	app = MainApp(0)
	app.MainLoop()

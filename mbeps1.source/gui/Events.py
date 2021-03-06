import wx
import wx.stc as stc
import numpy as NP

# Define notification event for thread completion
#Step 1, Define a new event ID
EVT_RESULT_ID = wx.NewId()
EVT_CONTROL_ID = wx.NewId()
EVT_CLOSEOP_ID = wx.NewId()
EVT_REFRESHGRAPH_ID = wx.NewId()
EVT_NEW_TIMESTEP_ID = wx.NewId()
EVT_RUNSTEP_ID = wx.NewId()
EVT_CLEARGRAPHSTACK_ID = wx.NewId()

#Step 2, create a new function to connect the event ID to a callback func
def EVT_REFRESHGRAPH(win, func):
	"""Define Result Event."""
	win.Connect(-1, -1, EVT_REFRESHGRAPH_ID, func)

def EVT_RESULT(win, func):
	"""Define Result Event."""
	win.Connect(-1, -1, EVT_RESULT_ID, func)

def EVT_CONTROL(win, func):
	"""Define Result Event."""
	win.Connect(-1, -1, EVT_CONTROL_ID, func)

def EVT_NEWTIME(win, func):
	"""Define Result Event."""
	win.Connect(-1, -1, EVT_NEW_TIMESTEP_ID, func)

def EVT_CLOSEOP(win, func):
	"""Define Result Event."""
	win.Connect(-1, -1, EVT_CLOSEOP_ID, func)

def EVT_RUNSTEP(win, func):
	win.Connect(-1, -1, EVT_RUNSTEP_ID, func)

def EVT_CLEARGRAPHSTACK(win, func):
	"""Define Result Event."""
	win.Connect(-1, -1, EVT_CLEARGRAPHSTACK_ID, func)

#Define An Event class
class RunStepEvent(wx.PyEvent):
	def __init__(self):
		wx.PyEvent.__init__(self)
		self.SetEventType(EVT_RUNSTEP_ID)

class SimTimeEvent(wx.PyEvent):
	def __init__(self, time):
		wx.PyEvent.__init__(self)
		self.SetEventType(EVT_NEW_TIMESTEP_ID)
		self.time = time

class ClearGraphStackEvent(wx.PyEvent):
	def __init__(self):
		wx.PyEvent.__init__(self)
		self.SetEventType(EVT_CLEARGRAPHSTACK_ID)

class ResultEvent(wx.PyEvent):
	"""Simple event to carry arbitrary result data."""
	def __init__(self, data, time):
		"""Init Result Event."""
		wx.PyEvent.__init__(self)
		self.SetEventType(EVT_RESULT_ID)
		self.data = data
		self.name = data.plottype
		self.data.simTime = time

class ControlEvent(wx.PyEvent):
	"""Simple event to carry arbitrary result data."""
	def __init__(self, data, time):
		"""Init Result Event."""
		wx.PyEvent.__init__(self)
		self.SetEventType(EVT_CONTROL_ID)
		self.data = data
		self.data.simTime = time

class CloseOptionsEvent(wx.PyEvent):
	"""Simple event to carry arbitrary result data."""
	def __init__(self):
		"""Init Result Event."""
		wx.PyEvent.__init__(self)
		self.SetEventType(EVT_CLOSEOP_ID)

class RefreshGraphEvent(wx.PyEvent):
	"""Simple event to carry arbitrary result data."""
	def __init__(self):
		"""Init Result Event."""
		wx.PyEvent.__init__(self)
		self.SetEventType(EVT_REFRESHGRAPH_ID)

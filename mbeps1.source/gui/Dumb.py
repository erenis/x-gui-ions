import wx
import wx.stc as stc

from LeftPanel import *

class Dumb(LeftPanel):
	def __init__(self, parent, dlist):
		self.parent = parent
		self.centralDispatcher = dlist
		#wx.Panel.__init__(self, parent, -1, wx.DefaultPosition, wx.DefaultSize, style=wx.SUNKEN_BORDER)
		LeftPanel.__init__(self,parent)
		self._mouseDownFlag = 0
		self.mycanvas.Bind(wx.EVT_RIGHT_DOWN, self.OnRightDown)
		self.Bind(wx.EVT_RIGHT_DOWN, self.OnRightDown)

	def OnRightDown(self,event):
		menu = wx.Menu()
		if self.isFrozen(): #If we are currently recording video, do not allow chaning of movies
			ti = menu.Append( 1, "Stop recording before changing graph type." )
			ti.Enable(False)
		else:  #Allow graph to change
			for (i,g) in enumerate(self.centralDispatcher):
				ti = menu.Append( i, g.description )
				self.Bind(wx.EVT_MENU,self.PopupHandler, ti)
		self.PopupMenu(menu, (event.GetX(), event.GetY()) )
		menu.Destroy()

	def PopupHandler(self,event):
		self.resetGraph()
		print event.GetId()
		for g in self.centralDispatcher:
			g.RemoveListener(self)
		self.centralDispatcher[event.GetId()].AddListener(self)
		re = self.centralDispatcher[event.GetId()].getRecent()
		if re != None:
			self.currentEvent = re
			self.DrawPlot()

	def setGraphByName(self,name):  #set the graph to display by name
		self.movieFileName = name + ".mp4"
		for g in self.centralDispatcher:
			g.RemoveListener(self)
			try:
				self.OnChangeGraph(None)
			except:
				True
		for g in self.centralDispatcher:
			if g.name == name:
				g.AddListener(self)
				re = g.getRecent()
				if re != None:
					self.currentEvent = re
					self.DrawPlot()
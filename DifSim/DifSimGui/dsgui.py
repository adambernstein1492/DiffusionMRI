#! The imports
#!-------------
#!
#! The MPLFigureEditor is imported from last example.

from __future__ import with_statement

import numpy as np
import math
from subprocess import Popen, PIPE, call, STDOUT
from threading import Thread, Lock
from time import sleep
from scipy import special, swapaxes
import sys
import nipy
import copy
from os import devnull
import os

import wxversion
wxversion.select("2.8")
import wx

#import pygtk
#import gtk

from traits.etsconfig.api import ETSConfig
if len(sys.argv) > 1 and (sys.argv[1] == 'wx' or sys.argv[1] == 'qt4'):
    print sys.argv[1]
    ETSConfig.toolkit = '%s' % (sys.argv[1])
else:
    ETSConfig.toolkit = 'wx'

from pyface.qt import QtGui, QtCore

from traits.api import HasTraits, Float, Int, Property, Str, Instance, Button, on_trait_change, Bool, Enum, Range
from traitsui.api import View, Item, Group, HSplit, Handler, VSplit, UItem, HGroup, VGroup, VGrid, HFlow, RangeEditor
from traitsui.file_dialog import open_file, FileInfo, TextInfo, ImageInfo
from traitsui.menu import NoButtons, Action
from traits.trait_base import not_event


from tvtk.api import tvtk

from mayavi.sources.api import ArraySource
from mayavi.modules.api import Outline, Surface, ContourGridPlane, IsoSurface, Vectors, Glyph
from mayavi.core.ui.api import MayaviScene, MlabSceneModel, SceneEditor
from mayavi.core.ui.engine_view import EngineView
from mayavi.sources.api import ParametricSurface
from mayavi.modules.scalar_cut_plane import ScalarCutPlane
from mayavi.mlab import view, pipeline, quiver3d, mesh, plot3d, triangular_mesh 

import matplotlib

if ETSConfig.toolkit == 'wx':
    matplotlib.use('WXAgg')
elif ETSConfig.toolkit == 'qt4':
    matplotlib.use('Qt4Agg')
    matplotlib.rcParams['backend.qt4']='PySide'

#matplotlib.use('GTKAgg')

from matplotlib.figure import Figure
from matplotlib.ticker import MultipleLocator
from matplotlib.colorbar import Colorbar
from matplotlib.pyplot import subplots
from matplotlib.ticker import MaxNLocator
from matplotlib.artist import setp
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d import axes3d
from matplotlib import gridspec
import matplotlib.cm as cm
#from mpl_figure_editor import MPLFigureEditor

if ETSConfig.toolkit == 'wx':
    from traitsui.wx.basic_editor_factory import BasicEditorFactory 
    from traitsui.wx.editor import Editor
    from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
elif ETSConfig.toolkit == 'qt4':
    from traitsui.qt4.basic_editor_factory import BasicEditorFactory 
    from traitsui.qt4.editor import Editor
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas

#! User interface objects
#!------------------------
#!
#! These objects store information for the program to interact with the
#! user via traitsUI.

class Geometry(HasTraits):
    """Object containing description of geometry"""
    spacing = Float(2., label="spacing", desc="spacing")
    periodic = Int(1, label="periodic", desc="periodic boundary condition (1=on, 0=off)")
    voxel_size = Int(80, label="voxel_size", desc="voxel size in microns, needs to be bigger than largest dim")
    x_dim = Property(depends_on=['spacing','cylinder_diameter'], label="x_dim", desc="x dimension in um") 
    y_dim = Property(depends_on=['spacing','cylinder_diameter'], label="y_dim", desc="y dimension in um") 
    z_dim = Property(depends_on=['spacing','cylinder_diameter'], label="z_dim", desc="z dimension in um") 
    cylinder_array = Int(1, label="cylinder_array", desc="cylinder array")
    cylinder_diameter = Float(7.992, label="cylinder_diameter", desc="cylinder diameter")
    cylinder_radius = Float(4.0, label="cylinder_radius", desc="cylinder radius")
    permeable = Float(1, label="permeable", desc="permeable")
    structure_file = Str("./cylinder-8.mdl", label="structure_file", desc="mcell file describing geometry")
    ncCurrent_strength = Float(1.0, label="ncCurrent_strength", desc="nc current strength")
    ncCurrent_output_file = Str("./ncoutput.txt", label="ncCurrent_output_file", desc="file to store nc output")
    
    def _get_x_dim(self):
        return self.cylinder_diameter + self.spacing
    def _get_y_dim(self):
        return np.sqrt(3)*(self.cylinder_diameter + self.spacing)
    def _get_z_dim(self):
        return 7.

class Simulation(HasTraits):
    """Object that contains the parameters that control the simulation."""
    #difsim_location = Str('/home/ben/src/DifSim/', label="difsim_location", desc="location of DifSim source folder")
    random_seed = Int(100, label="random_seed", desc="random seed")
    nstep = Int(300, label="nstep", desc="number of steps (not used)")
    step_size = Int(25, label="step_size", desc="time step length in microseconds")
    nparts = Int(1, label="nparts", desc="number of diffusing molecules")
    
class Stimulus(HasTraits):
    stimmaker_max_magnitude = Float(1., label="maximum magnitude", desc="maximum magnitude of stimulus")
    stimmaker_phase = Float(0., label="phase", desc="phase of stimulus")
    stimmaker_tottime = Float(1., label = "total stimulus time", desc="total time of stimulus")
    stimmaker_timepoints = Int(200, label = "number of timepoints", desc="total number of time points in stimulus")
    stimmaker_offset = Float(1., label = "offset", desc="offset from the origin")
    stimmaker_frequency = Float(1., label = "frequency", desc="frequency of stimulus")

    stimulus_magnitude = np.array([])
    stimulus_timecourse= np.array([])


    traits_view = View(VGrid(Item('stimmaker_max_magnitude'),
                             Item('stimmaker_offset'),
                             Item('stimmaker_tottime'),
                             Item('stimmaker_frequency'),
                             Item('stimmaker_timepoints'),
                             Item('stimmaker_phase'),))

    def set_timecourse(self):
        self.stimulus_timecourse = np.linspace(0, self.stimmaker_tottime, self.stimmaker_timepoints)

    def sinusoidal_stimulus(self):
        self.set_timecourse()
        self.stimulus_magnitude = self.stimmaker_max_magnitude*np.sin(self.stimulus_timecourse*2*np.pi*self.stimmaker_frequency - np.pi*self.stimmaker_phase) + self.stimmaker_offset

    def sawtooth_stimulus(self):
        self.set_timecourse()
        self.stimulus_magnitude = -2*self.stimmaker_max_magnitude/np.pi*np.arctan(special.cotdg(((180./np.pi)*self.stimulus_timecourse*np.pi*self.stimmaker_frequency - 2*np.pi*self.stimmaker_phase))) + self.stimmaker_offset

    def square_stimulus(self):
        self.set_timecourse()
        self.stimulus_magnitude = 0.5*self.stimmaker_max_magnitude*(np.sign(np.sin(self.stimulus_timecourse*2*np.pi*self.stimmaker_frequency - np.pi*self.stimmaker_phase)) + self.stimmaker_offset)
        self.stimulus_magnitude[0] = 0.

class Output():
    """Object holding image for now, later might carry more info"""
    
    def __init__(self):
        image = np.zeros((256,256, 1, 1))
        #*span is length/width/height
        xspan = 1.
        yspan = 1.
        zspan = 1.
        #*dim is dimension per unit length (ie 2mm per pixel)
        xdim = 1.
        ydim = 1.
        zdim = 1.
        units = r'$\mathrm{\mu m}$'

class AcquisitionThread(Thread):
    def run(self):
        with self.lock:
            image=self.work()
            self.output.image = image
            #self.output.vmax = image.max()
            #self.output.has_image=True
            #self.image_show(image)
                
    def work(self):
        if (self.workfile == ''):
            self.workfile = 'mfs-cylinder-s-%d.txt' % (self.geom.spacing)
        if (self.workfile[-4:] != '.txt'):
            self.workfile = self.workfile + '.txt'
        self.gensim()
        cmd = '../difsim %s' % (self.workfile)
        call(cmd, shell=True)
        try:
            image = np.genfromtxt(self.geom.ncCurrent_output_file)
        except IOError:
            image = np.zeros((256,256))
        image = self.fake_4d(image)
        return(image)

    def workold(self):
        if (self.workfile == ''):
            self.workfile = 'mfs-cylinder-s-%d.txt' % (self.geom.spacing)
        if (self.workfile[-4:] != '.txt'):
            self.workfile = self.workfile + '.txt'
        self.gensim()
        cmd = '../difsim %s' % (self.workfile)
        p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
        stdout, stderr = p.communicate()
        if not "not found" in stderr:
            out = stdout.split("\n")
            image = np.array([o.split(" ")[:-1] for o in out[3:259]],dtype=float)
        else:
            print stderr.split("\n")["not found" in stderr.split("\n") ]
            image = np.zeros((256,256)) 
        image = self.fake_4d(image)
        return(image)

    #Takes cylinder data and makes copies in z
    def fake_4d(self, image):
        outimage = np.zeros((image.shape[0], image.shape[1],7))
        for i in range(7):
            outimage[:,:,i] = image
        outimage = np.reshape(outimage, (outimage.shape+(1,)))
        return outimage

    def gensim(self):
        out = open(self.workfile, 'w')
        result = []
        for curclass in [self.geom, self.pulsepanel, self.simulation]:
            if (curclass == self.pulsepanel):
                names = self.pulsepanel.savedtraits
            else:
                names = curclass.trait_names(type=not_event)
            names.sort()
            for name in names:
                try:
                    value = getattr( curclass, name )
                except:
                    value = '<undefined>'
                desc = curclass.base_trait( name ).desc
                lname = '   ' + (name.replace('_',' ') + ' =')
                if name == 'q_value':
                    lname = '# ' + 'q = '
                    result.insert(0, '%s %s # %s' % ( lname, value, 'mm^(-1)') )
                elif name == 'pulsetype':
                    lname = '   pulse ='
                    if value == 'gradient echo':
                       value = 0 
                    if value == 'spin echo':
                       value = 1 
                    if value == 'DPFG':
                       value = 2 
                    result.append( '%s %s # %s' % ( lname, value, desc ) )
                else:
                    result.append( '%s %s # %s' % ( lname, value, desc ) )

        result.append( '%s %s # %s' % ( '   ncCurrent_strength_file =', getattr( self.stimpanel, 'ncCurrent_strength_file'), self.stimpanel.base_trait('ncCurrent_strength_file').desc))
            
        out.write('\n'.join(result))

class _SimpleMPLEditor(Editor):
    scrollable=True
    
    def init(self, parent):
        self.control = self._create_canvas(parent)
        self.set_tooltip()

    def update_editor(self):
        pass

    def _create_canvas(self, parent):
        if ETSConfig.toolkit == 'wx':
            panel = wx.Panel(parent, -1, style=wx.CLIP_CHILDREN)
            sizer = wx.BoxSizer(wx.VERTICAL)
            panel.SetSizer(sizer)
            mpl_control = FigureCanvas(panel, -1, self.value)
            sizer.Add(mpl_control, 1, wx.EXPAND)#wx.LEFT | wx.TOP | wx.GROW)
            self.value.canvas.SetMinSize((10,10))
            return panel
             
        if ETSConfig.toolkit == 'qt4':
            frame = QtGui.QWidget()
            mpl_canvas = FigureCanvas(self.value)
            mpl_canvas.setParent(frame)

            vbox = QtGui.QVBoxLayout()
            vbox.addWidget(mpl_canvas)
            #vbox.addWidget(mpl_toolbar)
            frame.setLayout(vbox)
            return frame 
    
class SimpleMPLEditor(BasicEditorFactory):
    klass = _SimpleMPLEditor

class PlotPanel(HasTraits):
    plotfig = Figure()
    #plotfig.add_axes([0.15, 0.15, 0.8, 0.8])

    traits_view = View(Group(Item('plotfig', editor=SimpleMPLEditor(), dock='vertical', show_label=False, height = 600, width=900)))

    def plot_data(self, data):
        self.plotfig.axes[0].clear()
        self.plotfig.axes[0].plot(np.arange(len(data)), data)
        if ETSConfig.toolkit == 'wx':
            wx.CallAfter(self.plotfig.canvas.draw)
        else:
            self.plotfig.canvas.draw()

    def plot_array(self, data):
        self.plotfig.clf()
        gs = gridspec.GridSpec(data.shape[0], data.shape[1])
        #ylimlow = np.min(data[data > 0.1])
        #ylimhigh = np.max(data)
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                print "Subplot %d %d" % (i,j)
                ax = self.plotfig.add_subplot(gs[i,j])
                #This doesn't really work, need to find way to plot on standard axes
                #if (i == 0) and (j == 0):
                #    ax = self.plotfig.add_subplot(gs[i,j])
                #    ax.set_ylim(ylimlow, ylimhigh)
                #    ax1 = ax
                #else:
                #    ax = self.plotfig.add_subplot(gs[i,j], sharex=ax1, sharey=ax1)
                ax.plot(np.arange(len(data[i,j,:])), data[i,j,:])
                setp(ax.get_xticklabels(), visible=False)
                setp(ax.get_yticklabels(), visible=False)
        if ETSConfig.toolkit == 'wx':
            wx.CallAfter(self.plotfig.canvas.draw)
        else:
            self.plotfig.canvas.draw()
        print "Subplot done"
    
        
class SelectRectangle():
    def __init__(self, ax):
        self.ax = ax
        self.rect = Rectangle((0,0), 0, 0, fill=False, color='yellow', lw = 2.)
        self.x0 = 0.
        self.y0 = 0.
        self.x1 = 0.
        self.y1 = 0.
        self.is_pressed = False
        self.ax.add_patch(self.rect)
        self.ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.ax.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.ax.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_press(self, event):
        if event.inaxes == self.ax:
            self.x0 = event.xdata
            self.y0 = event.ydata
            self.is_pressed = True
    
    def on_release(self, event):
        if event.inaxes == self.ax:
            self.is_pressed = False
            self.x1 = event.xdata
            self.y1 = event.ydata
            if (self.x0 == self.x1) and (self.y0 == self.y1):
                self.rect.set_width(2)    
                self.rect.set_height(2)
                self.rect.set_xy((self.x0, self.y0))
            else:
                self.rect.set_width(self.x1 - self.x0)
                self.rect.set_height(self.y1 - self.y0)
                self.rect.set_xy((self.x0, self.y0))
            self.ax.figure.canvas.draw()

    def on_motion(self, event):
        if (event.inaxes == self.ax) and (self.is_pressed == True):
            self.x1 = event.xdata
            self.y1 = event.ydata
            self.rect.set_width(self.x1 - self.x0)
            self.rect.set_height(self.y1 - self.y0)
            self.rect.set_xy((self.x0, self.y0))
            self.ax.figure.canvas.draw()

class VisualizePanel(HasTraits):
    """ Panel for handling all figure functions """
    #Assume image has shape [sagittal, coronal, transverse, time]
    tfigure = Figure()
    sfigure = Figure()
    cfigure = Figure()

    tslice = Int(0, label="transverse slice select", desc="transverse slice to visualize") 
    maxtslice = Int(0, label="max transverse slice", desc="max transverse slice index")
    sslice = Int(0, label="sagittal slice select", desc="sagittal slice to visualize") 
    maxsslice = Int(0, label="max sagittal slice", desc="max sagittal slice index")
    cslice = Int(0, label="coronal slice select", desc="coronal slice to visualize") 
    maxcslice = Int(0, label="max coronal slice", desc="max coronal slice index")

    curtime = Int(0, label="timepoint", desc="time point to visualize")
    maxtime = Int(0, label="max time", desc="max time point index")

    scene = Instance(MlabSceneModel, ())
    engine_view = Instance(EngineView)

    plotpanel = Instance(PlotPanel, ())

    plot_tline = Button("Plot time course")
    plot_sline = Button("Plot time course")
    plot_cline = Button("Plot time course")
    tx, ty, sy, sz, cx, cz  = (0.,)*6
    
    #currect = matplotlib.patches.Rectangle((0,0), 1,1)
    #selectrect = Instance(Annotate, (self.))
    
    view = View(VGroup(
                 Item('curtime', editor=RangeEditor(high_name='maxtime', format='%4d',is_float=False),show_label=True, springy=True),
                 Group(
                 HGroup(VGroup(Item('sfigure',  editor=SimpleMPLEditor(), dock='vertical', show_label=False,),
                               HGroup(Item('sslice', editor=RangeEditor(high_name='maxsslice', format='%2d',is_float=False),show_label=False),
                                      Item('plot_sline', show_label=False),),),
                        VGroup(Item('tfigure',  editor=SimpleMPLEditor(), dock='vertical', show_label=False,),
                               HGroup(Item('tslice', editor=RangeEditor(high_name='maxtslice', format='%2d',is_float=False),show_label=False),
                                      Item('plot_tline', show_label=False),),),),
                 HGroup(Item('scene', editor=SceneEditor(scene_class=MayaviScene), show_label=False,springy=True),
                        VGroup(Item('cfigure',  editor=SimpleMPLEditor(), dock='vertical', show_label=False,springy=True),
                               HGroup(Item('cslice', editor=RangeEditor(high_name='maxcslice', format='%2d',is_float=False),show_label=False),
                                      Item('plot_cline', show_label=False),),),)
                 ,columns=2, orientation='vertical')))

    def __init__(self, output, **traits):
        HasTraits.__init__(self, **traits)
        self.engine_view = EngineView(engine=self.scene.engine)
        self.output = output

    def onclick(self, event):
        if event.inaxes == self.tfigure.axes[0]:
            self.tx = np.floor((event.xdata/self.output.xdim)+np.ceil(self.maxsslice/2.))
            self.ty = np.floor((event.ydata/self.output.ydim)+np.ceil(self.maxcslice/2.))
            del self.tfigure.axes[0].collections[:]
            self.tfigure.axes[0].autoscale(enable=False)
            self.tfigure.axes[0].scatter(event.xdata, event.ydata, c='yellow', marker='s')
            if ETSConfig.toolkit == 'wx':
                wx.CallAfter(self.tfigure.canvas.draw)
            else:
                self.tfigure.canvas.draw()
        if event.inaxes == self.sfigure.axes[0]:
            self.sy = np.floor((event.xdata/self.output.ydim)+np.ceil(self.maxcslice/2.))
            self.sz = np.floor((event.ydata/self.output.zdim)+np.ceil(self.maxtslice/2.))
            del self.sfigure.axes[0].collections[:]
            self.sfigure.axes[0].autoscale(enable=False)
            self.sfigure.axes[0].scatter(event.xdata, event.ydata, c='yellow', marker='s')
            if ETSConfig.toolkit == 'wx':
                wx.CallAfter(self.sfigure.canvas.draw)
            else:
                self.sfigure.canvas.draw()
        if event.inaxes == self.cfigure.axes[0]:
            self.cx = np.floor((event.xdata/self.output.xdim)+np.ceil(self.maxsslice/2.))
            self.cz = np.floor((event.ydata/self.output.zdim)+np.ceil(self.maxtslice/2.))
            del self.cfigure.axes[0].collections[:]
            self.cfigure.axes[0].autoscale(enable=False)
            self.cfigure.axes[0].scatter(event.xdata, event.ydata, c='yellow', marker='s')
            if ETSConfig.toolkit == 'wx':
                wx.CallAfter(self.cfigure.canvas.draw)
            else:
                self.cfigure.canvas.draw()

    def update_mayavi_scene(self, output):
        self.scene.mlab.clf(figure=self.scene.mayavi_scene)
        e = self.scene.engine
        curscene = self.scene 
        curscene.scene.disable_render = True
        s = ArraySource()
        #image3d = self.fake_3d(output.image)
        image3d = output.image
        s.scalar_data = image3d[:,:,:,0]
        s.spacing = (output.xspan/256., output.yspan/256., 1.)
        e.add_source(s, scene = curscene)
        e.add_module(Outline())
        iso = IsoSurface()
        iso.contour.auto_contours=True
        iso.contour.number_of_contours=5
        iso.actor.property.opacity = 0.3
        e.add_module(iso)
        cp = ScalarCutPlane()
        e.add_module(cp)
        cp.implicit_plane.normal = 0,0,1
        curscene.scene.disable_render = False

    #Assume image has shape [sagittal, coronal, transverse, time]
    def plot_all_planes(self):
        self.plot_transverse()
        self.plot_coronal()
        self.plot_sagittal()
        #cidt = self.tfigure.canvas.mpl_connect('button_press_event', self.onclick)
        #cids = self.sfigure.canvas.mpl_connect('button_press_event', self.onclick)
        #cidc = self.cfigure.canvas.mpl_connect('button_press_event', self.onclick)

    def plot_transverse(self):
        self.tfigure.axes[0].clear()
        dx = self.output.xspan/2.
        dy = self.output.yspan/2.
        im = self.tfigure.axes[0].imshow(self.output.image[:,:,self.tslice,self.curtime].T, extent=[-dx,dx,-dy,dy], interpolation='none',  cmap=cm.Greys_r, origin='lower', picker=True)
        self.tfigure.axes[0].get_xaxis().set_major_locator(MaxNLocator(6, integer=True))
        self.tfigure.axes[0].get_yaxis().set_major_locator(MaxNLocator(6, integer=True))
        #self.tfigure.axes[0].invert_xaxis()
        #self.tfigure.axes[0].set_xticklabels(self.tfigure.axes[0].get_xticks()[::-1])
        #self.tfigure.axes[0].invert_yaxis()
        #self.tfigure.axes[0].set_yticklabels(self.tfigure.axes[0].get_yticks()[::-1])
        self.tfigure.axes[0].set_xlabel('X [%s]' % (self.output.units))
        self.tfigure.axes[0].set_ylabel('Y [%s]' % (self.output.units))
        if ETSConfig.toolkit == 'wx':
            wx.CallAfter(self.tfigure.canvas.draw)
        else:
            self.tfigure.canvas.draw()
        self.tselectrect = SelectRectangle(self.tfigure.axes[0])

    def plot_sagittal(self):
        self.sfigure.axes[0].clear()
        dy = self.output.yspan/2.
        dz = self.output.zspan/2.
        im = self.sfigure.axes[0].imshow(self.output.image[self.sslice,:,:,self.curtime].T, extent=[-dy,dy,-dz,dz], interpolation='none', cmap=cm.Greys_r, origin='lower')
        self.sfigure.axes[0].get_xaxis().set_major_locator(MaxNLocator(6, integer=True))
        self.sfigure.axes[0].get_yaxis().set_major_locator(MaxNLocator(6, integer=True))
        self.sfigure.axes[0].set_xlabel('Y [%s]' % (self.output.units))
        self.sfigure.axes[0].set_ylabel('Z [%s]' % (self.output.units))
        if ETSConfig.toolkit == 'wx':
            wx.CallAfter(self.sfigure.canvas.draw)
        else:
            self.sfigure.canvas.draw()
        self.sselectrect = SelectRectangle(self.sfigure.axes[0])

    def plot_coronal(self):
        self.cfigure.axes[0].clear()
        dx = self.output.xspan/2.
        dz = self.output.zspan/2.
        im = self.cfigure.axes[0].imshow(self.output.image[:,self.cslice,:,self.curtime].T, extent=[-dx,dx,-dz,dz], interpolation='none', cmap=cm.Greys_r, origin='lower')
        self.cfigure.axes[0].get_xaxis().set_major_locator(MaxNLocator(6, integer=True))
        self.cfigure.axes[0].get_yaxis().set_major_locator(MaxNLocator(6, integer=True))
        self.cfigure.axes[0].set_xlabel('X [%s]' % (self.output.units))
        self.cfigure.axes[0].set_ylabel('Z [%s]' % (self.output.units))
        if ETSConfig.toolkit == 'wx':
            wx.CallAfter(self.cfigure.canvas.draw)
        else:
            self.cfigure.canvas.draw()
        self.cselectrect = SelectRectangle(self.cfigure.axes[0])

    def update_plots(self):
        self.tfigure.axes[0].images[0].set_data(self.output.image[:,:,self.tslice,self.curtime].T)  
        self.sfigure.axes[0].images[0].set_data(self.output.image[self.sslice,:,:,self.curtime].T)  
        self.cfigure.axes[0].images[0].set_data(self.output.image[:,self.cslice,:,self.curtime].T)  
        if ETSConfig.toolkit == 'wx':
            wx.CallAfter(self.sfigure.canvas.draw)
            wx.CallAfter(self.tfigure.canvas.draw)
            wx.CallAfter(self.cfigure.canvas.draw)
        else:
            self.sfigure.canvas.draw()
            self.tfigure.canvas.draw()
            self.cfigure.canvas.draw()

    def _tslice_changed(self):
        if (self.tslice >= self.maxtslice):
            pass
        else:
            self.plot_transverse()
    
    def _cslice_changed(self):
        if (self.cslice >= self.maxcslice):
            pass
        else:
            self.plot_coronal()
    
    def _sslice_changed(self):
        if (self.sslice >= self.maxsslice):
            pass
        else:
            self.plot_sagittal()

    def _curtime_changed(self):
        if (self.curtime >= self.maxtime):
            pass
        else:
            self.plot_all_planes()

    def _plot_tline_fired(self):
        self.plotpanel.configure_traits()
        imdim = self.tfigure.axes[0].images[0].get_extent()
        xextent = imdim[1] - imdim[0]
        yextent = imdim[3] - imdim[2]
        dx0 = int((self.tselectrect.x0 + xextent/2.)/self.output.xdim)
        dy0 = int((self.tselectrect.y0 + yextent/2.)/self.output.ydim)
        dx1 = int((self.tselectrect.x1 + xextent/2.)/self.output.xdim)
        dy1 = int((self.tselectrect.y1 + yextent/2.)/self.output.ydim)
        if (dx1 < dx0):
            dx0, dx1 = dx1, dx0
        if (dy1 < dy0):
            dy0, dy1 = dy1, dy0
        if (dx0 == dx1) and (dy0 == dy1):
            outdat = self.output.image[dx0,dy0, self.tslice,:].reshape((1,1,self.maxtime+1))
        elif (dx0 == dx1):
            outdat = self.output.image[dx0,dy0:dy1, self.tslice,:].reshape((dy1-dy0,1,self.maxtime+1))
        elif (dy0 == dy1):
            outdat = self.output.image[dx0:dx1,dy0, self.tslice,:].reshape((1,dx1-dx0,self.maxtime+1))
        else:
            outdat = swapaxes(self.output.image[dx0:dx1,dy0:dy1, self.tslice,:], 0, 1)[::-1,:,:]
        self.plotpanel.plot_array(outdat)

    def _plot_sline_fired(self):
        self.plotpanel.configure_traits()
        imdim = self.sfigure.axes[0].images[0].get_extent()
        xextent = imdim[1] - imdim[0]
        yextent = imdim[3] - imdim[2]
        dx0 = int((self.sselectrect.x0 + xextent/2.)/self.output.ydim)
        dy0 = int((self.sselectrect.y0 + yextent/2.)/self.output.zdim)
        dx1 = int((self.sselectrect.x1 + xextent/2.)/self.output.ydim)
        dy1 = int((self.sselectrect.y1 + yextent/2.)/self.output.zdim)
        if (dx1 < dx0):
            dx0, dx1 = dx1, dx0
        if (dy1 < dy0):
            dy0, dy1 = dy1, dy0
        print("Select coords 0: %d %d 1: %d %d" % (self.sselectrect.x0, self.sselectrect.y0, self.sselectrect.x1, self.sselectrect.y1) )

        if (dx0 == dx1) and (dy0 == dy1):
            outdat = self.output.image[self.sslice,dx0,dy0,:].reshape((1,1,self.maxtime+1))
        elif (dx0 == dx1):
            outdat = self.output.image[self.sslice,dx0,dy0:dy1,:].reshape((dy1-dy0,1,self.maxtime+1))
        elif (dy0 == dy1):
            outdat = self.output.image[self.sslice,dx0:dx1,dy0,:].reshape((1,dx1-dx0,self.maxtime+1))
        else:
            outdat = swapaxes(self.output.image[self.sslice,dx0:dx1,dy0:dy1,:], 0, 1)[::-1,:,:]
        self.plotpanel.plot_array(outdat)

    def _plot_cline_fired(self):
        self.plotpanel.configure_traits()
        imdim = self.cfigure.axes[0].images[0].get_extent()
        xextent = imdim[1] - imdim[0]
        yextent = imdim[3] - imdim[2]
        dx0 = int((self.cselectrect.x0 + xextent/2.)/self.output.xdim)
        dy0 = int((self.cselectrect.y0 + yextent/2.)/self.output.zdim)
        dx1 = int((self.cselectrect.x1 + xextent/2.)/self.output.xdim)
        dy1 = int((self.cselectrect.y1 + yextent/2.)/self.output.zdim)
        if (dx1 < dx0):
            dx0, dx1 = dx1, dx0
        if (dy1 < dy0):
            dy0, dy1 = dy1, dy0

        if (dx0 == dx1) and (dy0 == dy1):
            outdat = self.output.image[dx0,self.cslice,dy0,:].reshape((1,1,self.maxtime+1))
        elif (dx0 == dx1):
            outdat = self.output.image[dx0,self.cslice,dy0:dy1,:].reshape((dy1-dy0,1,self.maxtime+1))
        elif (dy0 == dy1):
            outdat = self.output.image[dx0:dx1,self.cslice,dy0,:].reshape((1,dx1-dx0,self.maxtime+1))
        else:
            outdat = swapaxes(self.output.image[dx0:dx1,self.cslice,dy0:dy1,:], 0, 1)[::-1,:,:]
        self.plotpanel.plot_array(outdat)
        
class StimulusPanel(HasTraits):
    ncCurrent_strength_file = Str("./sinusoidal_stimulus.txt", label="ncCurrent strength file", desc="file containing stimulus strength over time")

    stimulus = Instance(Stimulus, ())
    load_stimulus = Button("Load Stimulus File")
    stim_sinusoid = Button("Sinusoidal")
    stim_sawtooth = Button("Sawtooth")
    stim_square = Button("Square")
    stim_configure_params = Button("Configure Stimulus Parameters")

    stimfig = Figure()
    stimfig.add_axes([0.15, 0.15, 0.8, 0.8])

    view = View(Group(
                    HGroup(
                        Item('ncCurrent_strength_file', springy=True, show_label=False,),
                        Item('load_stimulus', show_label=False),
                        ),
                    Item('stimulus', style='custom', show_label=False,),
                    HGroup(
                        Item('stim_sinusoid', show_label=False,),  
                        Item('stim_sawtooth', show_label=False,),  
                        Item('stim_square', show_label=False,),
                        ),
                    Item('stimfig', editor=SimpleMPLEditor(),  show_label=False,),
                    springy=True),)


    def _load_stimulus_fired(self):
        stimfile = np.genfromtxt(self.ncCurrent_strength_file)
        self.stimulus.stimulus_timecourse = stimfile[:,0]
        self.stimulus.stimulus_magnitude = stimfile[:,1]
        self.plot_stimulus()

    def plot_stimulus(self):
        self.stimfig.axes[0].clear()
        self.stimfig.axes[0].images=[]
        self.stimfig.axes[0].plot(self.stimulus.stimulus_timecourse, self.stimulus.stimulus_magnitude)
        if ETSConfig.toolkit == 'wx':
            wx.CallAfter(self.stimfig.canvas.draw)
        else:
            self.stimfig.canvas.draw()
    
    def _stim_sinusoid_fired(self):
        self.stimulus.sinusoidal_stimulus()
        self.plot_stimulus()
    def _stim_sawtooth_fired(self):
        self.stimulus.sawtooth_stimulus()
        self.plot_stimulus()
    def _stim_square_fired(self):
        self.stimulus.square_stimulus()
        self.plot_stimulus()
        self.stimfig.axes[0].set_ylim(-0.5, self.stimulus.stimmaker_max_magnitude + 1.)

class BlochSimParams(HasTraits):
    frame = Enum('PM', 'FM')
    t1 = Range(0, 2000, 250, label="t1", desc="T1 (ms) {0,2000}")
    t2 = Range(0, 500, 100, label="t2", desc="T2 (ms) {0,500}")
    b1 = Range(1, 5, 1, label="b1", desc="0.1456 x b1 (Gauss), {1,5}")
    pulsedur = Range(1, 100, 10, label="pulsedur", desc="pulse duration (ms) {1,100}")
    offsetfreq = Range(0, 10, 1, label="offsetfreq", desc="offset frequency (kHz)")
    rftype = Enum(['0 inversion','1 excitation','2 refocusing'],label="rf type", desc="rf type (0) Inversion, (1) excitation, (2) refocusing {0,2}")
    #rftype = Range(low=int(0),high=int(2),value=int(0),label="rf type", desc="rf type (0) Inversion, (1) excitation, (2) refocusing {0,2}", is_float=False)
    beta = Range(1, 1000, 400, label="beta", desc="beta (s^-1) {1,1000}")
    mu = Range(1, 20, 10, label="mu", desc="mu {1,20}")
    vrgon = Bool(0, label="vrgon", desc="use variable rate gradient? {0,1}")
    opgp = Bool(0, label="opgp", desc="composite gp pulse? {0,1}")
    gph = Range(0.,1., 0.5, label="gph", desc="gph (fraction of PI) {0,1}")

    traits_view = View(VGroup(Item('frame'), HGroup(Item('t1'), Item('t2')), Item('b1'), Item('pulsedur'), 
                              Item('offsetfreq'), Item('rftype')#, editor=RangeEditor(mode='enum',format='%d'))
                              , Item('beta'), Item('mu')), HGroup(Item('vrgon'), Item('opgp')), Item('gph', visible_when='opgp'))

class PulsePanel(HasTraits):
    pulsetype = Enum('DPFG', 'gradient echo', 'spin echo')

    ramp = Int(10, label="ramp", desc="ramp in microseconds (for both leading and trailing edges)")
    Delta = Int(200, label="Delta", desc="Delta in microseconds (Delta - delta using standard definitions)")
    delta = Int(100, label="delta", desc="delta in microseconds")
    mixing_time = Int(0, label="mixing_time", desc="mixing time in microseconds (only for DPFG pulses)")
    gradient_strength = Int(1, label="gradient_strength", desc="gradient strength in G/cm")    
    snr = Int(100, label="snr", desc="signal-to-noise ratio (percent)")
    DPFG_angle = Int(0, label="DPFG_angle", desc="DPFG angle (pulse=2)")
    q_value = Int(0, label="q_value", desc="q-value")
    G = Property(depends_on='q_value')
    ntess = Int(0, label="ntess", desc="tesselation number (not used for DPFG pulses)")
    ndifdir = Range(low=6,high=150, label="ndifdirn", desc="number of diffusion directions")
    ncCurrent_strength = Float(1.0, label="ncCurrent_strength", desc="nc current strength")

    savedtraits = ['pulsetype', 'ramp', 'Delta', 'delta', 'mixing_time', 'gradient_strength', 'snr', 'DPFG_angle', 'q_value', 'G', 'ntess', 'ncCurrent_strength']

    blochparams = Instance(BlochSimParams, ())
    
    is_dpfg = Property(depends_on='pulsetype')

    pulsefig = Figure()
    plot_pulse = Button("Plot pulse sequence")
    
    difdirscene = Instance(MlabSceneModel, ())
    difdir_engine_view = Instance(EngineView)
    plot_difdir = Button("Plot diffusion directions") 

    blochfig = Instance(MlabSceneModel, ())
    blochfig_engine_view = Instance(EngineView)
    plot_bloch = Button('Plot Bloch sphere')
    config_bloch = Button('Set Bloch sim parameters')


    traits_view = View(Item('pulsetype', show_label=False),
                       VGrid(Item('Delta'), Item('delta'), Item('ramp'), Item('snr'), Item('q_value'), Item('G'),                               Item('ntess'),
                       VGrid(Item('mixing_time'), Item('DPFG_angle'), visible_when = 'is_dpfg'),),
                       Group(
                         VGroup(
                            Item('plot_pulse', show_label=False,),
                            Item('pulsefig', editor=SimpleMPLEditor(),  show_label=False,springy=True,),label='Pulse Viewer', springy=True),
                         VGroup(
                            HGroup(Item('plot_difdir', show_label=False,),Item('ndifdir'),),
                            Item('difdirscene', editor=SceneEditor(scene_class=MayaviScene), show_label=False),label='Diffusion Directions'),
                         VGroup(
                            HGroup(Item('plot_bloch', show_label=False,), Item('config_bloch', show_label=False,springy=True),),
                            Item('blochfig', editor=SceneEditor(scene_class=MayaviScene), show_label=False),label='Bloch sphere'),
                         layout='tabbed',springy=True
                       ),)


    def __init__(self, **traits):
        HasTraits.__init__(self, **traits)
        self.difdir_engine_view = EngineView(engine=self.difdirscene.engine)

    def _get_is_dpfg(self):
        return self.pulsetype == 'DPFG'
 
    def _get_G(self):
        return 10000.0*math.pi*self.q_value/26752.0
    
    def parsetensor(self):
        dat = np.array(open('tensor.dat').readlines())
        testl = '%s\n' % (self.ndifdir)
        w = int(np.where(dat==testl)[0]+1)
        v = np.array([curd.split(' ') for curd in dat[w:w+self.ndifdir]], dtype=float)
        return v
    
    def _plot_difdir_fired(self):
        self.difdirscene.mlab.clf(figure=self.difdirscene.mayavi_scene)
        e = self.difdirscene.engine
        curscene = self.difdirscene
        curscene.scene.disable_render = True
        dif_vector = self.parsetensor()
        z = np.zeros((self.ndifdir, 3)) 
        origvec = quiver3d(z[:,0], z[:,1], z[:,2], dif_vector[:,0], dif_vector[:,1], dif_vector[:,2], figure=self.difdirscene.mayavi_scene) 
        negvec = quiver3d(z[:,0], z[:,1], z[:,2], -dif_vector[:,0], -dif_vector[:,1], -dif_vector[:,2], color=(0,0,1.), figure=self.difdirscene.mayavi_scene) 
        e.add_module(origvec)
        e.add_module(negvec)
        curscene.scene.disable_render = False
        curscene.scene.parallel_projection = True 

    def _plot_pulse_fired(self):
        if self.pulsetype=='gradient echo' or self.pulsetype=='spin echo':
            self.plot_combo()
        elif self.pulsetype=='DPFG':
            self.plot_dpfg()
        else:
            pass

    def make_ramppulse(self, Gstr, ramptime, ontime):
        return Gstr*np.concatenate((np.linspace(0,1,ramptime), np.ones(ontime), np.linspace(1,0, ramptime)))

    def plot_combo(self):
        rfoffset = 160 + 2*self.ramp
        hrfoffset = rfoffset/2.
        totaltime = 4*self.ramp + 2*self.delta + self.Delta
        afteracqtime = 2*self.ramp+self.delta
        time = np.linspace(-rfoffset,totaltime+afteracqtime, totaltime+rfoffset+afteracqtime+1)
        ramppulse = self.make_ramppulse(1., self.ramp, self.delta)
        
        if self.pulsetype=='gradient echo':
            cGS = np.concatenate((np.zeros(rfoffset+1), ramppulse, np.zeros(self.Delta), -ramppulse, np.zeros(afteracqtime)))
            rfpulse = np.exp(-0.02*(time+rfoffset/2.)**2)*np.cos(0.5*(time+rfoffset/2.))
            sliceselect = np.concatenate((np.zeros(rfoffset/4.-self.ramp), self.make_ramppulse(0.5,self.ramp, rfoffset/2.), np.zeros(rfoffset/4.-self.ramp), -self.make_ramppulse(0.25, self.ramp, self.delta), np.zeros(totaltime+1-2*self.ramp-self.delta), np.zeros(afteracqtime)))
            phaseencoding = np.concatenate((np.zeros(rfoffset+1), np.sin(np.pi*np.linspace(0,1, 2*self.ramp + self.delta)), np.zeros(2*self.ramp+self.delta+self.Delta), np.zeros(afteracqtime)))
            signal = np.concatenate((np.ones(rfoffset+1), np.linspace(1.,0.,2*self.ramp+self.delta)**2, np.zeros(self.Delta), np.linspace(0.,1., 2*self.ramp+self.delta)**2, np.ones(afteracqtime)))
        elif self.pulsetype=='spin echo':
            rfpulse = np.exp(-0.02*(time+rfoffset/2.)**2)*np.cos(0.5*(time+rfoffset/2.)) + 2*np.exp(-0.02*(time-totaltime/2.)**2)*np.cos(0.5*(time-totaltime/2.))
            cGS = np.concatenate((np.zeros(rfoffset+1), ramppulse, np.zeros(self.Delta), ramppulse, np.zeros(afteracqtime)))
            sliceselect = np.concatenate((np.zeros(rfoffset/4.-self.ramp), self.make_ramppulse(0.5,self.ramp, rfoffset/2.), np.zeros(totaltime+rfoffset/4.+1-self.ramp), np.zeros(afteracqtime)))
            phaseencoding = np.concatenate((np.zeros(rfoffset+1), np.sin(np.pi*np.linspace(0,1, 2*self.ramp + self.delta)), np.zeros(2*self.ramp+self.delta+self.Delta), np.zeros(afteracqtime)))
            sigtime = np.linspace(0,hrfoffset+totaltime,hrfoffset + totaltime)
            expsig = np.exp((sigtime-len(sigtime)/2.)/10)*(np.exp(-1j*sigtime))
            subsig = (2./expsig.max())*(expsig + np.exp(-(sigtime-len(sigtime)/2)/10)*np.exp(-1j*sigtime))
            signal = np.concatenate((np.zeros(hrfoffset+1), subsig, subsig[:afteracqtime]))


        ax1 = self.pulsefig.add_subplot(611)
        ax1.clear()
        ax1.plot(time, rfpulse)
        ax1.yaxis.set_major_locator(MaxNLocator(3))
        ax1.set_ylabel('RF pulse')
        setp(ax1.get_xticklabels(), visible=False)

        ax2 = self.pulsefig.add_subplot(612, sharex=ax1)
        ax2.clear()
        ax2.plot(time, sliceselect)
        ax2.set_ylim(-1., 1.)
        ax2.yaxis.set_major_locator(MaxNLocator(3))
        ax2.set_ylabel('Gs')
        setp(ax2.get_xticklabels(), visible=False)

        ax3 = self.pulsefig.add_subplot(613, sharex=ax1)
        ax3.clear()
        ax3.plot(time, cGS)
        ax3.set_ylim(-1.5*self.gradient_strength,1.5*self.gradient_strength)
        ax3.yaxis.set_major_locator(MaxNLocator(3))
        ax3.set_ylabel('Gx')
        setp(ax3.get_xticklabels(), visible=False)

        ax4 = self.pulsefig.add_subplot(614, sharex=ax1)
        ax4.clear()
        ax4.plot(time, 0.5*phaseencoding)
        ax4.set_ylim(-0.2, 1.)
        ax4.set_ylabel('Gy')
        setp(ax4.get_xticklabels(), visible=False)
        ax4.yaxis.set_major_locator(MaxNLocator(3))

        ax5 = self.pulsefig.add_subplot(615, sharex=ax1)
        ax5.clear()
        ax5.plot(time,signal)
        setp(ax5.get_xticklabels(), visible=False)
        if self.pulsetype=='gradient echo':
            ax5.set_ylim(-.2,2.)
        elif self.pulsetype=='spin echo':    
            ax5.set_ylim(-2.,2.)
        ax5.yaxis.set_major_locator(MaxNLocator(3))
        ax5.set_ylabel('signal')

        acqsignal = np.concatenate((np.zeros(rfoffset+totaltime), np.ones(1), np.zeros(afteracqtime)))
        ax6 = self.pulsefig.add_subplot(616, sharex=ax1)
        ax6.clear()
        ax6.plot(time, acqsignal)
        ax6.set_xlim(-rfoffset, totaltime+afteracqtime)
        ax6.set_ylim(-0.2, 1.2)
        ax6.set_yticks([0,0.5,1])
        ax6.xaxis.set_major_locator(MaxNLocator(8, integer=True))
        ax6.set_ylabel('A/D')
        ax6.set_xlabel('time [ms]')
        if ETSConfig.toolkit == 'wx':
            wx.CallAfter(self.pulsefig.canvas.draw)
        else:
            self.pulsefig.canvas.draw()

    def plot_dpfg(self):
        pass

    def _plot_bloch_fired(self): 
        self.plot_sphere()

    def _config_bloch_fired(self):
        self.blochparams.configure_traits()

    def run_bloch_sim(self):
        cmd = './bloch/bloch -c %s %d %d %d %d %d %d %d %d %d %d %d' % (self.blochparams.frame, self.blochparams.t1, self.blochparams.t2, self.blochparams.b1, self.blochparams.pulsedur, self.blochparams.offsetfreq, int(str.split(self.blochparams.rftype,' ')[0]), self.blochparams.beta, self.blochparams.mu, self.blochparams.vrgon, self.blochparams.opgp, self.blochparams.gph)
        #FNULL = open(devnull, 'w')
        #call(cmd, shell=True, stdout=FNULL, stderr=STDOUT)
        call(cmd, shell=True)

    def plot_sphere(self, n=100,rad=1):
        self.blochfig.mlab.clf(figure=self.blochfig.mayavi_scene)
        e = self.blochfig.engine
        curscene = self.blochfig
        curscene.scene.disable_render = True
        self.run_bloch_sim() 
        d=open('./Mxyz.dat').readlines()
        npts = len(d)
        x=map(lambda s: np.array(s.strip().split(' ')).astype(float),d)
        Mxyz = np.swapaxes(np.array(x),0,1)
        u = np.linspace(0, 2*np.pi, n)
        v = np.linspace(0, np.pi, n)
        x = rad*np.outer(np.cos(u), np.sin(v))
        y = rad*np.outer(np.sin(u), np.sin(v))
        z = rad*np.outer(np.ones(np.size(u)), np.cos(v))

        sphere = mesh(x,y,z, color=(1,0.8,0), opacity=0.5, figure=self.blochfig.mayavi_scene)
        e.add_module(sphere)
        traj = plot3d(Mxyz[0], Mxyz[1], Mxyz[2], tube_radius=None, line_width=2., color = (0.,0.,1.))
        e.add_module(traj)
        curscene.scene.disable_render = False
        curscene.scene.parallel_projection = True 

class GeometryPanel(HasTraits):
    """Object containing description of geometry"""
    spacing = Float(2., label="spacing", desc="spacing")
    periodic = Int(1, label="periodic", desc="periodic boundary condition (1=on, 0=off)")
    voxel_size = Int(80, label="voxel_size", desc="voxel size in microns, needs to be bigger than largest dim")
    x_dim = Property(depends_on=['spacing','cylinder_diameter'], label="x_dim", desc="x dimension in um") 
    y_dim = Property(depends_on=['spacing','cylinder_diameter'], label="y_dim", desc="y dimension in um") 
    z_dim = Property(depends_on=['spacing','cylinder_diameter'], label="z_dim", desc="z dimension in um") 
    cylinder_array = Int(1, label="cylinder_array", desc="cylinder array")
    cylinder_diameter = Float(7.992, label="cylinder_diameter", desc="cylinder diameter")
    cylinder_radius = Float(4.0, label="cylinder_radius", desc="cylinder radius")
    permeable = Float(1, label="permeable", desc="permeable")
    structure_file = Str("./cylinder-8.mdl", label="structure_file", desc="mcell file describing geometry")
    ncCurrent_strength = Float(1.0, label="ncCurrent_strength", desc="nc current strength")
    
    plot_mesh = Button('Plot mesh')

    def _get_x_dim(self):
        return self.cylinder_diameter + self.spacing
    def _get_y_dim(self):
        return np.sqrt(3)*(self.cylinder_diameter + self.spacing)
    def _get_z_dim(self):
        return 7.
    
    geomscene = Instance(MlabSceneModel, ())
    geom_engine_view = Instance(EngineView)

    traits_view = View(Group(Item('plot_mesh', show_label=False),
                      Item('geomscene', editor=SceneEditor(scene_class=MayaviScene), show_label=False)))

    def __init__(self, **traits):
        HasTraits.__init__(self, **traits)
        self.geom_engine_view = EngineView(engine=self.geomscene.engine)
   
    def _plot_mesh_fired(self):
        self.plot_geom()

    def plot_geom(self):
        self.geomscene.mlab.clf(figure=self.geomscene.mayavi_scene)
        e = self.geomscene.engine
        curscene = self.geomscene
        curscene.scene.disable_render = True

        filename = './cylinders.stl'
        stlr = tvtk.STLReader()
        stlr.file_name = filename
        stlr.update()
        stld = stlr.output

        surf = pipeline.surface(stld, opacity=1, color=(0, 1, 1), figure=self.geomscene.mayavi_scene)
        edge = pipeline.extract_edges(surf, figure=self.geomscene.mayavi_scene)
        edge_surf = pipeline.surface(edge, opacity=.1, color=(1,0,0), figure=self.geomscene.mayavi_scene)
        e.add_module(edge_surf)
        curscene.scene.disable_render = False
        curscene.scene.parallel_projection = True 

class ControlPanel(HasTraits):
    totalanimtime = Float(0.2, label="total animation time", desc="total time over which to animate (s)")

    #Objects holding all parameters
    geom = Instance(Geometry, ())
    simulation = Instance(Simulation, ())


    output = Instance(Output, ())
    eventlock = Lock() 
    
    start_simulation = Button("Start computation")
    start_stop_animation = Button("Animate")

    load_data = Button("Load data")

    workfile = Str('difsim_config', label='DifSim config', desc='file to save DifSim model')

    acquisition_thread = Instance(AcquisitionThread)
    animation_thread = Instance(Thread)
    animation_wants_abort = False
    vpanel = Instance(VisualizePanel, (output))
    stimpanel = Instance(StimulusPanel, ())

    pulsepanel = Instance(PulsePanel, ())

    geompanel = Instance(GeometryPanel, ())
        
    view = View(Group(
                    HGroup(
                        Item('load_data', show_label=False,),
                        Item('start_stop_animation', show_label=False,),
                        Item('totalanimtime', show_label=False,),),
                  Group(
                      Item('stimpanel', style='custom', show_label=False, label='Stimulus',),
                      Item('pulsepanel', style='custom', show_label=False, label='Pulse Sequence'),
                      VGroup(Item('simulation', style='custom', show_label=False,),
                             Item('start_simulation', show_label=False,),
                             Item('workfile',show_label=False),label="Simulation"),
                      Item('geom', style='custom', show_label=False,),
                      Item('geompanel', style='custom', show_label=False,),
                    dock="tab", layout='tabbed',springy=True),
                  ),#layout=''),
               )
         
    def _start_simulation_fired(self):
        self.acquisition_thread = AcquisitionThread()
        self.acquisition_thread.workfile = self.workfile
        self.acquisition_thread.geom= self.geom
        self.acquisition_thread.simulation = self.simulation
        self.pulsepanel.pulsetype = 'DPFG'
        self.acquisition_thread.pulsepanel = self.pulsepanel
        self.acquisition_thread.stimpanel = self.stimpanel
        self.acquisition_thread.output = self.output

        self.acquisition_thread.output.xspan = self.geom.x_dim
        self.acquisition_thread.output.yspan = self.geom.y_dim
        self.acquisition_thread.output.zspan = self.geom.z_dim
        self.acquisition_thread.output.xdim = 1.
        self.acquisition_thread.output.ydim = 1.
        self.acquisition_thread.output.zdim = 1.

        #self.acquisition_thread.output.xdir = -1.
        self.acquisition_thread.output.units = r'$\mathrm{\mu m}$'

        self.acquisition_thread.lock = self.eventlock
        self.acquisition_thread.start()

        with self.eventlock:
            self.vpanel.maxsslice = self.output.image.shape[0]-1
            self.vpanel.maxcslice = self.output.image.shape[1]-1
            self.vpanel.maxtslice = self.output.image.shape[2]-1
            self.vpanel.plot_all_planes()
            self.vpanel.update_mayavi_scene(self.acquisition_thread.output)
    
    def animate_plot_helper(self, curt):
        self.vpanel.curtime = curt
        self.vpanel.update_plots()

    def _start_stop_animation_fired(self):
        if self.animation_thread and self.animation_thread.isAlive():
            self.animation_wants_abort = True
            return
        self.animation_wants_abort=False
        self.animation_thread = Thread(target=self.animate_plot)
        self.animation_thread.daemon=True
        self.animation_thread.start()

    def animate_plot(self):
        #curim = self.acquisition_thread.output.image
        if (self.vpanel.curtime == self.vpanel.maxtime):
            self.vpanel.curtime = 0
        for curt in range(self.vpanel.curtime, self.vpanel.maxtime):
            if (self.animation_wants_abort):
                break
            else:
                self.animate_plot_helper(curt)
                sleep(self.totalanimtime)
        #for st in self.stimpanel.stimulus.stimulus_magnitude:
        #    if (self.animation_wants_abort):
        #       break 
        #    else:
        #       self.animate_plot_helper(st*curim)
        #       sleep(sleeptime)
            
    def _load_data_fired(self):
        file_name = ''
        if ETSConfig.toolkit == 'wx':
            dlg = wx.FileDialog(None, "Choose a NFITI file", os.getcwd(), "", "*.*", wx.OPEN)
            if dlg.ShowModal() == wx.ID_OK:
                file_name = dlg.GetPath()
        else: 
            file_name = QtGui.QFileDialog.getOpenFileName(None, 'Choose a NIFTI file', os.getcwd())[0]

        if file_name == '':
            pass 
        else:
            #Need to make sure data is in format [x,y,z,time] with x left to right
            niftidat = nipy.load_image(file_name)
            if len(niftidat.shape) < 4:
                self.output.image = np.reshape(niftidat.get_data(), (niftidat.shape+(1,)))
                self.vpanel.maxtime = 0
            else:
                self.output.image = niftidat.get_data()
                self.vpanel.maxtime   = niftidat.shape[3]-1
            niftidims = niftidat.header.get_zooms()
            self.output.units = 'mm'
            self.output.xspan = niftidims[0]*niftidat.shape[0]
            self.output.yspan = niftidims[1]*niftidat.shape[1]
            self.output.zspan = niftidims[2]*niftidat.shape[2]
            self.output.xdim = niftidims[0]
            self.output.ydim = niftidims[1]
            self.output.zdim = niftidims[2]
            # xdir determines frame of reference for plotting
            #Shouldn't need to worry about xdir, put data in right shape in first place
            if (np.sign(niftidat.header.get_sform()[0][0]) > 0): #LAS, needs to be fixed to RAS
                self.output.image = self.output.image[::-1, :, :, :] 
            self.vpanel.maxsslice = niftidat.shape[0]-1
            self.vpanel.maxcslice = niftidat.shape[1]-1
            self.vpanel.maxtslice = niftidat.shape[2]-1
            self.vpanel.plot_all_planes()

class MainWindowHandler(Handler):
    def close(self, info, is_OK):
        if ( info.object.panel.acquisition_thread 
                        and info.object.panel.acquisition_thread.isAlive() ):
            info.object.panel.animation_wants_abort = True
            while info.object.panel.acquisition_thread.isAlive():
                sleep(0.5)
            #wx.Yield()
        return True

class MainWindow(HasTraits):
    """ The main window, here go the instructions to create and destroy the application."""
    output = Instance(Output, ()) 
    vpanel = Instance(VisualizePanel, (output))
    panel = Instance(ControlPanel)
    

    def _vpanel_default(self):
        vpanel = VisualizePanel(self.output)
        vpanel.tfigure.add_axes([0.17, 0.14, 0.78, 0.78])
        vpanel.cfigure.add_axes([0.17, 0.14, 0.78, 0.78])
        vpanel.sfigure.add_axes([0.17, 0.14, 0.78, 0.78])
        
        return vpanel 

    def _panel_default(self):
        return ControlPanel(vpanel=self.vpanel, output=self.output)

    
    if ETSConfig.toolkit == 'wx':
        panelwidth = 450
    else:
        panelwidth = 500

    view = View(HSplit( Item('vpanel', style="custom", springy=True),
                        Item('panel', style="custom", width = panelwidth),
                    show_labels=False, 
                    ),
                resizable=True, 
                height=0.85, width=0.80,
                handler=MainWindowHandler(),
                buttons=NoButtons, title="DifSim")
                
if __name__ == '__main__':
    MainWindow().configure_traits()

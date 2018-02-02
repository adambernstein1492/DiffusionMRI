# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

# <pep8 compliant>

"""
This file contains the classes for Difsim GUI

"""

# blender imports
import bpy
from bpy.props import BoolProperty, CollectionProperty, EnumProperty, \
                      FloatProperty, FloatVectorProperty, IntProperty, \
                      IntVectorProperty, PointerProperty, StringProperty
from bl_operators.presets import AddPresetBase
from bpy.app.handlers import persistent
from bpy.types import Menu

# python imports
import sys
import numpy 
import os
import shutil
import re
import subprocess
import difsim
from cellblender.cellblender_utils import mcell_files_path
import cellblender
import fileinput 
#import mathutils
#from . import mathviz_draw


# We use per module class registration/unregistration
def register():
    bpy.utils.register_module(__name__)


def unregister():
    bpy.utils.unregister_module(__name__)

class DIFSIM_OT_set_difsim_binary(bpy.types.Operator):
    bl_idname = "difsimgui.set_difsim_binary"
    bl_label = "Set Difsim Binary"
    bl_description = "Set Difsim Binary"
    bl_options = {'REGISTER'}

    filepath = bpy.props.StringProperty(subtype='FILE_PATH', default="")

    def execute(self, context):
        context.scene.difsim.difsim_binary = self.filepath
        bpy.ops.difsim.preferences_save(name='Ds')
        return {'FINISHED'}

    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

class DIFSIM_MT_presets(Menu):
    bl_label = "DifSim Blender Presets"
    preset_subdir = "difsim"
    preset_operator = "script.execute_preset"
    draw = Menu.draw_preset

class DIFSIM_OT_save_preferences(AddPresetBase, bpy.types.Operator):
    """Save DifSimBlender Preferences."""

    bl_idname = "difsim.preferences_save"
    bl_label = "Save Preferences"
    preset_menu = "DIFSIM_MT_presets"

    preset_defines = [ "scene = bpy.context.scene"]

    preset_values = [
        "scene.difsim.difsim_binary",
    ]

    preset_subdir = "difsim"

class DIFSIM_OT_update_molecule_list(bpy.types.Operator):
    bl_idname = "difsimgui.update_molecule_list"
    bl_label = "Update Molecule List"
    bl_description = "Update Molecule List"
    bl_options = {'REGISTER'}

    def execute(self, context):
        if context.scene.difsim.molecules.molecule_list:
            context.scene.difsim.molecules.remove_properties(context) 
        for mcell_mol in context.scene.mcell.molecules.molecule_list:
            context.scene.difsim.molecules.add_molecule(self, mcell_mol.name, float(mcell_mol.diffusion_constant.get_expr()))
        return {'FINISHED'}

class DIFSIM_OT_update_bounds(bpy.types.Operator):
    bl_idname = "difsimgui.update_bounds"
    bl_label = "Update Voxel Boundaries from MCell Model"
    bl_description = "Update Voxel Boundaries from MCell Model"
    bl_options = {'REGISTER'}

    def execute(self, context):
        bpy.ops.mcell.auto_generate_boundaries()
        #Must be centered on 0
        context.scene.difsim.xDim = 2*numpy.max((context.scene.mcell.partitions.x_start, context.scene.mcell.partitions.x_end))
        context.scene.difsim.yDim = 2*numpy.max((context.scene.mcell.partitions.y_start, context.scene.mcell.partitions.y_end))
        context.scene.difsim.zDim = 2*numpy.max((context.scene.mcell.partitions.z_start, context.scene.mcell.partitions.z_end))
        return {'FINISHED'}

# Operators can't be callbacks, so we need this function.
def save_preferences(self, context):
    bpy.ops.difsim.preferences_save(name='Ds')
    return None

@persistent
def load_preferences(context):
    """ Load DS preferences using preset on startup. """
    print ( "load post handler: difsim_operators.load_preferences() called" )

    #active_preset = bpy.types.MCELL_MT_presets.bl_label
    #bpy.ops.script.execute_preset(
    #    filepath=preset_path, menu_idname="MCELL_MT_presets")

    # Cant use execute_preset here because of incorrect context
    # Look for existing preset named "Cb" in "cellblender" preset subdir
    # Presets are named with leading capital, but modules are lowercase
    # Therefore the preset "Cb" corresponds to the file "cb.py"
    try:
        preset_path = bpy.utils.preset_find(
            "Ds", 'difsim', display_name=True)
        with open(preset_path, 'r') as preset_file:
            preset_string = preset_file.read()
            exec(preset_string)
    except TypeError:
        print('No preset file found')


class DIFSIM_OT_export_difsim(bpy.types.Operator):
    bl_idname = "difsimgui.export_difsim"
    bl_label = "Output Difsim Model File"
    bl_description = "Output Difsim Model File"
    bl_options = {'REGISTER'}

    def execute(self, context):
        mcell_project_dir = mcell_files_path()
        config_file = os.path.join(mcell_project_dir, 'difsim_params.txt') 
        bpy.ops.mcell.export_project()
        self.finalize_mcellfiles(context)
        context.scene.difsim.export_difsim(context, config_file)
        return {'FINISHED'}

    def finalize_mcellfiles(self, context):
        mcell_project_dir = mcell_files_path()
        base_name = context.scene.mcell.project_settings.base_name
        counter = -1
        for line in fileinput.input(os.path.join(mcell_project_dir, base_name) + '.main.mdl', inplace=True):
            if fileinput.isfirstline(): 
                if (context.scene.difsim.pulsetype == 'DWSSFP'):
                    sys.stdout.write("ITERATIONS = %d\n" % ((1000*context.scene.difsim.Delta,1000*context.scene.difsim.repeat_time)[context.scene.difsim.Delta < context.scene.difsim.repeat_time]))
                else:
                    sys.stdout.write("ITERATIONS = %d\n" % (1000*context.scene.difsim.Delta + 2000*context.scene.difsim.delta))
                    
            elif (context.scene.difsim.periodic and line.startswith("INSTANTIATE")):
                sys.stdout.write("\n")
                sys.stdout.write("DEFINE_SURFACE_CLASSES {\n")
                sys.stdout.write("  periodic_boundary {\n")
                sys.stdout.write("    PERIODIC = ALL_MOLECULES\n")
                sys.stdout.write("  }\n")
                sys.stdout.write("}\n")
                sys.stdout.write("\n")
                sys.stdout.write("periodic_box BOX {\n")
                sys.stdout.write("  CORNERS = [%0.6f, %0.6f, %0.6f], [ %0.6f, %0.6f, %0.6f]\n" % (-context.scene.difsim.xDim/2., -context.scene.difsim.yDim/2.,-context.scene.difsim.zDim/2., context.scene.difsim.xDim/2., context.scene.difsim.yDim/2., context.scene.difsim.zDim/2.))
                sys.stdout.write("  DEFINE_SURFACE_REGIONS {\n")
                sys.stdout.write("    interface {\n")
                sys.stdout.write("      ELEMENT_LIST = [ALL_ELEMENTS]\n")
                sys.stdout.write("      SURFACE_CLASS = periodic_boundary\n")
                sys.stdout.write("    }\n")
                sys.stdout.write("  }\n")
                sys.stdout.write("}\n")
                sys.stdout.write("\n")
                sys.stdout.write(line)
                counter = 1
            elif (context.scene.difsim.periodic and counter == 0):
                sys.stdout.write("  periodic_box OBJECT periodic_box {}\n")
                counter += -1
                sys.stdout.write(line)
            else:
                counter += -1
                sys.stdout.write(line)

        if context.scene.difsim.use_relaxation and context.scene.difsim.molecules.molecule_list:
            # If using relaxation and have a molecule list use that to write .molecules.mdl
            with open(os.path.join(mcell_project_dir, base_name) + '.molecules.mdl', 'w') as molout:
                molout.write("DEFINE_MOLECULES\n")
                molout.write("{\n")
                for mol in context.scene.difsim.molecules.molecule_list:
                    molout.write("  %s\n" % mol.name)
                    molout.write("  {\n")
                    molout.write("      DIFFUSION_CONSTANT_3D = %0.8g\n" % mol.diff_const)
                    print(mol.relax_t1)
                    molout.write("      T1 = %d\n" % (1000.*mol.relax_t1))
                    molout.write("      T2 = %d\n" % (1000.*mol.relax_t2))
                    molout.write("  }\n")
                molout.write("}")



class DIFSIM_OT_rundifsim(bpy.types.Operator):
    bl_idname = "difsimgui.rundifsim"
    bl_label = "Run Difsim on Current Model"
    bl_description = "Run Difsim on Current Model"
    bl_options = {'REGISTER'}


    def execute(self, context):
        mcell_project_dir = mcell_files_path()
        config_file = os.path.join(mcell_project_dir, 'difsim_params.txt') 
        context.scene.difsim.export_difsim(context, config_file)
        difsim_binary_path = context.scene.difsim.difsim_binary
        sp = subprocess.Popen([difsim_binary_path,config_file])#, stdout=None, stderr=None)
        return {'FINISHED'}

class DIFSIM_OT_plot_signal(bpy.types.Operator):
    bl_idname = "difsimgui.plot_signal"
    bl_label = "Display Plot of Signal"
    bl_description = "Display Plot of Signal"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        
    
        # Look up the plotting module by its name
        for plot_module in cellblender.cellblender_info[
                'cellblender_plotting_modules']:
            mod_name = plot_module.get_name()
            if mod_name == "MatPlotLib Plotter":
                break

        #Reorganize signal file to have better format for plotter
        #data_path = os.path.join(os.path.dirname(bpy.data.filepath),"signal.dat")
        data_path = os.path.join(mcell_files_path(),"signal")
        dat = numpy.genfromtxt(os.path.join(data_path, context.scene.difsim.signal_file))
        plot_spec_string = ""
        for i in range(dat.shape[1]):
            column = dat[:,i]
            x = numpy.arange(len(column))
            y = column

            numpy.savetxt(os.path.join(data_path,"plotsignal_%d" % i), numpy.vstack([x,y]).T)
            plot_spec_string += " f=plotsignal_%d" % i
            
        #plot_spec_string = "f=signal.dat"
        plot_module.plot(data_path, plot_spec_string)
        return{'FINISHED'}

#Eventually show gradient directions over model
class DIFSIM_OT_show_tensor_ball(bpy.types.Operator):
    bl_idname = "difsimgui.show_tensor_ball"
    bl_label = "Display Gradient Direction Vectors"
    bl_description = "Display Gradient Direction Vectors"
    bl_options = {'REGISTER', 'UNDO'}

    def _parsetensor(self, gradient_directions):
        tensorfile = os.path.join(os.path.dirname(__file__), 'tensor.dat')
        dat = numpy.array(open(tensorfile).readlines())
        testl = '%s\n' % (gradient_directions)
        w = int(numpy.where(dat==testl)[0]+1)
        v = numpy.array([curd.split(' ') for curd in dat[w:w+gradient_directions]], dtype=float)
        return v

    def execute(self, context):
        gradient_directions = context.scene.difsim.gradient_directions
        all_vec = self._parsetensor(gradient_directions)
        #glColor3f(0.5, 0.5, 1)
        #glLineWidth(2.0)

        vector_data = {}
        for i in range(len(all_vec)):
            vec = all_vec[i]
            print( "Creating gradient vector: " + str(vec))
            pvec = mathutils.Vector((vec[0], vec[1], vec[2]))
            vector_data[i] = pvec
            #glBegin(GL_LINE_STRIP)
            #glVertex3f(*mathutils.Vector((0.0,0.0,0.0)))
            #glVertex3f(*(100*pvec))
            #glEnd()
            #glPointSize(1.0)
        mathviz_draw.callback_enable()

        return{'FINISHED'}

class DIFSIM_OT_calc_diffstats(bpy.types.Operator):
    bl_idname = "difsimgui.calc_diffstats"
    bl_label = "Calculate Diffusion Statistics"
    bl_description = "Calculate Diffusion Statistics"
    bl_options = {'REGISTER', 'UNDO'}

    def parsetensor(self, gradient_directions):
        tensorfile = os.path.join(os.path.dirname(__file__), 'tensor.dat')
        dat = numpy.array(open(tensorfile).readlines())
        testl = '%s\n' % (gradient_directions)
        w = int(numpy.where(dat==testl)[0]+1)
        v = numpy.array([curd.split(' ') for curd in dat[w:w+gradient_directions]], dtype=float)
        return v
        
    def execute(self, context):
        allsig = []
        allbvec = []
        #for i in range(context.scene.difsim.signals.next_id):
        if not context.scene.difsim.signals.signal_list:
            if not os.path.exists(context.scene.difsim.signal_file):
                print("No simulations have been run and no signal files in list.")
                return{'FINISHED'}
            else:
                context.scene.difsim.signals.add_signal(context, context.scene.difsim.signal_file)

        #If no files selected, add most recent simulation, else break and analyze
        add_to_list = 1
        for sig in context.scene.difsim.signals.signal_list:
            if sig.analyze:
                add_to_list = 0
                break
        if add_to_list:
            context.scene.difsim.signals.add_signal(context, context.scene.difsim.signal_file)

        for sig in context.scene.difsim.signals.signal_list:
            if sig.analyze:
                for line in numpy.genfromtxt(sig.filename):
                    allsig.append(line) 
                grad_dir = self.parsetensor(sig.nDirs)
                for gd in grad_dir:
                    cbv = sig.bvalue*numpy.array([gd[0]**2, gd[1]**2, gd[2]**2, gd[0]*gd[1],gd[0]*gd[2], gd[1]*gd[2]])
                    allbvec.append(cbv)
        if not allsig:
            print("No signal files are selected to be analyzed")
            return{'FINISHED'}
        allbvec = numpy.array(allbvec, dtype=float)
        bmatrix_fn = os.path.join(mcell_files_path(),"signal/bmatrix.txt")
        numpy.savetxt(bmatrix_fn, allbvec) 
        signal_fn = os.path.join(mcell_files_path(), "signal/signal.dat")
        numpy.savetxt(signal_fn, allsig)
        signal_nifti_fn = os.path.join(mcell_files_path(), "signal/signal.nii")
        dtoutput_fn = os.path.join(mcell_files_path(),"signal/output.1D")
        eigsoutput_fn =os.path.join(mcell_files_path(),"signal/eig.1D")

        python_path = context.scene.mcell.cellblender_preferences.python_binary
        save_nifti_script = os.path.join(os.path.dirname(__file__), 'save_nifti.py')
        cmd = [python_path, save_nifti_script, signal_fn, signal_nifti_fn]
        sp = subprocess.Popen(cmd)
        sp.communicate()

        afnicmd = "3dDWItoDT"
        prefix = "-prefix"
        bmarg = "-bmatrix_NZ"
        
        cmd = [afnicmd,prefix,dtoutput_fn,bmarg,bmatrix_fn,signal_nifti_fn] 
        print( " ".join(cmd))
        sp = subprocess.Popen(cmd)
        sp.communicate()

        eigcmd = "3dDTeig"
        cmd = [eigcmd,prefix,eigsoutput_fn,dtoutput_fn]
        sp = subprocess.Popen(cmd)
        sp.communicate()

        return{'FINISHED'}

class DIFSIM_PT_DifSimGui(bpy.types.Panel):
    bl_label = "Difsim Model"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "Difsim"
    bl_options = {'DEFAULT_CLOSED'}

    def draw(self, context):
        if context.scene != None:
            context.scene.difsim.draw_panel(context, panel=self)

class DIFSIM_OT_signal_add(bpy.types.Operator):
    bl_idname = "difsim.signal_add"
    bl_label = "Add Signal"
    bl_description = "Add a new signal file to the analysis pipeline"
    bl_options = {'REGISTER', 'UNDO'}

    filepath = bpy.props.StringProperty(subtype='FILE_PATH', default="")
    
    def execute(self, context):
        context.scene.difsim.signals.add_signal(context, self.filepath)
        return {'FINISHED'}

    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

class DIFSIM_OT_signal_remove(bpy.types.Operator):
    bl_idname = "difsim.signal_remove"
    bl_label = "Remove Signal"
    bl_description = "Remove selected signal from analysis pipeline"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        context.scene.difsim.signals.remove_active_signal(context)
        self.report({'INFO'}, "Deleted Signal")
        return {'FINISHED'}

class DifsimSignalProperty(bpy.types.PropertyGroup):
    name = StringProperty(name="Signal Name", default="Signal",description="The user specified signal name")
    filename = StringProperty(name="Signal File Path", default="", description="Location where signal was loaded from")
    id = IntProperty(name="Signal ID", default=0)
    bvalue = IntProperty(name="b-factor", default=0)
    nDirs = IntProperty(name="number of directions", default = 1)
    analyze = BoolProperty(default=True)
    status = StringProperty(name="Status")

    def init_properties ( self, fn ):
        self.name = str.split(fn,"/")[-1]
        self.filename = fn 
        #self.bvalue = int(re.search(r'\d+', fn).group())
        self.bvalue = int(str.split(str.split(fn,".")[0],"_")[-1])
        self.nDirs = sum(1 for line in open(fn)) - 2   #First line is header, second line is times

    def draw_props ( self, layout):
        layout.prop ( self, "name" )
        layout.prop ( self, "bvalue")
        layout.prop ( self, "nDirs")
        layout.prop (self, "filename") 
        layout.prop (self, "analyze")
  
class DifsimSignalListProperty(bpy.types.PropertyGroup):
    signal_list = CollectionProperty(type=DifsimSignalProperty, name="Signal List")
    active_sig_index = IntProperty(name="Active Signal Index", default=0)
    #show_display = bpy.props.BoolProperty(default=False)  # If Some Properties are not shown, they may not exist!!!
    #show_advanced = bpy.props.BoolProperty(default=False)  # If Some Properties are not shown, they may not exist!!!
    next_id = IntProperty(name="Counter for Unique Signal IDs", default=1)  # Start ID's at 1 to confirm initialization

    def init_properties ( self ):
        if self.signal_list:
            for sig in self.signal_list:
                sig.init_properties("")

    def remove_properties ( self, context ):
        print ( "Removing all Signal List Properties..." )
        for item in self.signal_list:
            item.remove_properties(context)
        self.signal_list.clear()
        self.active_sig_index = 0
        self.next_id = 1
        print ( "Done removing all Signal List Properties." )
        
    
    def add_signal( self, context, fn ):
        """ Add a new molecule to the list of signals and set as the active signal """
        new_sig = self.signal_list.add()
        new_sig.id = self.allocate_available_id()
        new_sig.init_properties(fn)
        self.active_sig_index = len(self.signal_list)-1

    def remove_active_signal ( self, context ):
        """ Remove the active molecule from the list of signals"""
        self.signal_list.remove ( self.active_sig_index )
        self.active_sig_index -= 1
        if self.active_sig_index < 0:
            self.active_sig_index = 0
        if self.signal_list:
            self.check(context)
    
    def check ( self, context ):
        """Checks for duplicate or illegal molecule name"""
        # Note: Some of the list-oriented functionality is appropriate here (since this
        #        is a list), but some of the molecule-specific checks (like name legality)
        #        could be handled by the the signals themselves. They're left here for now.

        sig = self.signal_list[self.active_sig_index]

        status = ""

        # Check for duplicate molecule name
        sig_keys = self.signal_list.keys()
        if sig_keys.count(sig.name) > 1:
            status = "Duplicate signal: %s" % (sig.name)

        # Check for illegal names (Starts with a letter. No special characters.)
        sig_filter = r"(^[A-Za-z]+[0-9A-Za-z_.]*$)"
        m = re.match(sig_filter, sig.name)
        if m is None:
            status = "Signal name error: %s" % (sig.name)

        sig.status = status

        return


    def allocate_available_id ( self ):
        """ Return a unique molecule ID for a new molecule """
        if len(self.signal_list) <= 0:
            # Reset the ID to 1 when there are no more signals 
            self.next_id = 1
        self.next_id += 1
        return ( self.next_id - 1 )


    def draw_layout ( self, context, layout ):
        """ Draw the signal "panel" within the layout """
        difsim = context.scene.difsim
        box = layout.box()
        row = box.row()
        row.label("Signal List")
        row = box.row()
        col = row.column()
        col.template_list("DIFSIM_UL_check_signal", "define_signals",
                          self, "signal_list",
                          self, "active_sig_index",
                          rows=2)
        col = row.column(align=True)
        col.operator("difsim.signal_add", icon='ZOOMIN', text="")
        col.operator("difsim.signal_remove", icon='ZOOMOUT', text="")
        if self.signal_list:
            sig = self.signal_list[self.active_sig_index]
            # The self is needed to pass the "advanced" flag to the molecule
            sig.draw_props ( box )


    def draw_panel ( self, context, panel ):
        """ Create a layout from the panel and draw into it """
        layout = panel.layout
        self.draw_layout ( context, layout )

class DifsimMoleculeProperty(bpy.types.PropertyGroup):
    name = StringProperty(name="Molecule Name", default="Molecule",description="The user specified signal name")
    id = IntProperty(name="Molecule ID", default=0)
    diff_const = FloatProperty(name = "DIFFUSION_CONSTANT_3D", default = 0.001)
    relax_t1 = FloatProperty(name = "T1 [ms]", default = 4000)
    relax_t2 = FloatProperty(name = "T2 [ms]", default = 2000)

    def init_properties ( self, mcell_mname, mcell_diff_const):
        self.name = mcell_mname
        self.diff_const = mcell_diff_const

    def draw_props ( self, layout):
        layout.prop ( self, "name" )
        layout.prop ( self, "diff_const")
        layout.prop ( self, "relax_t1")
        layout.prop ( self, "relax_t2")

    def remove_properties ( self, context ):
        print ( "Removing all Molecule Properties ... not implemented yet!" )

class DifsimMoleculeListProperty(bpy.types.PropertyGroup):
    molecule_list = CollectionProperty(type=DifsimMoleculeProperty, name="Molecule List")
    active_mol_index = IntProperty(name="Active Molecule Index", default=0)
    next_id = IntProperty(name="Counter for Unique Molecule IDs", default=1)  # Start ID's at 1 to confirm initialization

    def init_properties ( self ):
        if self.molecule_list:
            for mol in self.molecule_list:
                mol.init_properties()

    def add_molecule( self, context, mcell_mname, mcell_diff_const):
        """ Add a new molecule to the list of molecules and set as the active molecule """
        new_mol = self.molecule_list.add()
        #new_mol.id = self.allocate_available_id()
        new_mol.init_properties(mcell_mname, mcell_diff_const)
        #self.active_mol_index = len(self.molecule_list)-1

    def remove_properties ( self, context ):
        print ( "Removing all Signal List Properties..." )
        for item in self.molecule_list:
            item.remove_properties(context)
        self.molecule_list.clear()
        self.active_sig_index = 0
        self.next_id = 1
        print ( "Done removing all Molecule List Properties." )
        
        
    def draw_layout ( self, context, layout ):
        """ Draw the signal "panel" within the layout """
        difsim = context.scene.difsim
        box = layout.box()
        row = box.row()
        row.label("Molecule List")
        row = box.row()
        col = row.column()
        col.template_list("DIFSIM_UL_check_molecule", "define_molecules",
                          self, "molecule_list",
                          self, "active_mol_index",
                          rows=2)
        col = row.column(align=True)
        #for mcell_mol in context.scene.mcell.molecules.molecule_list:
        #    self.add_molecule(self, mcell_mol.name)
        if self.molecule_list:
            mol = self.molecule_list[self.active_mol_index]
            # The self is needed to pass the "advanced" flag to the molecule
            mol.draw_props ( box )

    def draw_panel ( self, context, panel ):
        """ Create a layout from the panel and draw into it """
        layout = panel.layout
        self.draw_layout ( context, layout )

class DIFSIM_UL_check_signal(bpy.types.UIList):
    def draw_item(self, context, layout, data, item, icon, active_data, active_propname, index):
        if item.analyze:
            layout.label(item.name, icon='FILE_TICK')
        else:
            layout.label(item.name)#, icon='')

class DIFSIM_UL_check_molecule(bpy.types.UIList):
    def draw_item(self, context, layout, data, item, icon, active_data, active_propname, index):
            layout.label(item.name)#, icon='')

class DIFSIM_PT_define_signals(bpy.types.Panel):
    bl_label = "Difsim - Signals"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "scene"
    bl_options = {'DEFAULT_CLOSED'}

    def draw ( self, context ):
        # Call the draw function for the instance being drawn in this panel
        context.scene.difsim.signals.draw_panel ( context, self )

class DIFSIM_PT_define_molecules(bpy.types.Panel):
    bl_label = "Difsim - Molecules"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "scene"
    bl_options = {'DEFAULT_CLOSED'}

    def draw ( self, context ):
        # Call the draw function for the instance being drawn in this panel
        context.scene.difsim.molecules.draw_panel ( context, self )

class DifSimObjectProperty(bpy.types.PropertyGroup):
    periodic = BoolProperty(name="periodic", default=False)
    xDim = FloatProperty(name="x dim", default = 50.)
    yDim = FloatProperty(name="y dim", default = 50.)
    zDim = FloatProperty(name="z dim", default = 50.)
    
    pulsetype_enum = [
                          ('gradient echo','gradient echo',''),
                          ('spin echo','spin echo',''),
                          ('dPFG','dPFG',''),
                          ('DWSSFP','DWSSFP',''),
                          ('multiple spin echo', 'multiple spin echo', '')
                        ]

    pulsetype = EnumProperty(items=pulsetype_enum, name="Pulse type", description="Gradient pulse protocol")
    ramp = FloatProperty(name="ramp [ms]", default=0.)
    Delta = FloatProperty(name="Delta [ms]", default=30.)
    delta = FloatProperty(name="delta [ms]", default=2.)
    gradient_strength = FloatProperty(name="gradient strength", default=4.)
    gradient_directions = IntProperty(name="gradient directions", default=32)
    workfile = StringProperty(name="work file", default='difsim_param.txt')

    difsim_binary = StringProperty(name="difsim binary", default='')
    mcell_file = StringProperty(name="mcell file", default='difsim.main.mdl')
    signal_file = StringProperty(name="signal file", default='signal.dat')
    tensor_file = StringProperty(name="tensor file", default='tensor.dat')

    use_relaxation = BoolProperty(name = "Relaxation", default=False)
    repeat_time = IntProperty(name="TR", default = 0)
    repeat_number = IntProperty(name="repeat number", default = 1)
    flip_angle = FloatProperty(name="flip angle", default = 90)

    #Only used for ncCurrent sims
    calculate_ncCurrent = BoolProperty(name="Calculate ncCurrent", default=False)
    ncCurrent_strength = FloatProperty(name="ncCurrent Strength", default=1.0)
    ncCurrent_output_file = StringProperty(name="ncCurrent Output", default="./ncoutput.txt")
    #only show for DPFG
    is_dpfg = BoolProperty(name="Using dPFG", default=False)
    mixing_time = IntProperty(name="mixing time [ms]", default=50)
    DPFG_angle = IntProperty(name="DPFG angle", default=0)
    DPFG_dde_72 = BoolProperty(name="DPFG 5-design (72 directions)", default=False)

    echo_time = IntProperty(name="echo time (TE) [ms]", default = 0)
    echo_number = IntProperty(name="echo number (NE)", default = 1)
    echo_spacing = IntProperty(name="echo spacing (ESP) [ms]", default = 0)

    show_DT = BoolProperty(name="Show DT", default=False)
    use_ge_tensor = BoolProperty(name="Use GE Tensor Gradient Directions", default=True)

    signals = PointerProperty(type=DifsimSignalListProperty, name="Available Signals")
    molecules = PointerProperty(type=DifsimMoleculeListProperty, name = "Available Molecules")


    def get_numpulsetype(self):
        if self.pulsetype == 'gradient echo': return 0 
        elif self.pulsetype == 'spin echo': return 1
        elif self.pulsetype == 'dPFG': return 2
        elif self.pulsetype == 'DWSSFP': return 3
        elif self.pulsetype == 'multiple spin echo': return 4
        else: return -1

    def export_difsim(self, context, workfile):
        print( "Exporting Difsim Config File" )
        out = open(workfile, 'w')
        result = []
        if (self.periodic):
            result.append('periodic = 1')
        result.append('')
        result.append('x dim = %0.6f' % (self.xDim))
        result.append('y dim = %0.6f' % (self.yDim))
        result.append('z dim = %0.6f' % (self.zDim))
        result.append('')
        if (self.use_relaxation):
            result.append('relaxation = 1')
            result.append('')
        result.append('pulse = %d' % (self.get_numpulsetype()))
        result.append('delta = %d' % (1000.*self.delta))
        result.append('Delta = %d' % (1000.*self.Delta))
        result.append('gradient strength = %0.1f' % (self.gradient_strength))
        result.append('gradient directions = %d' % (self.gradient_directions))
        result.append('ramp = %d' % (self.ramp))
        ts = bpy.context.scene.mcell.initialization.time_step.get_value()
        result.append('step size = %0.1f' % (1e6*bpy.context.scene.mcell.initialization.time_step.get_value()))
        if (self.calculate_ncCurrent):
            result.append('')
            result.append('ncCurrent strength = %0.1f' % (self.ncCurrent_strength))
            result.append('ncCurrent output file = %s' % (self.ncCurrent_output_file))
        if (self.pulsetype == "DWSSFP"):
            result.append('')
            result.append('TR = %d' % (1000.*self.repeat_time))
            result.append('repeat number = %d' % (self.repeat_number))
            result.append('flip angle = %0.3f' % (self.flip_angle))
        if (self.pulsetype == "dPFG"):
            result.append('')
            result.append('mixing time = %d' % (1000.*self.mixing_time))
            if (self.DPFG_dde_72):
                result.append('tensor file DPFG = 1')
                result.append('tensor file = %s' % (os.path.join(os.path.dirname(__file__),'tensor_dde_72.dat')))
            else:
                result.append('DPFG angle = %d' % (self.DPFG_angle))
        if (self.pulsetype == "multiple spin echo"):
            result.append('')
            result.append('TE = %d' % (1000.*self.echo_time))
            result.append('NE = %d' % (self.echo_number))
            result.append('ESP = %d' % (1000.*self.echo_spacing))
        result.append('')

        mcell_project_dir = mcell_files_path()
        base_name = context.scene.mcell.project_settings.base_name
        self.mcell_file = os.path.join(mcell_project_dir, base_name) + '.main.mdl'
        result.append('mcell file = %s' % (self.mcell_file))
        if (self.use_ge_tensor and self.pulsetype != "dPFG"):
            self.tensor_file = os.path.join(os.path.dirname(__file__), 'tensor.dat')
            result.append('tensor file = %s' % (self.tensor_file))
        signal_path = os.path.join(mcell_files_path(),"signal")
        if not os.path.exists(signal_path):
            os.makedirs(signal_path)
        self.signal_file = os.path.join(signal_path, "signal.dat")
        
        bfactor = 26752*26752 * context.scene.difsim.gradient_strength/10.0 * context.scene.difsim.gradient_strength/10.0 *(context.scene.difsim.delta/1000.0) * (context.scene.difsim.delta/1000.0) * (context.scene.difsim.Delta + 2.0*context.scene.difsim.delta/3.0)/1000.0
       
        sig_with_bfactor = "signal_%d.dat" % (bfactor)
        self.signal_file = os.path.join(signal_path, sig_with_bfactor)
        context.scene.difsim.signal_file = self.signal_file

        result.append('signal file = %s' % (self.signal_file))
        out.write('\n'.join(result))
        print("Finished Exporting")

    def draw_panel(self, context, panel):
        layout = panel.layout
        self.draw_layout(context, panel)

    def draw_DT(self, context, panel):
        layout = panel.layout
        box = layout.box()

        dtoutput_fn = os.path.join(mcell_files_path(),"signal/output.1D")
        eigsoutput_fn = os.path.join(mcell_files_path(),"signal/eig.1D")
        if (os.path.exists(dtoutput_fn)):
            dtcor = numpy.genfromtxt(dtoutput_fn )
        else:
            dtcor = numpy.zeros((6))
        if (os.path.exists(eigsoutput_fn )):
            eigs = numpy.genfromtxt(eigsoutput_fn)
        else:
            eigs = numpy.zeros((3))
        row = box.row()
        row.label("D = ") 
        inbox = row.box()
        
        split = inbox.split()
        col = split.column()
        col.label("%f" % (dtcor[0]))
        col = split.column()
        col.label("%f" % (dtcor[1]))
        col = split.column()
        col.label("%f" % (dtcor[3]))

        row = inbox.row()
        split = row.split()
        col = split.column()
        col.label("%f" % (dtcor[1]))
        col = split.column()
        col.label("%f" % (dtcor[2]))
        col = split.column()
        col.label("%f" % (dtcor[4]))

        row = inbox.row()
        split = row.split()
        col = split.column()
        col.label("%f" % (dtcor[3]))
        col = split.column()
        col.label("%f" % (dtcor[4]))
        col = split.column()
        col.label("%f" % (dtcor[5]))

        row = box.row()
        row.label("lambda = ")

        inbox = row.box()
        row = inbox.row()
        split = row.split()
        col = split.column()
        col.label("%f" % (eigs[0]))
        col = split.column()
        col.label("%f" % (eigs[1]))
        col = split.column()
        col.label("%f" % (eigs[2]))

        row = box.row()
        split = row.split()
        col = split.column()
        col.label("FA = %f" % (eigs[-2]))
        col = split.column()
        col.label("MD %f" % (eigs[-1]))

    def draw_layout(self, context, panel):
        mcell = context.scene.mcell
        if not mcell.initialized:
            mcell.draw_uninitialized ( layout )
        else:
            layout = panel.layout
            box = layout.box()
            box.label(text="Difsim Parameters", icon='MESH_CUBE')
            col = box.column()
            
            box.prop(self,"periodic")

            box.prop(self,"xDim")
            box.prop(self,"yDim")
            box.prop(self,"zDim")
            box.operator("difsimgui.update_bounds", text="Update Voxel Dimensions from MCell Model")

            box = layout.box()
            box.label(text="Gradient Pulse Sequence", icon='RNA')
            box.prop(self,"pulsetype")
            row = box.row()
            bval = 26752*26752 * self.gradient_strength/10.0 * self.gradient_strength/10.0 *(self.delta/1000.0) * (self.delta/1000.0) * (self.Delta + 2.0*self.delta/3.0)/1000.0
            row.label(text="b-value: %0.2f mm^-1" % bval)
            q = self.gradient_strength*self.delta*1e-3*26752/1e4
            row = box.row()
            row.label(text="q: %0.2f um^-1" % q) 
            box.prop(self, "gradient_strength", )
            box.prop(self, "ramp")
            box.prop(self, "delta")
            box.prop(self, "Delta")
            box.prop(self, "gradient_directions")
            box.prop(self, "use_relaxation")
            if (self.use_relaxation):
                box.operator("difsimgui.update_molecule_list", text="Update Molecule List") 
                self.molecules.draw_panel(context, panel)
                #TODO Figure out how to put molecules here so T1, T2 can be specified
                    
                
            if (self.pulsetype != "dPFG"):
                box.prop(self, "use_ge_tensor")
            if (self.pulsetype == "DWSSFP"):
                box.prop(self, "repeat_time")
                box.prop(self, "repeat_number")
                box.prop(self, "flip_angle")
            if (self.pulsetype == "dPFG"):
                box.prop(self, "mixing_time")
                box.prop(self, "DPFG_dde_72")
                if (self.DPFG_dde_72 == False):
                    box.prop(self, "DPFG_angle")
            if (self.pulsetype == "multiple spin echo"):
                box.prop(self, "echo_time")
                box.prop(self, "echo_number")
                box.prop(self, "echo_spacing")
            box = layout.box()
            row = box.row()
            row.prop(self,"calculate_ncCurrent")#, icon='FORCE_MAGNETIC')
            if (self.calculate_ncCurrent):
                row = box.row()
                row.prop(self,"ncCurrent_strength")
                row = box.row()
                row.prop(self,"ncCurrent_output_file")
            row = box.row()
            box = layout.box()
            row = box.row()
            if not os.path.dirname(bpy.data.filepath):
                row.label(text="Open or save a .blend file to set the project directory",icon='ERROR')
            else:
                box.operator("difsimgui.set_difsim_binary", text="Set Difsim Binary")
                box.operator("difsimgui.export_difsim", text="Export Difsim Config File")
                box.operator("difsimgui.rundifsim", text="Run Difsim Simulation of Current Model")
                box.operator("difsimgui.plot_signal", text="Plot Signal File")
                #signalbox.prop(self, "signals")
                self.signals.draw_panel(context, panel)
                #signalbox.operator("difsim.signal_add", text="Add Signal")
                #signalbox.operator("difsim.signal_remove", text="Remove Signal")
                python_path = context.scene.mcell.cellblender_preferences.python_binary
                if not python_path:
                    box.label("Set python path in CellBlender to compute statistics", icon='UNPINNED')
                else: 
                    box.operator("difsimgui.calc_diffstats", text="Calculate Diffusion Statistics")
                box.prop(self,"show_DT", text="Show Diffusion Statistics")
                if (self.show_DT):
                    self.draw_DT(context,panel)




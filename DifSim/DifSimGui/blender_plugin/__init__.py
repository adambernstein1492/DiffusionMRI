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

import os
import sys

bl_info = {
    "name": "Difsim Model Builder",
    "author": "Ben Regner",
    "version": (1, 0, 0),
    "blender": (2, 66, 1),
    "api": 55057,
    "location": "View3D > Tools > Difsim",
    "description": "Model Builder for Difsim Simulator",
    "warning": "",
    "wiki_url": "http://www.mcell.org",
    "tracker_url": "http://code.google.com/p/cellblender/issues/list",
    "category": "Cell Modeling"
}


# To support reload properly, try to access a package var.
# If it's there, reload everything
if "bpy" in locals():
    print("Reloading Difsim")
    import imp
    imp.reload(difsimgui)
else:
    print("Importing Difsim")
    from . import \
        difsimgui

import bpy

difsim_added_handlers = []

def add_handler ( handler_list, handler_function ):
    """ Only add a handler if it's not already in the list """
    if not (handler_function in handler_list):
        handler_list.append ( handler_function )

        difsim_added_handlers

def add_handler ( handler_list, handler_function ):
    """ Only add a handler if it's not already in the list """
    if not (handler_function in handler_list):
        handler_list.append ( handler_function )

        difsim_added_handlers

def register():
    bpy.utils.register_module(__name__)
    bpy.types.Scene.difsim = bpy.props.PointerProperty(
        type=difsimgui.DifSimObjectProperty)
    add_handler ( bpy.app.handlers.load_post, difsimgui.load_preferences )

    print("DifSim registered")


def unregister():
    bpy.utils.unregister_module(__name__)
    remove_handler ( bpy.app.handlers.load_post, difsimgui.load_preferences )

    print("DifSim unregistered")



#!/usr/bin/env python
import xml.etree.ElementTree as ET
import bp
import command
import vis
import grid2
import contour2
import btypes


class Framework(command.CommandReceiver):
    """ Presents data collection for single work flow.
    """
    _visCls = vis.QtFrameworkVis

    def __init__(self):
        #grids list
        self.grids2 = bp.NamedList()
        #user contours list
        self.contours2 = bp.NamedList()
        #boundary types
        self.boundary_types = btypes.BndTypesList()
        #visualizer set
        self._visualiser = self._visCls()
        self._visualiser.set_framework(self)

    def view_update(self):
        'updates visualisation'
        self._visualiser.update()

    def post_proc(self):
        'callback after each command invocation '
        self.view_update()

    #grids manager
    def add_grid(self, name, grid):
        "-> (index, name, grid). add a grid to a grid list"
        self.grids2[name] = grid
        return self.grids2.get(val=grid)

    def remove_grid(self, name):
        'removes a grid with the certain name from a grid list'
        self.grids2.pop(name)

    def get_grid(self, ind=None, name=None, grid=None):
        "-> (index, name, grid). Get a grid by name or reference"
        return self.grids2.get(ind, name, grid)

    def get_grid_names(self):
        '-> [list of grid names]'
        return self.grids2.keys()

    def get_checked_grid_names(self):
        '-> [list of user checked grid names]'
        return self._visualiser.ask_checked_grid_names()

    def get_unchecked_grid_names(self):
        '-> [list of user unchecked grid names]'
        return [x for x in self.get_grid_names()
                if x not in self.get_checked_grid_names()]

    #boundary types manager
    def get_bnd_types(self):
        '-> btypes.BndTypesList'
        return self.boundary_types

    #contours manager
    def add_user_contour(self, name, cont):
        "-> (index, name, contour). add a user contour"
        self.contours2[name] = cont
        return self.contours2.get(val=cont)

    def get_user_contour(self, ind=None, name=None, cont=None):
        "-> (index, name, grid). Get a contour by name or reference"
        return self.contours2.get(ind, name, cont)

    def get_any_contour(self, name):
        '-> Grid or user contour by its name'
        if name in self.contours2:
            return self.contours2[name]
        else:
            return self.grids2[name].cont

    def remove_user_contour(self, name):
        "removes contour from framework"
        del self.contours2[name]

    def get_contour_names(self):
        '-> [list of user contour names]'
        return self.contours2.keys()

    def get_checked_contour_names(self):
        '-> [list of user checked contour names]'
        return self._visualiser.ask_checked_contours_names()

    def get_checked_gbnd_names(self):
        '-> [list of user checked grid boundaries names]'
        return self._visualiser.ask_checked_gbnd_names()

    def get_checked_any_contour_names(self):
        '-> [list of user checked grid and user contours]'
        return self.get_checked_contour_names() + \
                self.get_checked_gbnd_names()

    def get_all_names(self):
        '-> [list of all grid and user contours]'
        return self.get_contour_names() + \
                self.get_grid_names()
    #user requests
    def ask_for_new_contours_bnd(self, cont_names):
        """ ->[{boundry-index: [list of contour edges]}, {}, ...]
            request for new boundary types for contours
        """
        return self._visualiser.ask_new_contours_bnd(cont_names,
                self.boundary_types)

    #overriden from CommandReceiver
    def to_zero_state(self):
        'deletes all grids and contours'
        self.grids2.clear()
        self.contours2.clear()
        self.boundary_types.clear()

    def save_state(self, xmlnode):
        ' writes data to @xmlNode of xml.etree.ElementTree.Element class'
        #Grids
        gnode = ET.SubElement(xmlnode, "GRIDS")
        for k, v in self.grids2.items():
            nd = ET.SubElement(gnode, "GRID2")
            nd.attrib["name"] = k
            v.xml_save(nd)
        #User contours
        cnode = ET.SubElement(xmlnode, "CONTOURS")
        for k, v in self.contours2.items():
            nd = ET.SubElement(cnode, "CONTOUR2")
            nd.attrib["name"] = k
            v.xml_save(nd)
        #boundary types
        bnode = ET.SubElement(xmlnode, "BTYPES")
        self.boundary_types.xml_save(bnode)

    def load_state(self, xmlnode):
        ' adds data from ElementTree Element object '
        #boundary types
        bnode = xmlnode.find("BTYPES")
        if bnode is not None:
            self.boundary_types.xml_load(bnode)
        #Grids
        gnodes = xmlnode.findall("GRIDS/GRID2")
        for nd in gnodes:
            name = nd.attrib["name"]
            grd = grid2.Grid2.create_from_xml(nd)
            self.grids2[name] = grd
        #User contours
        cnodes = xmlnode.findall("CONTOURS/CONTOUR2")
        for nd in cnodes:
            name = nd.attrib['name']
            cont = contour2.AbstractContour2.create_from_xml(nd)
            self.contours2[name] = cont
        self.view_update()

    def deep_copy(self):
        ret = Framework()
        #boundarys
        ret.boundary_types.add_data(self.boundary_types)
        #grids
        for name, g in self.grids2.items():
            gnew = g.deepcopy()
            ret.grids2[name] = gnew
        #contours
        for name, g in self.contours2.items():
            gnew = g.deepcopy()
            ret.contours2[name] = gnew
        return ret

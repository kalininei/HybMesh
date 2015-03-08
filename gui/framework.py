#!/usr/bin/env python
import xml.etree.ElementTree as ET
import bp
import command
import vis
import grid2


class Framework(command.CommandReceiver):
    """
        Presents data collection for single work flow.
    """
    _visCls = vis.QtFrameworkVis

    def __init__(self):
        self._visualiser = self._visCls()
        self._visualiser.set_framework(self)
        self.grids2 = bp.NamedList()

    def view_update(self):
        'updates visualisation'
        self._visualiser.update()

    def post_proc(self):
        'callback after each command invokation '
        self.view_update()

    def get_grid_names(self):
        return self.grids2.keys()

    def get_checked_grid_names(self):
        return self._visualiser.get_checked_grid_names()

    def get_unchecked_grid_names(self):
        return [x for x in self.get_grid_names()
                if x not in self.get_checked_grid_names()]

    #overriden from CommandReceiver
    def to_zero_state(self):
        ' deletes all data '
        self.grids2 = bp.NamedList()

    def save_state(self, xmlnode):
        ' writes data to @xmlNode of xml.etree.ElementTree.Element class'
        #Grids
        gnode = ET.SubElement(xmlnode, "GRIDS")
        for k, v in self.grids2.items():
            nd = ET.SubElement(gnode, "GRID2")
            nd.attrib["name"] = k
            v.xml_save(nd)

    def load_state(self, xmlnode):
        ' adds data from ElementTree Element object '
        #Grids
        gnodes = xmlnode.findall("GRIDS/GRID2")
        for nd in gnodes:
            name = nd.attrib["name"]
            grd = grid2.Factory().xml_create(nd)
            self.grids2[name] = grd
        self.view_update()

    def deep_copy(self):
        ret = Framework()
        for name, g in self.grids2.items():
            gnew = g.deepcopy()
            ret.grids2[name] = gnew
        return ret



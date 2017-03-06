import xml.etree.ElementTree as ET
from hybmeshpack.com import factory
from hybmeshpack.hmcore import hmxml
from hybmeshpack.gdata import framework
import native_import
from hybmeshpack.gdata.cont2 import Contour2
from hybmeshpack.gdata.srf3 import Surface3
from hybmeshpack.gdata.grid2 import Grid2
from hybmeshpack.gdata.grid3 import Grid3


def _load_command_flow(comflow, xmlnode):
    'adds commands to flow.CommandFlow comflow from xmlnode tagged FLOW'
    # Commands list
    cfnd = xmlnode.findall("COMMANDS/COM")
    # reset all commands because they can not be undone after this procedure
    if len(cfnd) > 0:
        comflow.purge_flow()
    for cline in cfnd:
        try:
            com_title = cline.attrib["name"]
        except:
            continue
        try:
            com_string = cline.find("LINE").text
        except:
            com_string = ""
        try:
            com_comment = cline.find("COMMENT").text
        except:
            com_comment = ""
        if com_title != 'Start':
            c = factory.create_from_string(com_title, com_string)
            comflow.append_command(c)
        else:
            c = comflow._commands[0]
        c.set_comment(com_comment)
        if 'current' in cline.attrib:
            comflow._startpos = len(comflow._commands) - 1
            comflow._curpos = comflow._startpos
        if 'end_load_here' in cline.attrib:
            break


def _load_btypes(fw, node):
    for n in node.findall("BTYPE"):
        ind = int(n.attrib["index"])
        nm = n.find("NAME").text
        fw.set_zone_type(ind, nm)


def flow_and_framework_fromfile(filename, flow, cb=None):
    doc, root, statenodes = 0, 0, []
    try:
        cb.pycall("Loading commands", "Reading xml", 0, 0)
        doc, root = hmxml.open_doc(filename)
        pstring = hmxml.purged_string(doc)
        pyside_root = ET.fromstring(pstring)
        pyside_flownode = pyside_root.findall('FLOW')
        if len(pyside_flownode) == 0:
            raise Exception('No proper data in xml file')
        else:
            pyside_flownode = pyside_flownode[0]

        # load commands
        cb.pycall("Loading commands", "Parsing", 0.5, 0)
        _load_command_flow(flow, pyside_flownode)
        cb.pycall("Loading commands", "Done", 1, 1)

        # state
        statenodes = hmxml.query(root, "FLOW/STATE")
        pyside_statenode = pyside_flownode.findall('STATE')
        if len(statenodes) > 0:
            statenode = statenodes[0]
            pyside_statenode = pyside_statenode[0]
        else:
            statenode = None
            pyside_statenode = None

        data = framework.Framework()
        if statenode is not None:
            # load boundary types
            _load_btypes(data, pyside_statenode)

            # load data
            c, g, s3, g3, cn, gn, s3n, g3n =\
                native_import.all_geom(doc, statenode, cb)
            for k, v in zip(cn, c):
                data.append_contour2(Contour2(v), k)
            for k, v in zip(gn, g):
                data.append_grid2(Grid2(v), k)
            for k, v in zip(s3n, s3):
                data.append_surface3(Surface3(v), k)
            for k, v in zip(g3n, g3):
                data.append_grid3(Grid3(v), k)
        flow.set_receiver(data)
    except:
        raise
    finally:
        hmxml.close_doc(doc, [root] + statenodes)

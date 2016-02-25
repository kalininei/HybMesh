'module imlements all write from xml routines'
import xml.etree.ElementTree as ET
from hybmeshpack import progdata


def _root_xml():
    'returns root xml node tagged HybMeshProject'
    ret = ET.Element('HybMeshProject')
    ret.attrib['ver'] = str(progdata.HybMeshVersion.current())
    return ret


def xmlindent(elem, level=0):
    """ http://effbot.org/zone/element-lib.htm#prettyprint.
        It basically walks your tree and adds spaces and newlines so the tree i
        printed in a nice way
    """
    tabsym = "    "
    i = "\n" + level * tabsym
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + tabsym
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            xmlindent(elem, level + 1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


def writexml(r, filename):
    xmlindent(r)
    tree = ET.ElementTree(r)
    tree.write(filename, xml_declaration=True, encoding='utf-8')


# ------------------- Commands
def write_command_flow(comflow, xmlnode, flowname=None):
    """ creates xmlnode/FLOW and writes flow commands there """
    rt = ET.SubElement(xmlnode, "FLOW")
    if flowname is not None:
        rt.attrib['name'] = flowname
    #Commands list
    cl = ET.SubElement(rt, "COMMANDS")
    for c in comflow._commands:
        nd = ET.SubElement(cl, "ENTRY")
        if c is comflow._commands[comflow._curpos]:
            nd.attrib['current'] = '1'
        ET.SubElement(nd, "TITLE").text = c.method_code()
        ET.SubElement(nd, "LINE").text = c.opt_line()
        if c.get_comment() != "":
            ET.SubElement(nd, "COMMENT").text = c.get_comment()


# ---------------------- Framework
def write_framework(fw, xmlnode):
    ' writes framework data to xmlnode/STATE'
    #Grids
    gnode = ET.SubElement(xmlnode, "GRIDS")
    for k, v in fw.grids2.items():
        nd = ET.SubElement(gnode, "GRID2")
        nd.attrib["name"] = k
        v.xml_save(nd)
    #User contours
    cnode = ET.SubElement(xmlnode, "CONTOURS")
    for k, v in fw.contours2.items():
        nd = ET.SubElement(cnode, "CONTOUR2")
        nd.attrib["name"] = k
        v.xml_save(nd)
    #boundary types
    bnode = ET.SubElement(xmlnode, "BTYPES")
    fw.boundary_types.xml_save(bnode)

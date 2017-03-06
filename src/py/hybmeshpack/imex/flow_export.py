import xml.etree.ElementTree as ET
from hybmeshpack import progdata
import hybmeshpack.hmcore.hmxml as hmxml
import native_export as natex


def xmlindent(elem, level=0):
    """ http://effbot.org/zone/element-lib.htm#prettyprint.
        It basically walks your tree and adds spaces and newlines so the tree
        is printed in a nice way
    """
    tabsym = "  "
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


def _root_xml():
    'returns root xml node tagged HybMeshProject'
    ret = ET.Element('HybMeshData')
    ret.attrib['ver'] = str(progdata.HybMeshVersion.current())
    return ret


def writexml(r, filename, formating=True):
    if formating:
        xmlindent(r)
    tree = ET.ElementTree(r)
    tree.write(filename, xml_declaration=True, encoding='utf-8')


# ------------------- Commands
def write_command_flow(comflow, xmlnode, flowname=None):
    """ creates xmlnode/FLOW and writes flow commands there """
    rt = ET.SubElement(xmlnode, "FLOW")
    if flowname is not None:
        rt.attrib['name'] = flowname
    # Commands list
    cl = ET.SubElement(rt, "COMMANDS")
    for c in comflow._commands:
        nd = ET.SubElement(cl, "COM")
        nd.attrib['name'] = c.method_code()
        if c is comflow._commands[comflow._curpos]:
            nd.attrib['current'] = '1'
        ET.SubElement(nd, "LINE").text = c.opt_line()
        if c.get_comment() != "":
            ET.SubElement(nd, "COMMENT").text = c.get_comment()


# ---------------------- Framework
def write_framework(fw, xmlnode):
    ' writes framework non-geometry data to xmlnode/STATE'
    zt = fw.get_zone_types()
    for k in sorted(zt.keys()):
        v = zt[k]
        nd = ET.SubElement(xmlnode, "BTYPE")
        nd.attrib["index"] = str(k)
        ET.SubElement(nd, "NAME").text = str(v)


# ---------------------- Everything
def flow_and_framework_tofile(filename, comflow, fmt="ascii"):
    'writes flow.CommandFlow object to xml file'
    cb = comflow.interface.ask_for_callback()
    cb1 = cb.subcallback(0, 3)

    cb1.pycall("Writing command flow", "", 0, 0)
    # write everithing except geometry data data
    r = _root_xml()
    write_command_flow(comflow, r)
    st = ET.SubElement(r.find('FLOW'), "STATE")
    write_framework(comflow.receiver, st)
    cb1.pycall("Writing command flow", "", 0.8, 0)
    writexml(r, filename, False)
    cb1.pycall("Writing command flow", "Done", 1, 1)

    # add geometry data and write to file
    doc, root = 0, 0
    try:
        cb2 = cb.subcallback(1, 3)
        doc, root = hmxml.open_doc(filename)
        [nd] = hmxml.query(root, "FLOW/STATE", required="=1")
        natex.export_all(doc, nd, comflow.receiver, fmt, cb2)
    except:
        raise
    finally:
        cb3 = cb.subcallback(2, 3)
        cb3.pycall("Write to file", "", 0.0, 0)
        hmxml.doc_to_file(doc, filename)
        hmxml.close_doc(doc, [root, nd])
        cb3.pycall("Write to file", "Done", 1, 1)

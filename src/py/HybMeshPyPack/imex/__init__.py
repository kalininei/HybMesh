import xml.etree.ElementTree as ET
import writexml
import readxml
import gridexport
from HybMeshPyPack import com
from HybMeshPyPack import gdata


def write_flow_and_framework_to_file(comflow, filename):
    'writes flow.CommandFlow comflow to xml file'
    r = writexml._root_xml()
    writexml.write_command_flow(comflow, r)
    writexml.write_framework(comflow.get_receiver(), r.find('FLOW'))
    writexml.writexml(r, filename)


def read_flow_and_framework_from_file(filename):
    '-> CommandFlow'
    root = ET.parse(filename).getroot()
    flow_nodes = root.findall('.//FLOW')
    if len(flow_nodes) == 0:
        raise Exception('No proper data in %s' % filename)
    flow_node = flow_nodes[0]
    f = com.flow.CommandFlow()
    readxml.load_command_flow(f, flow_node)
    state_node = flow_node.find('STATE')
    if state_node is not None:
        data = gdata.Framework()
        readxml.load_framework_state(data, state_node)
    else:
        data = gdata.Framework()
    f.set_receiver(data)
    return f


def export_grid(fmt, fn, name, fw=None, flow=None):
    """ exports grid from framework fw or flow receiver to
        filename fn using format fmt. Possible formats:
            vtk, hmg, msh, ggen, gmsh
    """
    # 1. Find grid
    try:
        if fw is None:
            fw = flow.get_receiver()
        _, _, grid = fw.get_grid(name=name)
    except:
        raise Exception('Can not find grid for exporting')
    # 2. Export regarding to format
    if fmt == 'vtk':
        gridexport.vtk(grid, fn)
    elif fmt == 'hmg':
        gridexport.hmg(grid, fn)
    elif fmt == 'msh':
        gridexport.msh(grid, fn, fw.boundary_types)
    elif fmt == 'ggen':
        gridexport.ggen(grid, fn)
    elif fmt == "gmsh":
        gridexport.gmsh(grid, fn, fw.boundary_types)
    else:
        raise Exception('Unknown grid format %s' % fmt)

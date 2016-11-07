import xml.etree.ElementTree as ET
import writexml
import readxml
import gridexport
import contexport
import grid3export
from hybmeshpack import com
from hybmeshpack import gdata
from hybmeshpack.basic.interf import Callback
from hybmeshpack.gdata import contour2


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


def export_grid(fmt, fn, name, fw=None, flow=None, adata=None):
    """exports grid from framework fw or flow receiver to
    filename fn using format fmt. Possible formats:
        vtk, hmg, msh, ggen, gmsh, tecplot,
        vtk3d, msh3d

    adata - additional data which dependes on format:
      * msh, msh3d: define periodic conditions

        msh   - [btype_per, btype_shadow, is_reversed,
                 .....]
        msh3d - [btype_per, btype_shadow, Point btype_per, Point btype_shadow,
                 .....]
    """
    # 1. Find grid
    try:
        if fw is None:
            fw = flow.get_receiver()
        if (fmt[-2:] == '3d'):
            _, _, grid = fw.get_grid3(name=name)
        else:
            _, _, grid = fw.get_grid(name=name)
        if flow is not None:
            callb = flow.get_interface().ask_for_callback(Callback.CB_CANCEL2)
        else:
            callb = None
    except:
        raise Exception('Can not find grid for exporting')
    # 2. Export regarding to format
    if fmt == 'vtk':
        gridexport.vtk(grid, fn)
    elif fmt == 'hmg':
        gridexport.hmg(grid, fn)
    elif fmt == 'msh':
        gridexport.msh(grid, fn, fw.boundary_types, adata)
    elif fmt == 'ggen':
        gridexport.ggen(grid, fn)
    elif fmt == "gmsh":
        gridexport.gmsh(grid, fn, fw.boundary_types)
    elif fmt == "tecplot":
        gridexport.tecplot(grid, fn, fw.boundary_types)
    elif fmt == "vtk3d":
        grid3export.vtk(grid, fn, callb)
    elif fmt == "msh3d":
        grid3export.msh(grid, fn, fw.boundary_types, callb, adata)
    elif fmt == "tecplot3d":
        grid3export.tecplot(grid, fn, fw.boundary_types, callb)
    else:
        raise Exception('Unknown grid format %s' % fmt)


def export_contour(fmt, fn, name, fw=None, flow=None):
    """exports contour from framework 'fw' or 'flow' receiver to
    filename fn using format fmt. Possible formats:
        vtk, tecplot, hmc
    """
    # Find contour
    try:
        if fw is None:
            fw = flow.get_receiver()
        if isinstance(name, list):
            cont = contour2.Contour2()
            for nm in name:
                try:
                    _, _, c = fw.get_ucontour(name=nm)
                except KeyError:
                    c = fw.get_grid(name=nm)[2].cont
                cont.add_from_abstract(c)
        else:
            try:
                _, _, cont = fw.get_ucontour(name=name)
            except KeyError:
                cont = fw.get_grid(name=name)[2].cont
    except:
        raise Exception('Can not find contour for exporting')
    # 2. Export regarding to format
    if fmt == 'vtk':
        contexport.vtk(cont, fn)
    elif fmt == 'hmc':
        contexport.hmg(cont, fn)
    elif fmt == 'tecplot':
        contexport.tecplot(cont, fn, fw.boundary_types)
    else:
        raise Exception('Unknown contour format %s' % fmt)


def export_grid3_surface(fmt, fn, name, fw=None, flow=None):
    # Find grid
    try:
        if fw is None:
            fw = flow.get_receiver()
        _, _, grid = fw.get_grid3(name=name)
        if flow is not None:
            callb = flow.get_interface().ask_for_callback(Callback.CB_CANCEL2)
        else:
            callb = None
    except:
        raise Exception('Can not find contour for exporting')
    # 2. Export regarding to format
    if fmt == 'vtk':
        grid3export.vtk_surface(grid, fn, callb)
    else:
        raise Exception('Unknown format %s' % fmt)

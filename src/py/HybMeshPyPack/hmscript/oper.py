from HybMeshPyPack import com, basic
from HybMeshPyPack.hmscript import flow
import HybMeshPyPack.com.gridcom
from HybMeshPyPack.basic.geom import Point2


def ExcludeContours(grid, conts, exclude_outer=True):
    """ Builds a new grid by excluding contour area from existing grid.
        grid - source grid identifier
        conts - contour or list of contours identifiers for exclusion
        exclude_outer - exclude inner or outer region of contours
        returns new grid identifier
    """
    if not isinstance(conts, list):
        conts = [conts]
    c = com.gridcom.ExcludeContours({"grid_name": grid,
        "cont_names": conts,
        "is_inner": not exclude_outer})
    flow.exec_command(c)
    return c._get_added_names()[0][0]


def UniteGrids(base_grid, imp_grids, empty_holes=False, fix_bnd=False):
    """ Makes grids impositions
        base_grid - basic grid identifier
        imp_grids - list of grids for imposition as [(grid_id, buffer), () ...]
            where grid_id is imposed grid identifier,
                  buffer - size of the buffer for current imposition
            Each next grid will be imposed on the result of previous imposition
        empty_holes - keep all empty zone in imposed grids
        fix_bnd - whether to fix all boundary nodes
        Returns identifier of the newly created grid
    """
    args = {"base": base_grid, "empty_holes": empty_holes,
            "fix_bnd": fix_bnd, "plus": []}
    for ig in imp_grids:
        args["plus"].append({"name": ig[0], "buf": ig[1], "den": 7})
    c = com.gridcom.UniteGrids(args)
    flow.exec_command(c)
    return c._get_added_names()[0][0]


class BoundaryGridOption(object):
    def __init__(self, contour_id=None,
                partition=[0],
                direction="left",
                bnd_stepping="no",
                bnd_step=0.1,
                range_angles=[30, 135, 225, 275],
                force_conformal=False,
                start_point=None,
                end_point=None):
        """ Create boundary grid options object
            contour_id - string id of user or grid contour
            partition - partition in perpendicular direction. List of
                        ascending floats starting with zero
                        which represents the distance from contour to grid
                        layer
            direction - "left"/"right". Grid will be built to the left/right
                        from the contour (using positive tracing)
            bnd_stepping - partition along the contour.
                        "no" - no artificial stepping. Only contour edges will
                               be used as grid nodes.
                        "const" - use artificial stepping and ignore contour
                               edges
                        "keep_shape" - use stepping and keep significant
                               contour edges
                        "keep_all" - use stepping and keep all contour edges
            bnd_step - float size of artificial stepping along the contour
            range_angles - list of 4 angle (deg) values which define algorithms
                               for corners treatment:
                           [0,     ra[0]]: acute angle algo
                           [ra[0], ra[1]]: right angle algo
                           [ra[1], ra[2]]: straight angle algo
                           [ra[2], ra[3]]: reentrant angle algo
                           [ra[3], 360]:  round algo
            start_point, end_point - points in [x, y] format which define
                 segment of the contour for building grid.
                 If Both are None -> whole contour will be used.
            force_conformal - if true then algorithm will use
                 strictly conformal mappings
        """
        self.contour_id = contour_id
        self.partition = partition
        self.direction = direction
        self.bnd_stepping = bnd_stepping
        self.bnd_step = bnd_step
        self.range_angles = range_angles
        self.force_conformal = force_conformal
        self.start_point = start_point
        self.end_point = end_point


def BuildBoundaryGrid(opts):
    """ Builds a boundary grid near contour
        opts - BoundaryGridOption objects (list or single object)
    """
    inp = []
    if not isinstance(opts, list):
        opts2 = [opts]
    else:
        opts2 = opts
    for op in opts2:
        d = {}
        d['source'] = op.contour_id
        d['partition'] = op.partition
        d['direction'] = op.direction
        if op.bnd_stepping == "no":
            d['mesh_cont'] = 0
        elif op.bnd_stepping == "const":
            d['mesh_cont'] = 3
        elif op.bnd_stepping == "keep_shape":
            d['mesh_cont'] = 2
        else:
            d['mesh_cont'] = 1
        d['mesh_cont_step'] = op.bnd_step
        d['algo_acute'] = op.range_angles[0]
        d['algo_right'] = op.range_angles[1]
        d['algo_straight'] = op.range_angles[2]
        d['algo_reentr'] = op.range_angles[3]
        if (op.start_point is not None and op.end_point is not None):
            d['start'] = Point2(*op.start_point)
            d['end'] = Point2(*op.end_point)
        else:
            d['start'] = d['end'] = Point2(0, 0)
        d['force_conf'] = op.force_conformal
        inp.append(d)

    c = com.gridcom.BuildBoundaryGrid({"opt": inp})
    flow.exec_command(c)
    return c._get_added_names()[0][0]

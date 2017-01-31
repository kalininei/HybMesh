import hybmeshpack.hmcore.g2 as g2core


def _parse(self, fn):
    """ ->
        {bindex: bname}
        [ ... points ... ],
        [ ... 2d cells ... ],
        [ [p1, p2, btype], ... ]
    """
    r1, r2, r3, r4 = {}, [], [], []
    r2m, r3m = {}, {}
    dt = file(fn).readlines()
    dt = map(lambda x: x.strip(), dt)
    dt = [item for item in dt if len(item) > 0]
    i = 0
    while i < len(dt):
        # 1
        if dt[i] == "$PhysicalNames":
            k = int(dt[i + 1])
            for d in dt[i + 2:i + 2 + k]:
                ds = d.split(' ', 2)
                r1[int(ds[1])] = ds[2].strip('"')
            i += k + 1
        # 2
        if dt[i] == "$Nodes":
            k = int(dt[i + 1])
            for d in dt[i + 2:i + 2 + k]:
                ds = d.split()
                r2m[int(ds[0]) - 1] = (float(ds[1]), float(ds[2]))
            i += k + 1
        # 3
        if dt[i] == "$Elements":
            k = int(dt[i + 1])
            for d in dt[i + 2:i + 2 + k]:
                ds = map(int, d.split())
                index = ds[0] - 1
                if ds[1] == 1:
                    # edge
                    tp = ds[3] if ds[2] > 0 else 0
                    if tp > 0:
                        r4.append([ds[-2] - 1, ds[-1] - 1, tp])
                elif ds[1] == 2:
                    # triangle cell
                    r3m[index] = [ds[-3] - 1, ds[-2] - 1, ds[-1] - 1]
                elif ds[1] == 3:
                    # quad cell
                    r3m[index] = [ds[-4] - 1, ds[-3] - 1,
                                  ds[-2] - 1, ds[-1] - 1]
                elif ds[1] == 15:
                    # ignore point cell
                    pass
                else:
                    raise Exception("Not supported gmsh element type")
            i += k + 1
        i += 1

    # point renumbering for successive indexing starting from zero
    need_point_renumbering = len(r2m) != max(r2m.keys()) + 1
    if need_point_renumbering:
        oldnew = {}
        i = 0
        for ind in sorted(r2m.keys()):
            oldnew[ind] = i
            r2.append(r2m[ind])
            i += 1
        for a in r4:
            a[0] = oldnew[a[0]]
            a[1] = oldnew[a[1]]
        for a in r3m.values():
            for i in range(len(a)):
                a[i] = oldnew[a[i]]
    else:
        for ind in sorted(r2m.keys()):
            r2.append(r2m[ind])

    for ind in sorted(r3m.keys()):
        r3.append(r3m[ind])

    return r1, r2, r3, r4


def _parse_bt(btnames, allbt):
    """ {boundary index: boundary name} """
    ret = {}
    for b in allbt:
        if b in btnames:
            ret[b] = btnames[b]
        else:
            nm = '-'.join(["gmsh", "boundary", str(b)])
            ret[b] = nm
    return ret


def grid2(self, fn):
    """ ->( Grid2.grid.cdata, {bindex: bname, ...}
    """
    bt, pts, cls, bedges = _parse(fn)
    abt = set()
    for b in bedges:
        abt.add(b[2])
    usedbt = _parse_bt(bt, abt)
    ret = g2core.from_points_cells(pts, cls, bedges)
    return ret, usedbt

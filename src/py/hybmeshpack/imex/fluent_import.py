import re
from hybmeshpack.hmcore import g2 as g2core


def _parse_hex(txt):
    """ -> [list of integers] from text """
    return [int(x, 16) for x in txt.split()]


def _parse_floats(txt):
    """ -> [list of floats] """
    return map(float, txt.split())


def _parse_parant(txt):
    """ returns list of strings within first level paranthesis """
    ret = []
    parant = {}
    f = -1
    while 1:
        f = txt.find('(', f + 1)
        if f == -1:
            break
        parant[f] = '('
    while 1:
        f = txt.find(')', f + 1)
        if f == -1:
            break
        parant[f] = ')'

    oparant = []
    for k in sorted(parant.keys()):
        oparant.append((k, parant[k]))
    st = -1
    counter = 0
    for p in oparant:
        if st == -1 and p[1] == '(':
            st = p[0]
            counter = 0
        else:
            if p[1] == '(':
                counter = counter + 1
            elif p[1] == ')':
                if counter > 0:
                    counter = counter - 1
                else:
                    ret.append(txt[st + 1:p[0]])
                    st = -1
    return ret


def _parse(fn):
    """ -> [(x0, y0), ....] - points list
           [p0, p1, cell left, cell right, btype] - edges connectivity
           [bindex: b name] - boundaries
    """
    dt = file(fn).read()
    dt = _parse_parant(dt)
    dtp = {}
    for d in dt:
        k = re.search('\d+', d).end(0)
        w1, w2 = d[:k], d[k + 1:]
        w1 = int(w1)
        if w1 != 0:
            if w1 not in dtp:
                dtp[w1] = []
            dtp[w1].append(w2)
    # check dimension
    if (int(dtp[2][0]) != 2):
        raise Exception("Invalid msh file dimenstion")

    # read nodes
    points = []
    for plines in dtp[10]:
        ln = _parse_parant(plines)
        if len(ln) < 2:
            continue
        info = _parse_hex(ln[0])
        index0, index1 = info[1], info[2]
        if (len(points) < index1):
            points.extend([None] * (index1 - len(points)))
        coords = _parse_floats(ln[1])
        for i in range(index1 - index0 + 1):
            p = [coords[2 * i], coords[2 * i + 1]]
            points[i + index0 - 1] = p

    # read edges as [p0, p1, cell_right, cell_left, zone_type]
    eds = []
    for clines in dtp[13]:
        ln = _parse_parant(clines)
        info = _parse_hex(ln[0])
        if info[0] == 0:
            continue
        index0, index1 = info[1], info[2]
        if (len(eds) < index1):
            eds.extend([None] * (index1 - len(eds)))
        conn = _parse_hex(ln[1])
        for i in range(index1 - index0 + 1):
            if info[4] == 2:
                line = conn[4 * i:4 * i + 4]
            elif info[4] == 0:
                line = conn[5 * i + 1: 5 * i + 5]
            eds[i + index0 - 1] = [line[0] - 1,
                                   line[1] - 1,
                                   line[2] - 1,
                                   line[3] - 1,
                                   info[0]]

    # read boundary names
    ztps = {}
    if 45 not in dtp and 39 not in dtp:
        raise Exception("ZONE section (39, 45) was not found")
    if 45 not in dtp:
        dtp[45] = []
    if 39 in dtp:
        dtp[45] += dtp[39]
    for zline in dtp[45]:
        w = _parse_parant(zline)[0].split()
        if w[1] != "interior":
            ztps[int(w[0])] = w[2]
    # filter out bc which are not present
    ztps2 = {}
    for e in eds:
        if e[4] in ztps:
            if e[4] not in ztps2:
                ztps2[e[4]] = ztps[e[4]]
        else:
            e[4] = 0  # set interior zone type to zero
    ztps = ztps2

    return points, eds, ztps


def grid2(fn):
    """ Grid2.grid2.cdata, {bindex: bname} """
    p, e, b = _parse(fn)
    return g2core.from_points_edges(p, e), b

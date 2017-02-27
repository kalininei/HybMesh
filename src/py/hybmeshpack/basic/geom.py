#!/usr/bin/env python
"basic geometry"
import math


def angle_3pnt(p1, p2, p3):
    '-> [0, 2*pi]. Get the angle from 3 consecutive points'
    v1x, v1y = p1[0] - p2[0], p1[1] - p2[1]
    v2x, v2y = p3[0] - p2[0], p3[1] - p2[1]
    dot = v1x * v2x + v1y * v2y
    cross = v1x * v2y - v1y * v2x
    a = math.atan2(cross, dot)
    if a < 0:
        a += 2 * math.pi
    return a


def build_seg(x0, x1, nx, cust):
    """ modifies cust float list so that:
          if len(cust) = 1 -> uniform partition of [x0, x1] with step=cust[0]
          if len(cust) = 0 -> uniform partition of [x0, x1] with nx steps
    """
    if x1 < x0:
        x0, x1 = x1, x0
    if len(cust) > 1:
        return
    if len(cust) == 1:
        nx = max(1, int(round((x1 - x0) / cust[0])))
        del cust[0]
    for i in range(nx + 1):
        cust.append(x0 + (x1 - x0) * float(i) / nx)


def div_range(a, b, n, k=1):
    """ -> [a, ...., b] - list of floats

    Divides [a, b] range into n subsections.
    k is the refinement coefficient which works like:
        h[i] = k * h[i-1]
    """
    if a != 0:
        return [x + a for x in div_range(0, b - a, n, k)]
    # now a = 0, b = length of the section
    st = [k ** i for i in range(n)]
    a0 = float(b) / sum(st)
    ret = [0]
    for x in st:
        ret.append(ret[-1] + a0 * x)
    return ret

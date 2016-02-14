from HybMeshPyPack import hmscript  # NOQA

b1 = hmscript.AddBoundaryType(1, "Wall1")
b2 = hmscript.AddBoundaryType(2, "Wall2")
b3 = hmscript.AddBoundaryType(3, "Wall3")
b4 = hmscript.AddBoundaryType(4, "Wall4")


c1 = hmscript.AddRectCont([-10, -10], [15, 15], [b1, b2, b3, b4])
g1 = hmscript.AddUnfRectGrid([0, 0], [10, 10], 4, 5)
c2 = hmscript.GridBndToCont(g1, True)


def btfun(x0, y0, x1, y1, b):
    global b1, b2, b3, b4
    if y0 == 0 and y1 == 0:
        return b1
    if y0 == 10 and y1 == 10:
        return b3
    if x0 == 0 and x0 == 0:
        return b2
    if x0 == 10 and x1 == 10:
        return b4
hmscript.SetBTypeToContour(g1, bfun=btfun)

c3 = hmscript.CreateContour([[0, 0], [2, 0], [3, 0], [8, 0.1], [4, 4],
    [0, 0]], b2)
c4 = hmscript.SimplifyContour(c3, simplify=True, angle=10)
c5 = hmscript.UniteContours([c1, c4])


def btfun(x0, y0, x1, y1, b):
    global b1, b2
    return b1 if x0 > 5 else b2

hmscript.SetBTypeToContour(c5, bfun=btfun)
c6 = hmscript.ImportContourASCII("asciicont.txt", force_closed=True)
hmscript.ExportHMG(g1, "hmggrid.hmg")
hmscript.ExportMSH(g1, "hmggrid.msh")
hmscript.ExportGMSH(g1, "gmshgrid.msh")
g2 = hmscript.ImportGridHMG("hmggrid.hmg")
g3 = hmscript.ExcludeContours(g2, [c5], exclude_outer=True)
g4 = hmscript.AddUnfRectGrid([-5, -5], [5, 5], 20, 20)
g5 = hmscript.AddUnfCircGrid([0, 0], 1, 30, 10, 1.2, True)
g6 = hmscript.UniteGrids(g4, [(g5, 1.0)])
o1 = hmscript.BoundaryGridOption(g4)
o1.partition = [0, 0.1, 0.3, 0.6]
o1.direction = -1
o1.bnd_stepping = "keep_shape"
o1.bnd_step = 0.1
g7 = hmscript.BuildBoundaryGrid(o1)
g8 = hmscript.ImportGridMSH("hmggrid.msh")
g9 = hmscript.ImportGridGMSH("gmshgrid.msh")
hmscript.SaveProject("scr1.hmp")

from HybMeshPyPack import hmscript  # NOQA

g1 = hmscript.ImportGridMSH("testgambit.msh")
hmscript.ExportMSH(g1, "exmsh1.msh")
hmscript.ExportGMSH(g1, "exmsh2.msh")

hmscript.SaveProject("scr2.hmp")

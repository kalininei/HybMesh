import sys
import os.path

try:
    # this is used for installed version of HybMesh
    # config_installed is generated during installation procedure
    import config_installed as config
except Exception:
    import config_dbg as config

# preload all libraries in certain order for windows runs
# this is needed to avoid errors due to implicit library loads
if sys.platform.startswith('win'):
    _liblist = [
        'libgpc.dll',
        'libwinpthread-1.dll',
        'libgcc_s_seh-1.dll',
        'libquadmath-0.dll',
        'libgfortran-3.dll',
        'libgomp-1.dll',
        'libstdc++-6.dll',
        'libiconv-2.dll',
        'libxml2-2.dll',
        'libblas.dll',
        'liblapack.dll',
        'libmetis.dll',
        'libsuitesparseconfig.4.5.1.dll',
        'libcamd.2.4.4.dll',
        'libccolamd.2.9.4.dll',
        'libcolamd.2.9.4.dll',
        'libamd.2.4.4.dll',
        'libcholmod.3.0.9.dll',
        'libspqr.2.0.5.dll',
        'libtetgen.dll',
        'libGmshLocal.dll',
        'libpolyclipping.dll',
        'libhmproject.dll',
        'libbgeom2d.dll',
        'libhmmath.dll',
        'libhybmesh_contours2d.dll',
        'libhybmesh_surfaces3d.dll',
        'libhmgrid2d.dll',
        'libhmgrid3d.dll',
        'libhmnumeric.dll',
        'libcrossgrid.dll',
        'libhmmapping.dll',
        'libhmblay.dll',
        'libhmcport.dll',
    ]
    import ctypes as ct
    for name in _liblist:
        try:
            ct.cdll.LoadLibrary(os.path.join(config.libdir, name))
        except Exception as e:
            sys.exit("error in library linking: %s, with message %s" %
                     (name, str(e)))

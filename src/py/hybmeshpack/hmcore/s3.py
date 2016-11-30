import ctypes as ct  # NOQA
from . import libhmcport


def extract_grid3_surface(c_g3):
    "-> c pointer to deep copied surface"
    ret = libhmcport.extract_grid3_surface(c_g3)
    if ret == 0:
        raise Exception("Failed to extract grid surface")
    return ret


# s3core.surface_from_hmxml(reader, c_s3reader[i])

from hybmeshpack.hmcore import s3 as s3core
import hybmeshpack.hmcore as hmcore


def hmc(surfs, names, fname, fmt, wr=None):
    c_writer, c_swriter = 0, 0
    try:
        c_writer = hmcore.hmxml_new() if wr is None else wr
        for s, nm in zip(surfs, names):
            # write
            c_swriter = s3core.swriter_create(nm, s.cdata, c_writer,
                                              c_writer, fmt)
            # free data
            s3core.free_swriter(c_swriter) if c_swriter != 0 else None
            c_swriter = 0
    except:
        raise
    finally:
        s3core.free_swriter(c_swriter) if c_swriter != 0 else None
        if wr is None and c_writer != 0:
            hmcore.hmxml_finalize(c_writer, fname)

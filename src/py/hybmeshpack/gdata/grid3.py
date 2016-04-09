class Grid3(object):
    def __init__(self):
        # pointer to data stored at c-side
        self.cpnt = None

    def __del__(self):
        if self.cpnt:
            import hmcore.g3
            hmcore.g3.free_g3(self.cdata)

    @staticmethod
    def from_cdata(cdata):
        ret = Grid3()
        ret.cdata = cdata
        return ret

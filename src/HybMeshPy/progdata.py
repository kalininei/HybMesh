try:
    import config_installed as config
except:
    import config


class HybMeshVersion:
    def __init__(self, s):
        """initializes version for string of "1.2.3" type"""
        def __only_nums(n):
            return ''.join(k for k in n if k.isdigit())
        self.nums = map(__only_nums, s.split('.'))
        self.nums = map(int, [n for n in self.nums if len(n) > 0])

    def __str__(self):
        return '.'.join(map(str, self.nums))

    @classmethod
    def str_compare(cls, v1, v2):
        '(str, str) -> [-1 or 0 or 1]. Compares v1 and v2 string'
        [v1, v2] = [HybMeshVersion(v) for v in [v1, v2]]
        for i in range(min(len(v1.nums), len(v2.nums))):
            if v1.nums[i] < v2.nums[i]:
                return -1
            elif v1.nums[i] > v2.nums[i]:
                return 1
        return 0


def check_for_updates():
    '->(str, str, -1/0/1) current version, latest version, what is larger'
    print "DUMMY updates check"
    return (config.version, None, None)


def program_version():
    '-> HybMeshVersion'
    return HybMeshVersion(config.version)


def project_url():
    '->str'
    return config.project_url


def get_lib_fn(string_id):
    """ string_id:
            crossgrid,
            hybmesh_contours2d
    """
    import glob
    if string_id == 'crossgrid':
        c = glob.glob("%s/*crossgrid*" % config.libdir)
        if len(c) == 0:
            raise Exception("libcrossgrid was not found at %s"
                    % config.libdir)
        return c[0]
    elif string_id == 'hybmesh_contours2d':
        c = glob.glob("%s/*hybmesh_contours2d*" % config.libdir)
        if len(c) == 0:
            raise Exception("libhybmesh_contours2d was not found at %s"
                    % config.libdir)
        return c[0]
    else:
        raise Exception("Unknown library %s" % string_id)





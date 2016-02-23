try:
    # this is used for installed version of HybMesh
    # config_installed is generated during installation procedure
    import config_installed as config
except Exception:
    import config as config


class HybMeshVersion:
    def __init__(self, s):
        """initializes version for string of "1.2.3" type"""
        try:
            def only_nums(n):
                return ''.join(k for k in n if k.isdigit())

            self.nums = map(only_nums, s.split('.'))
            self.nums = map(int, [n for n in self.nums if len(n) > 0])
            if len(self.nums) != 3:
                raise
        except Exception:
            raise Exception("Invalid version string: %s" % str(s))

    def __str__(self):
        return '.'.join(map(str, self.nums))

    def __cmp__(self, other):
        if not isinstance(other, HybMeshVersion):
            raise Exception("Invalid comparison")
        return cmp(self.nums, other.nums)

    @staticmethod
    def last():
        import urllib2
        import json
        try:
            response = urllib2.urlopen(config.last_ver_url, timeout=1)
            r = json.loads(response.read())
            return HybMeshVersion(r['tag_name'])
        except:
            return None

    @staticmethod
    def current():
        '-> HybMeshVersion'
        return HybMeshVersion(config.version)


def check_for_updates():
    '->(str, str, -1/0/1): current version, latest version, cur<lat / == / >'
    cv = HybMeshVersion.current()
    lv = HybMeshVersion.last()
    if lv is not None:
        cmpres = cmp(cv, lv)
        lvs = str(lv)
    else:
        cmpres = None
        lvs = None
    return (str(cv), lvs, cmpres)


def project_url():
    '-> github url of the hybmesh project'
    return config.proj_url


def get_lib_fn(string_id):
    """ string_id:
            cport
    """
    import glob
    if string_id == 'cport':
        c = glob.glob("%s/*hmcport*" % config.libdir)
        if len(c) == 0:
            raise Exception("libhmcport was not found at %s"
                            % config.libdir)
        return c[0]
    else:
        raise Exception("Unknown library %s" % string_id)

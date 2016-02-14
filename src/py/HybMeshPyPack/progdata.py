try:
    # this is used for installed version of HybMesh
    # config_installed is generated during installation procedure
    import config_installed as config
except Exception:
    import config as config

  
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
                return 1
            elif v1.nums[i] > v2.nums[i]:
                return -1
        return 0
    
    @staticmethod
    def last_version():
        import urllib2
        import json
        try:
            response = urllib2.urlopen(config.last_ver_url, timeout=1)
            r = json.loads(response.read())
            return r['tag_name']
        except:
            return None


def check_for_updates():
    '->(str, str, -1/0/1) current version, latest version, what is larger'
    lv = HybMeshVersion.last_version()
    if lv is not None: 
        cmp = HybMeshVersion.str_compare(config.version, lv)
    else:
        cmp = None
    return (config.version, lv, cmp)


def program_version():
    '-> HybMeshVersion'
    return HybMeshVersion(config.version)

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





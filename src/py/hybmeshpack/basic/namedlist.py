""" present NamedList class which is used for storage of all
geometric objects
"""
import copy


class _NamedListEntry(object):
    """ entry with object and its name in NamedList"""
    def __init__(self, name, obj):
        self.name = name
        self.obj = obj


class _StemsDict(object):
    def __init__(self):
        self._stems = {}  # {"stem": set of used indices}

    def newind(self, stem, ind):
        "adds an index to the stem and returns it"
        try:
            a = self._stems[stem]
            if ind is None and ind in a:
                ind = 1
            while ind in a:
                ind += 1
        except KeyError:
            pass
        self.forceadd(stem, ind)
        return ind

    def forceadd(self, stem, ind):
        if stem in self._stems:
            self._stems[stem].add(ind)
        else:
            self._stems[stem] = set({ind})

    def clear(self):
        self._stems.clear()


class NamedListPool(object):
    """ pool of unique names. Helper class for NamedList object. """
    def __init__(self):
        self._data = set()
        self._stems = _StemsDict()

    def _decompose(self, name):
        index = len(name) - 1
        # last non-digit index
        while index >= 0 and name[index].isdigit():
            index -= 1
        if index < 0:
            # all values are digits
            raise Exception("Invalid name")
        elif index == len(name) - 1:
            return name, None
        else:
            return name[:index + 1], int(name[index + 1:])

    def _analyze(self, name):
        """ returns stem and positive index for given name """
        nm, ind = self._decompose(name)
        if name in self._data:
            ind = self._stems.newind(nm, ind)
        self._stems.forceadd(nm, ind)
        return nm, ind

    def _add_new_name(self, name):
        """ adds a name to the pool. If name already exists modifies it
        by adding positive index to its end. Returns added name.
        """
        newname, ind = self._analyze(name)
        if ind is not None:
            newname = newname + str(ind)
        self._data.add(newname)
        return newname

    def _remove(self, name):
        """ removes name from the pool """
        self._data.remove(name)

    def clear(self):
        """ frees everything including addition history"""
        self._data.clear()
        self._stems.clear()

    def fillfrom(self, pool):
        """ deep copies all dictionaries from pool """
        self._data = copy.deepcopy(pool._data)
        self._stems = copy.deepcopy(pool._stems)


class NamedList(object):
    """ Represent a named list of objects. All objects and names are unique.
        Addition of new entry with already existed object leads to ValueError.
        Addition of new entry with already existed name leads to
        modification of the name by adding a positive index to it.
    """

    def __init__(self, prefix, pool):
        """ prefix (str) - default name, which will be used if no name is
                           defined;
            pool (NamedListPool) - pool of names which could be shared
                                   amoung multiple NamedLists
        """
        self._prefix = prefix
        self._names_pool = pool
        self._data = []  # list of _NamedListEntry

    def __str__(self):
        s = []
        for e in self._data:
            s.append("  " + e.name + ": " + str(e.obj))
        if len(s) == 0:
            return '{}'
        s[0] = '{' + s[0][1:]
        return '\n'.join(s) + '}'

    def _find_entry_by_name(self, name):
        for e in self._data:
            if e.name == name:
                return e
        else:
            raise KeyError

    def _find_entry_by_obj(self, obj):
        for e in self._data:
            if e.obj == obj:
                return e
        else:
            raise KeyError

    def _find_entry(self, obj=None, name=None):
        if obj is not None and name is None:
            return self._find_entry_by_obj(obj)
        elif name is not None and obj is None:
            return self._find_entry_by_name(name)
        elif name is not None and obj is not None:
            e = self._find_entry_by_name(name)
            if e.obj != obj:
                raise KeyError
            return e
        else:
            raise KeyError

    def _name_exists(self, name):
        try:
            self._find_entry_by_name(name)
            return True
        except KeyError:
            return False

    def _obj_exists(self, obj):
        try:
            self._find_entry_by_obj(obj)
            return True
        except KeyError:
            return False

    def _set_name(self, name):
        name2 = name if name is not None else self._prefix
        name3 = self._names_pool._add_new_name(name2)
        return name3

    def _remove_name(self, name):
        self._names_pool._remove(name)

    def append(self, obj, name=None):
        """ adds obj with name to the end of data list.
            If name==None predifined prefix will be used.
            Returns added name.
        """
        if self._obj_exists(obj):
            raise ValueError
        name3 = self._set_name(name)
        entry = _NamedListEntry(name3, obj)
        self._data.append(entry)
        return name3

    def get(self, obj=None, name=None):
        """ get and object by pointer of name. -> (obj, name).
            Raises KeyError if fails.
            If both obj and name are defined checks whether
            name corresponds to obj and raises on fail.
        """
        e = self._find_entry(obj, name)
        return e.obj, e.name

    def remove(self, obj=None, name=None):
        " removes object by value or name from data"
        e = self._find_entry(obj, name)
        self._remove_name(e.name)
        self._data.remove(e)

    def change(self, oldobj=None, oldname=None, newobj=None, newname=None):
        "changes object name or value"
        e = self._find_entry(oldobj, oldname)
        oldobj = e.obj
        oldname = e.name
        if newobj is not None and oldobj != newobj:
            if self._obj_exists(newobj):
                raise ValueError
            e.obj = newobj
        if newname is not None and newname != oldname:
            self._remove_name(oldname)
            e.name = self._set_name(newname)

    def obj_by_name(self, name):
        " gives object by its name"
        e = self._find_entry(None, name)
        return e.obj

    def name_by_obj(self, obj):
        " gives name by object pointer"
        e = self._find_entry(obj, None)
        return e.name

    def all_obj(self):
        " -> [list of all objects] "
        return [a.obj for a in self._data]

    def all_names(self):
        " -> [list of all object names]"
        return [a.name for a in self._data]

    def all_data(self):
        " -> [(obj, name)] for each entry "
        return [(a.obj, a.name) for a in self._data]

    def len(self):
        " number of data entries "
        return len(self._data)

    def clear(self):
        """ removes all data and all respective names from the pool.
            For complete removal of all history call clear() procedure
            for pool object explicitly after this call.
        """
        for e in self._data:
            self._names_pool._remove(e.name)
        del self._data[:]

    def shallow_fillfrom(self, nlist):
        """ Copies data from nlist. Do not affect pool and prefix data"""
        del self._data[:]
        for e in nlist._data:
            self._data.append(_NamedListEntry(e.name, e.obj))

import h5py
import scipy as sp

class hdf5data():
    """Minimal class dummy"""
    pass


def get_data(IN):
    """This is a wrapper function that uses low level h5py to dereference
       data in a (most likely) matlab generated HDF5."""
    
    ### do we need to go deeper?
    if isinstance(IN, h5py._hl.group.Group):
        ret = hdf5data() 
        for k in IN:
            setattr(ret, k, get_data(IN[k]))
    elif isinstance(IN, h5py.h5g.GroupID):
        ret = hdf5data()
        for k in IN:
            dset = h5py.h5o.open(IN, k)
            tmp = sp.empty(dset.shape, dtype=dset.dtype)
            dset.read(h5py.h5s.ALL, h5py.h5s.ALL, tmp)
            setattr(ret, k, tmp.T)
    else:
        ret = sp.empty(shape=IN.shape, dtype='object')
        if len(IN.shape) > 1:
            for i in xrange(IN.shape[0]):
                for j in xrange(IN.shape[1]):
                    ref = IN[i, j]
                    dref = h5py.h5r.dereference(ref, IN._id)
                    if isinstance(dref, h5py.h5g.GroupID):
                        ret[i, j] = get_data(dref)
                    else:
                        ret[i, j] = sp.empty(dref.shape, dtype=dref.dtype)
                        dref.read(h5py.h5s.ALL, h5py.h5s.ALL, ret[i, j])
                        ret[i, j] = ret[i, j].T
        else:
            for i in xrange(IN.shape[0]):
                ref = IN[i]
                dref = h5py.h5r.dereference(ref, IN._id)
                ret[i] = sp.empty(dref.shape, dtype=dref.dtype)
                dref.read(h5py.h5s.ALL, h5py.h5s.ALL, ret[i])
                ret[i] = ret[i].T
    return ret
        

def loadmat(fname):
    """Replacement for the scipy.io loadmat function in case matlab has 
       stored the file in HDF5 format due to total size."""
    
    ret = dict()

    IN = h5py.File(fname, 'r')
    for k in IN:
        if k == '#refs#':
            continue
        ret[k] = get_data(IN[k])
    IN.close()

    return ret


def fix_string_representation(data, fields):
    """Due to a misinterpratation of data types, all strings in the matlab HDF5
        are represented as INT arrays. This function aims on fixing that."""

    for field in fields:
        if '/' in field:
            _field1, _field2 = field.split('/')
            shp = getattr(data, _field1).shape
            if len(shp) > 1:
                for i in xrange(shp[0]):
                    for j in xrange(shp[1]):
                        setattr(getattr(data, _field1)[i, j], _field2, ''.join([unichr(x) for x in getattr(getattr(data, _field1)[i, j], _field2)[0, :]]))
            else:
                for i in xrange(shp[0]):
                    setattr(getattr(data, _field1)[i], _field2, ''.join([unichr(x) for x in getattr(getattr(data, _field1)[i], _field2)[0, :]]))
        else:
            if isinstance(data, sp.ndarray):
                shp = data.shape
                if len(shp) > 1:
                    for i in xrange(shp[0]):
                        for j in xrange(shp[1]):
                            setattr(data[i, j], field, ''.join([unichr(x) for x in getattr(data[i, j], field)[0, :]]))
                else:
                    for i in xrange(shp[0]):
                        setattr(data[i], field, ''.join([unichr(x) for x in getattr(data[i], field)[0, :]]))
            else:
                shp = getattr(data, field).shape
                if len(shp) > 1:
                    for i in xrange(shp[0]):
                        for j in xrange(shp[1]):
                            getattr(data, field)[i, j] = ''.join([unichr(x) for x in getattr(data, field)[i, j][0, :]])
                else:
                    for i in xrange(shp[0]):
                        getattr(data, field)[i] = ''.join([unichr(x) for x in getattr(data, field)[i][0, :]])


def fix_array_structure(data):

    fields = filter(lambda x: x[:2] != '__', dir(data))
    shp = getattr(data, fields[0]).shape
    _data = sp.empty(shp, dtype='object')

    for f, field in enumerate(fields):
        if len(shp) > 1:
            for i in xrange(shp[0]):
                for j in xrange(shp[1]):
                    if f == 0:
                        _data[i, j] = hdf5data()
                    setattr(_data[i, j], field, getattr(data, field)[i, j])
        else:
            for i in xrange(shp[0]):
                if f == 0:
                    _data[i] = hdf5data()
                setattr(_data[i], field, getattr(data, field)[i])
    return _data

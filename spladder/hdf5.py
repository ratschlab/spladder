import sys
import scipy as sp

def appendToHDF5(file, data, name, faxis=0, daxis=0):
    """
    Goal of this function is to append more data to and
    existing HDF5 data entry.
    The dimensions other than the appending dimension have to match.

    
    """

    ### get shape of file object we append to
    fshape = file[name].shape
    ### get shape of data we want to append
    dshape = data.shape

    ### force one dimensional input data into the right form
    if len(fshape) > 1 and len(dshape) == 1 and faxis == 0:
        data = data[sp.newaxis, :]
        dshape = data.shape

    shapediff = len(fshape) - len(dshape)
    assert shapediff in [0, 1]

    ### check whether axes have been chosen correctly
    assert faxis < len(fshape)
    assert daxis < len(dshape)

    ### check that axes are compatible
    cfaxis = [x for i, x in enumerate(fshape) if i != faxis]
    if shapediff == 0:
        cdaxis = [x for i, x in enumerate(dshape) if i != daxis]
    else:
        cdaxis = daxis
    assert(sp.all(cfaxis == cdaxis))

    ### compute new shape, resize and add
    newshape = [x if i != faxis else x + data.shape[daxis] for i, x in enumerate(fshape)]
    file[name].resize(newshape)
    append_string = 'file[name][' + ''.join([':,' if i != faxis else 'fshape[%i]:,' % i for i, x in enumerate(fshape)]).strip(',') + '] = data'
    exec(append_string)



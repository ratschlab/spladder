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
    if fshape > 1 and dshape == 1 and faxis == 0:
        data = data[sp.newaxis, :]
        dshape = data.shape

    shapediff = len(fshape) - len(dshape)
    assert shapediff in [0, 1]

    ### check whether axes have been chosen correctly
    assert faxis < fshape
    assert daxis < dshape

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


#def appendToHDF5_(file, data, name, axis=0):
#
#    ### get current shape 
#    tmp = file[name].shape
#    ### resize
#    if len(tmp) == 1:
#        file[name].resize((tmp[0] + data.shape[0],))
#        file[name][tmp[0]:] = data
#    elif len(tmp) == 2:
#        if axis == 0:
#            file[name].resize((tmp[0] + data.shape[0], tmp[1]))
#            file[name][tmp[0]:, :] = data
#        else:
#            if len(data.shape) < 2:
#                file[name].resize((tmp[0], tmp[1] + 1))
#                file[name][:, tmp[1]:] = data[:, sp.newaxis]
#            else:
#                file[name].resize((tmp[0], tmp[1] + data.shape[1]))
#                file[name][:, tmp[1]:] = data
#    elif len(tmp) == 3:
#        if axis == 0:
#            file[name].resize((tmp[0] + data.shape[0], tmp[1], tmp[2]))
#            file[name][tmp[0]:, :, :] = data
#        elif axis == 1:
#            if len(data.shape) < 3:
#                file[name].resize((tmp[0], tmp[1] + 1, tmp[2]))
#                file[name][:, tmp[1]:, :] = data[:, sp.newaxis, :]
#            else:
#                file[name].resize((tmp[0], tmp[1] + data.shape[1], tmp[2]))
#                file[name][:, tmp[1]:, :] = data
#        else:
#            if len(data.shape) < 3:
#                file[name].resize((tmp[0], tmp[1], tmp[2] + 1))
#                file[name][:, :, tmp[2]:] = data[:, :, sp.newaxis]
#            else:
#                file[name].resize((tmp[0], tmp[1], tmp[2] + data.shape[1]))
#                file[name][:, :, tmp[2]:] = data
#
#    else:
#        print >> sys.stderr, "cannot append data to HDF5 with more than 3 dimensions" 
#        sys.exit(-1)



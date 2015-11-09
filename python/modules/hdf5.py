
def appendToHDF5(file, data, name, axis=0):
    
    ### get current shape 
    tmp = file[name].shape
    ### resize
    if len(tmp) == 1:
        file[name].resize((tmp[0] + data.shape[0],))
        file[name][tmp[0]:] = data
    elif len(tmp) == 2:
        if axis == 0:
            file[name].resize((tmp[0] + data.shape[0], tmp[1]))
            file[name][tmp[0]:, :] = data
        else:
            file[name].resize((tmp[0], tmp[1] + 1))
            if len(data.shape) < 2:
                file[name][:, tmp[1]:] = data[:, sp.newaxis]
            else:
                file[name][:, tmp[1]:] = data
    else:
        print >> sys.stderr, "cannot append data to HDF5 with more than 2 dimensions" 
        sys.exit(-1)



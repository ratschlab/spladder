import scipy as sp

def unique_rows(array, index = None):
    """Make array unique by rows"""

    if array.shape[0] == 0:
        return array
    if index == True:
        array_u, index = sp.unique(array.view([('', array.dtype)] * array.shape[1]), return_index=True)
        return (array_u.view(array.dtype).reshape(array_u.shape[0], array.shape[1]), index)
    else:
        array_u = sp.unique(array.view([('', array.dtype)] * array.shape[1]))
        return array_u.view(array.dtype).reshape(array_u.shape[0], array.shape[1])

def sort_rows(array, index = None):
    """Sort array by rows"""

    if array.shape[0] == 0:
        return array

    array_v = array.view([('', array.dtype)] * array.shape[1])
    if index == True:
        return (sp.sort(array_v, axis = 0).view(array.dtype).reshape(array_u.shape[0], array.shape[1]), sp.argsort(array_v, axis = 0).view(array.dtype)[:, 0])
    else:
        return sp.sort(array_v, axis = 0).view(array.dtype).reshape(array_u.shape[0], array.shape[1])

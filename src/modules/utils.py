import scipy as sp

def isequal(A, B):
    
    if not sp.all(A.shape == B.shape):
        return False
    else:
        return sp.all(A == B)

def issubset(A, B):
    # function found = issubset(A,B)
    #
    # returns true if A is a subset of B, where 
    # both A and B are vectors with 1 where elements exists and 0 otherwise.
    return sp.all(sp.where(A)[0] == sp.in1d(sp.where(A)[0], sp.where(B)[0]))


def intersect_rows(array1, array2, index = None):
    """Return intersection of rows"""

    if (array1.shape[0] == 0) or (array2.shape[0] == 0):
        if index == True:
            return (array, sp.zeros((0,)), sp.zeros((0,)))
        else:
            return array

    array1_v = array1.view([('', array1.dtype)] * array1.shape[1])
    array2_v = array2.view([('', array2.dtype)] * array2.shape[1])
    array_i = sp.intersect1d(array1_v, array2_v)

    if index == True:
        a1_i = sp.where(sp.in1d(array1_v, array_i))[0]
        a2_i = sp.where(sp.in1d(array2_v, array_i))[0]
        return (array_i.view(array1.dtype).reshape(array_i.shape[0], array1.shape[1]), a1_i, a2_i)
    else:
        return array_i.view(array1.dtype).reshape(array_i.shape[0], array1.shape[1])


def unique_rows(array, index = None):
    """Make array unique by rows"""

    if array.shape[0] == 0:
        if index == True:
            return (array, [])
        else:
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
        return (sp.sort(array_v, axis = 0).view(array.dtype).reshape(array_v.shape[0], array.shape[1]), sp.argsort(array_v, axis = 0).view(array.dtype)[:, 0])
    else:
        return sp.sort(array_v, axis = 0).view(array.dtype).reshape(array_v.shape[0], array.shape[1])

def ismember(element, array, rows=False):
    """Check if element is member of array"""

    if rows:
        return sp.any([sp.all(array[x, :] == element) for x in range(array.shape[0])])
    else:
        return sp.all([element[i] in array for i in element.shape[0]])

        ### TODO I think this is not quite right yet


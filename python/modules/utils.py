import scipy as sp

def isequal(A, B):
    
    if A is None or B is None:
        return False

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


#def intersect_rows(array1, array2, index = None):
#    """Return intersection of rows"""
#
#    if (array1.shape[0] == 0):
#        if index == True:
#            return (array1, sp.zeros((0,)), sp.zeros((0,)))
#        else:
#            return array1
#    if (array2.shape[0] == 0):
#        if index == True:
#            return (array2, sp.zeros((0,)), sp.zeros((0,)))
#        else:
#            return array2
#
#    array1_v = array1.view([('', array1.dtype)] * array1.shape[1])
#    array2_v = array2.view([('', array2.dtype)] * array2.shape[1])
#    array_i = sp.intersect1d(array1_v, array2_v)
#
#    if index == True:
#        a1_i = sp.where(sp.in1d(array1_v, array_i))[0]
#        a2_i = sp.where(sp.in1d(array2_v, array_i))[0]
#        return (array_i.view(array1.dtype).reshape(array_i.shape[0], array1.shape[1]), a1_i, a2_i)
#    else:
#        return array_i.view(array1.dtype).reshape(array_i.shape[0], array1.shape[1])


def unique_rows(array, index = None):
    """Make array unique by rows"""

    if array.shape[0] == 0:
        if index == True:
            return (array, [])
        else:
            return (array)

    if len(array.shape) == 1:
        if index == True:
            return (array, [0])
        else:
            return array

    (array_s, s_idx) = sort_rows(array, True)
    tmp = [False]
    tmp.extend([sp.all(array_s[i-1, :] == array_s[i, :]) for i in range(1, array.shape[0])])
    k_idx = sp.where(~sp.array(tmp, dtype='bool'))[0]

    if index == True:
        return (array[s_idx[k_idx], :], s_idx[k_idx])
    else:
        return array[s_idx[k_idx], :]


def intersect_rows(array1, array2, index=False):
    """Return array with rows that intersect between array1 and array2"""

    tmp1 = sp.array(['-'.join(array1[i, :].astype('str')) for i in range(array1.shape[0])])
    tmp2 = sp.array(['-'.join(array2[i, :].astype('str')) for i in range(array2.shape[0])])
    
    idx = sp.where(sp.in1d(tmp1, tmp2))[0]
    if index:
        idx2 = sp.where(sp.in1d(tmp2, tmp1))[0]

    if index:
        return (array1[idx, :], idx, idx2)
    else:
        return (array1[idx, :], None, None)

def sort_rows(array, index = None):
    """Sort array by rows"""

    ### empty array
    if array.shape[0] == 0:
        if index == True:
            return (array, [])
        else:
            return (array)

    ### only one row
    if len(array.shape) == 1:
        if index == True:
            return (array, [0])
        else:
            return (array)

    ### more than one row
    s_idx = sp.lexsort([array[:, -i] for i in range(1, array.shape[1] + 1)])

    if index == True:
        return (array[s_idx, :], s_idx)
    else:
        return array[s_idx, :]


def ismember(element, array, rows=False):
    """Check if element is member of array"""

    if rows:
        return sp.any([sp.all(array[x, :] == element) for x in range(array.shape[0])])
    else:
        return sp.all([element[i] in array for i in element.shape[0]])

        ### TODO I think this is not quite right yet

def replace_sub_matrix(mat_in, idx, mat_put):
    """Replaces the values in mat_in in rows and cols idx with values of mat_put"""
    
    assert((idx.shape[0] * idx.shape[0]) == mat_put.ravel().shape[0])

    sp.put(mat_in, sp.ravel_multi_index([[x for x in idx for _ in idx], [x for _ in idx for x in idx]], (mat_in.shape[0], mat_in.shape[1])), mat_put.ravel())

    return mat_in



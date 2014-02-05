def unique_rows(array):
    """Make array unique by rows"""

    if array.shape[0] == 0:
        return array

    array_u = np.unique(array.view([('', array.dtype)] * array.shape[1]))
    return array_u.view(array.dtype).reshape(array_u.shape[0], array.shape[1])

# TODO: Eventually look into ctypes or numba for streamlining the memory
#       of this class
# NOTE: This implementation of deuterater does not actually need the 'i'
#       attribute of this class. Consult the program logic in the extract
#       function in extract.py to see why.


class Peak(object):
    '''Contains the data for one peak

    Attributes
    ----------
    mz : float
        The mass-to-charge ratio of the peak
    ab : float
        The abundance value of the peak (also referred to as the intensity)
    i : int
        The peak number (i.e. m0, m1, m2)
    '''
    __slots__ = (
        # defining '__slots__' lets the python interpreter know what fields
        # we will be defining in the class
        'mz',
        'ab',
        'i',
    )

    def __init__(self, mz, abundance, i=None):
        self.mz = mz
        self.ab = abundance
        self.i = i

    # Defining the __repr__ function allows python to call repr()
    # on this object. This is usually much less formatted than the related
    # '__str__' function
    def __repr__(self):
        return 'Peak(mz={0},ab={1},i={2})'.format(
            self.mz, self.ab, self.i
        )

    # Defining the __str__ function allows python to call str()
    # on this object. This is usually the best way to define a 'toString'
    # or similar function
    def __str__(self):
        return 'Peak(mz={:10.2f}, ab={:10.2f}, i={:3d})'.format(
            self.mz, self.ab, self.i
        )

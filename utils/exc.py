'''
Module containing all of Deuterater's custom errors and warnings

'''
import warnings  # noqa: F401


class InvalidHeaderError(RuntimeError):
    '''Exception thrown when headers do not match Deuterater's specifications

    Parameters
    ----------
    msg : str
        Human readable string describing the exception.
    original_e : :obj:`Exception`
        original exception thrown, in order to carry through

    '''
    def __init__(self, msg, original_e=None):
        super(InvalidHeaderError, self).__init__(msg + (': %s' % original_e))
        self.original_exception = original_e


class PeakIndexError(IndexError):
    '''Exception thrown when headers do not match Deuterater's specifications

    Parameters
    ----------
    msg : str
        Human readable string describing the exception.
    original_e : :obj:`Exception`
        original exception thrown, in order to carry through

    '''
    def __init__(self, msg, original_e=None):
        super(PeakIndexError, self).__init__(msg + (': %s' % original_e))
        self.original_exception = original_e


class EmptyIdChunkWarning(RuntimeWarning):
    '''Warning raised when a chunk has no ids

    Parameters
    ----------
    msg : str
        Human readable string describing the exception.
    original_e : :obj:`Exception`
        original exception thrown, in order to carry through

    '''
    def __init__(self, msg, original_e=None):
        super(EmptyIdChunkWarning, self).__init__(msg + (': %s' % original_e))
        self.original_exception = original_e


# class SingleElementIdChunkWarning(RuntimeWarning):
#     '''Exception thrown when a chunk has only one id

#     Parameters
#     ----------
#     msg : str
#         Human readable string describing the exception.
#     original_e : :obj:`Exception`
#         original exception thrown, in order to carry through

#     '''
#     def __init__(self, msg, original_e=None):
#         super(SingleElementIdChunkWarning, self).__init__(
#             msg + (': %s' % original_e)
#         )
#         self.original_exception = original_e


class OutOfMZMLTimeBoundsWarning(RuntimeWarning):
    '''Warning raised when retention time is not in the span of the mzml

    Parameters
    ----------
    msg : str
        Human readable string describing the exception.
    original_e : :obj:`Exception`
        original exception thrown, in order to carry through

    '''
    def __init__(self, msg, original_e=None):
        super(OutOfMZMLTimeBoundsWarning, self).__init__(
            msg + (': %s' % original_e)
        )
        self.original_exception = original_e


class InvalidSettingsError(RuntimeError):
    '''Exception thrown when a setting has been entered incorrectly

    Parameters
    ----------
    msg : str
        Human readable string describing the exception.
    original_e : :obj:`Exception`
        original exception thrown, in order to carry through

    '''
    def __init__(self, msg, original_e=None):
        super(InvalidSettingsError, self).__init__(
            msg + (': %s' % original_e)
        )
        self.original_exception = original_e


class InvalidSettingsWarning(RuntimeWarning):
    '''Exception thrown when a setting has been entered incorrectly

    Parameters
    ----------
    msg : str
        Human readable string describing the exception.
    original_e : :obj:`Exception`
        original exception thrown, in order to carry through

    '''
    def __init__(self, msg, original_e=None):
        super(InvalidSettingsWarning, self).__init__(
            msg + (': %s' % original_e)
        )
        self.original_exception = original_e

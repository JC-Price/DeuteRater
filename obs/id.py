from obs.envelope import Envelope  # noqa: 401
from obs.peak import Peak
from utils.math import angle_between  # noqa: 401
import deuterater.settings as settings

from numpy import mean, median, std, argsort

class ID(object):
    '''Contains isotoptic envelope data

    TODO: we need to discuss how technical we want this definition to be

    Attributes
    ----------
    _envelopes : :obj:`list` of :obj:`Envelope`
        The isotopic envelopes corresponding to this identification
    _settings : str
        the path of the settings file
    condensed_envelope : :obj:`Envelope`
        The envelope representing the aggregated data of the envelopes
    rt : float
        The retention time of the identification
    rt_min : float
        The lowest retention time of the identification
    rt_max : float
        The highest retention time of the identification
    mz : float
        The mass-to-charge ratio given by the identification file
    mass : float
        The neutral mass of the unlabeled identification
    z : int
        Charge of the identification
    n_isos : int
        The number of isotopes (or isotopic peaks) to look for
    max_m0_abundance : float
        The maximum m0 abundance of the identification's isotopic envelopes
    mads : :obj:`list` of float
        The median absolute deviations of abundances at each peak location
    '''
    # defining '__slots__' lets the python interpreter know what fields
    # we will be defining in the class
    __slots__ = (
        '_envelopes',
        '_settings',
        'condensed_envelope',
        'rt',
        'rt_min',
        'rt_max',
        'mz',
        'mass',
        'z',
        'n_isos',
        #'cf',
        # 'envelopes_before_angle_filter',
        # 'envelopes_after_angle_filter',
        # 'deviation_before_angle_filter',
        # 'deviation_after_angle_filter',
        # 'deviation_outside_angle_filter',
        # TODO: Do we need the max m0 abundance for anything else?
        # 'max_m0_abundance',
        'mads'
    )

    def __init__(self, rt, mz, mass, z, n_isos):#, cf):
        # NOTE: passing in an empty list as default constructor makes
        #       all ids in the chunk share the list internally this is a result
        #       of quirks of how python handles objects and default
        #       constructors
        self._envelopes = []
        self.condensed_envelope = None

        self.rt = rt
        self.mz = mz
        self.mass = mass
        self.z = z
        self.n_isos = n_isos
        self.mads = []
        #self.cf = cf

    # Defining the __repr__ function allows python to call repr()
    # on this object. This is usually much less formatted than the related
    # '__str__' function
    def __repr__(self):
        # TODO: clean this up and add metainformation
        return 'ID({})'.format(
            ','.join([repr(envelope) for envelope in self._envelopes])
        )

    # Defining the __str__ function allows python to call str()
    # on this object. This is usually the best way to define a 'toString'
    # or similar function
    def __str__(self):
        # TODO: clean this up and add metainformation

        return 'ID(rt={:.3f}, mz={}, mass={:.3f}, z={}, n_isos={})'.format(
            self.rt, self.mz, self.mass, self.z, self.n_isos
        )

    def append_envelope(self, envelope):
        # TODO: add docstring
        if not isinstance(envelope, Envelope):
            raise TypeError('"{}" is not of type Envelope'.format(envelope))
        self._envelopes.append(envelope)

    def aggregate_envelopes(self):
        # TODO: add docstring

        # First, we see if there are even enough envelopes to warrant analyzing
        #   this identification
        if len(self._envelopes) < settings.min_envelopes_to_combine:
            # TODO: not enough data. What should be logged?
            return

        def get_max_M0_vector(self):
            vector_list = [[peak.ab for peak in envelope.get_peaks()]
                           for envelope
                           in self._envelopes]

            m0_list = [v[0] for v in vector_list]
            max_m0_ab, i_max = max([(v, i) for i, v in enumerate(m0_list)])

            max_vector = vector_list[i_max]

            return max_m0_ab

        def cut_finger_filter(self):
            valid_indices = list()
            invalid_indices = list()
            def local_min(list_of_floats):
                for i in range(1, len(list_of_floats) - 1):
                    if list_of_floats[i - 1] > list_of_floats[i] and list_of_floats[i + 1] > list_of_floats[i]:
                        return False
                return True

            vector_list = [[peak.ab for peak in envelope.get_peaks()]
                           for envelope
                           in self._envelopes]

            cut_finger_mask = [local_min(vector) for vector in vector_list]
            filtered_vector_list = [vector_list[i] for i in range(len(vector_list)) if cut_finger_mask[i]]
            for i in range(len(cut_finger_mask)):
                if cut_finger_mask[i]:
                    valid_indices.append(i)
                else:
                    invalid_indices.append(i)

            self._envelopes = [self._envelopes[i] for i in valid_indices]

            return self._envelopes

        # relative angle filter (What is a better way to state this?)
        def DR_3_5_angle_filter(self, max_valid_angle=1.2):
            # minLen is removed because the length of each envelope should all be the exact same
            filtered_envelopes = []
            angle_list = []
            minAngle = 1000000.0  # we want to subtract the smallest angle from all the rest in a sort of normalization

            # Access the abundances in each of the envelopes in this identification
            #   and orginize as a list of lists. (Uses nested list comprehensions)
            vector_list = [[peak.ab for peak in envelope.get_peaks()]
                           for envelope
                           in self._envelopes]

            # this loop finds the smallest angle and puts all the angles in an array
            for i in range(1, len(vector_list)):
                prev = vector_list[i - 1]
                cur = vector_list[i]
                angle = angle_between(prev, cur)
                angle_list.append(angle)

                if angle < minAngle:
                    minAngle = angle

            # now we subtract that smallest angle and antilog transform to increase spread between good and bad data
            for i in range(len(angle_list)):
                angle_list[i] = 10.0 ** (angle_list[i] - minAngle)

            # Note that the construction of this for loop means that the 1st item is
            #   compared to the 2nd to determine how good the 1st item is
            # this does result in automatically discarding the last set of data in envelopes
            #   (granted it was probably bad anyway, the ends normally are)
            # an angle of 1-1.2 is generally good and we grab only those isotopes that correspond as such

            valid_indices = []
            invalid_indices = []
            for i in range(0, len(angle_list)):
                if angle_list[i] < max_valid_angle:
                    valid_indices.append(i)
                else:
                    invalid_indices.append(i)

            self._envelopes = [self._envelopes[i] for i in valid_indices]

            return self._envelopes

            # return angle_list, invalid_indices, valid_indices

        # reference vector angle filter
        def DR_4_0_angle_filter(self, max_valid_angle=3.14159, reference_vector=None):
            # Angles are in Radians
            # Determine the reference vector

            # Access the abundances in each of the envelopes in this identification
            #   and organize as a list of lists. (Uses nested list comprehensions)
            vector_list = [[peak.ab for peak in envelope.get_peaks()]
                           for envelope
                           in self._envelopes]

            def reference_vector_max_m0(vector_list):
                m0_list = [v[0] for v in vector_list]
                max_m0_ab, i_max = max([(v, i) for i, v in enumerate(m0_list)])

                max_vector = vector_list[i_max]

                return max_vector

            if reference_vector is None:
                reference_vector = reference_vector_max_m0(vector_list)
                print(reference_vector)

            # TODO: check if its faster to use generator expression
            angle_list = [angle_between(reference_vector, vector)
                          for vector in vector_list]

            # self.envelopes_before_angle_filter = len(self._envelopes)
            # self.deviation_before_angle_filter = std(angle_list)

            # Throw away bad envelopes
            valid_indices = []
            invalid_indices = []
            for i, angle in enumerate(angle_list):
                if angle < max_valid_angle:
                    valid_indices.append(i)
                else:
                    invalid_indices.append(i)

            return angle_list, invalid_indices, valid_indices

        def highest_intensity_scans_filter(self, MAX_SCANS=60, MAX_SCANS_TYPE="M0"):
            if len(self._envelopes) >= MAX_SCANS:
                if MAX_SCANS_TYPE == "angle":
                    sorted_angles = argsort(angle_list)
                    self._envelopes = [self._envelopes[sorted_angles[i]] for i in range(MAX_SCANS)]
                elif MAX_SCANS_TYPE == "M0":
                    vector_list = [[peak.ab for peak in envelope.get_peaks()]
                                   for envelope
                                   in self._envelopes]
                    m0_list = [v[0] for v in vector_list]
                    sorted_M0_scans = argsort(m0_list)
                    self._envelopes = [self._envelopes[sorted_M0_scans[-(i + 1)]] for i in range(MAX_SCANS)]

            return self._envelopes

        # from deuterater.emass import emass
        # emass_output = emass(self.cf, 0, 1, 0, 0, self.n_isos)
        # reference_vector = emass_output[1].iloc[0].to_list()[1:]
        # angle_list, invalid_indices, valid_indices = DR_4_0_angle_filter(self, settings.max_valid_angle, reference_vector)

        self._envelopes = cut_finger_filter(self)

        if len(self._envelopes) < settings.min_envelopes_to_combine:
            return

        self._envelopes = DR_3_5_angle_filter(self)

        if len(self._envelopes) < settings.min_envelopes_to_combine:
            return

        self._envelopes = highest_intensity_scans_filter(self)

        if len(self._envelopes) < settings.min_envelopes_to_combine:
            return

        # NOTE: this threshold was defined in the original deuterater
        threshold = 1 / settings.peak_ratio_denominator
        self._envelopes = [
            envelope for envelope in self._envelopes
            if envelope.get_m(0).ab > threshold * get_max_M0_vector(self)
               # and envelope.get_m(0).ab > envelope.baseline
        ]

        if len(self._envelopes) < settings.min_envelopes_to_combine:
            return


        # NOTE: we would perform any normaliztions here. We decided not to.

        # TODO: make this happen in only one pass
        self.rt_min = min(e.rt for e in self._envelopes)
        self.rt_max = max(e.rt for e in self._envelopes)
        # self.max_m0_abundance = max_m0_ab

        # After performing all of this filtration, aggregate all of the data
        #   in the remaining envelopes in to a 'condensed envelope'
        self.condense_envelopes()

    def condense_envelopes(self):
        # TODO: add docstring
        def condense_peak(peak_num):
            mzs = [envelope.get_m(peak_num).mz
                   for envelope in self._envelopes]
            abundances = [envelope.get_m(peak_num).ab
                          for envelope in self._envelopes]
            mz_median = median(mzs)
            ab_median = median(abundances)
            absolute_mz_deviations = [abs(mz - mz_median) for mz in mzs]
            mad = median(absolute_mz_deviations)
            if mad == 0:
                peak = Peak(mz_median, ab_median, peak_num)
            else:
                # TODO: What is this magic number from?
                z_scores = [0.6745 * dev / mad
                            for dev in absolute_mz_deviations]
                filtered_mzs = [mzs[i]
                                for i, z_score
                                in enumerate(z_scores)
                                if z_score < settings.zscore_cutoff]
                # TODO: is it ok to still take the intensity medians?
                peak = Peak(mean(filtered_mzs), ab_median, peak_num)
            return peak, mad

        # TODO: sum the data together instead of a median
        # TODO: add a number of scans that passed

        # mads as optional column

        peak_list = []
        for peak_num in range(0 - settings.peak_lookback, 0):
            # TODO should we keep all the mads instead?
            peak, _ = condense_peak(peak_num)
            peak_list.append(peak)
        for peak_num in range(0, self.n_isos):
            peak, mad = condense_peak(peak_num)
            peak_list.append(peak)
            self.mads.append(mad)
        for peak_num in range(self.n_isos, self.n_isos + settings.peak_lookahead):
            # TODO should we keep all the mads instead?
            peak, _ = condense_peak(peak_num)
            peak_list.append(peak)
        rt_median = median([envelope.rt for envelope in self._envelopes])
        self.condensed_envelope = Envelope(
            peaks=peak_list,
            rt=rt_median,
            n_lookback=settings.peak_lookback,
            n_lookahead=settings.peak_lookahead
        )

        # TODO: mean or median?
        # NOTE: if we sum the signal then we need to sum the baseline
        self.condensed_envelope.baseline = \
            median([envelope.baseline for envelope in self._envelopes])

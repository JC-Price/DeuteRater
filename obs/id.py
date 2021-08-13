try:
    from obs.envelope import Envelope  # noqa: 401
    from obs.peak import Peak
    from utils.math import angle_between  # noqa: 401
    import deuterater.settings as settings
except:
    from DeuteRater.obs.envelope import Envelope  # noqa: 401
    from DeuteRater.obs.peak import Peak
    from DeuteRater.utils.math import angle_between  # noqa: 401
    import DeuteRater.deuterater.settings as settings

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
        # TODO: Do we need the max m0 abundance for anything else?
        # 'max_m0_abundance',
        'mads',
        '_unfiltered_envelopes',
        "rt_windows",
        "rt_peak_index",
        "neutromer_peak_maximums",
        "is_valid",
        "signal_noise"
    )

    def __init__(self, rt, mz, mass, z, n_isos):#, cf):
        # NOTE: passing in an empty list as default constructor makes
        #       all ids in the chunk share the list internally this is a result
        #       of quirks of how python handles objects and default
        #       constructors
        self._envelopes = []
        self._unfiltered_envelopes = None
        self.condensed_envelope = None

        self.rt = rt
        self.mz = mz
        self.mass = mass
        self.z = z
        self.n_isos = n_isos
        self.mads = []
        self.is_valid = False
        self.rt_windows = []
        self.neutromer_peak_maximums = []
        self.rt_peak_index = []
        self.signal_noise = []
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

        from copy import deepcopy
        
        valid_envelopes = [envelope for envelope in self._envelopes if envelope.is_valid]

        # First, we see if there are even enough envelopes to warrant analyzing
        #   this identification
        if len(valid_envelopes) < settings.min_envelopes_to_combine:
            # TODO: not enough data. What should be logged?
            self._unfiltered_envelopes = None
            return
        
        # [e.baseline for e in self._envelopes]  # Calculate the baseline for any used scans.
        
        if settings.use_chromatography_division != "No":
            self._unfiltered_envelopes = deepcopy(self._envelopes)
        
        self._envelopes = valid_envelopes
        
        def get_max_M0_vector(current_id):
            vector_list = [[peak.ab for peak in envelope.get_peaks()]
                           for envelope
                           in current_id._envelopes]

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
        def DR_3_5_angle_filter(current_id, max_valid_angle=1.2):
            # minLen is removed because the length of each envelope should all be the exact same
            filtered_envelopes = []
            angle_list = []
            minAngle = 1000000.0  # we want to subtract the smallest angle from all the rest in a sort of normalization

            # Access the abundances in each of the envelopes in this identification
            #   and organize as a list of lists. (Uses nested list comprehensions)
            vector_list = [[peak.ab for peak in envelope.get_peaks()]
                           for envelope
                           in current_id._envelopes]

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

            current_id._envelopes = [current_id._envelopes[i] for i in valid_indices]

            return current_id._envelopes

            # return angle_list, invalid_indices, valid_indices

        def highest_intensity_scans_filter(self, MAX_SCANS=60):
            if len(self._envelopes) >= MAX_SCANS:
                vector_list = [[peak.ab for peak in envelope.get_peaks()]
                               for envelope
                               in self._envelopes]
                m0_list = [v[0] for v in vector_list]
                sorted_M0_scans = argsort(m0_list)
                self._envelopes = [self._envelopes[sorted_M0_scans[-(i + 1)]] for i in range(MAX_SCANS)]

            return self._envelopes

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

        rt_list = self._get_rt_list()
        self.rt_min = min(rt_list)
        self.rt_max = max(rt_list)
        
        # Calculate Signal Noise:
        from scipy.signal import savgol_filter

        smoothing_width = settings.smoothing_width
        smoothing_order = settings.smoothing_order

        from numpy import array, sqrt, std

        def mad(values):
            m = median(values)
            return median([abs(a - m) for a in values])

        def rmse(actual, estimated):
            s = sum((actual - estimated) ** 2)
            return sqrt(s / len(actual))

        temp = smoothing_width
        for a in self._get_peak_list():
            if len(a) < smoothing_width:
                smoothing_width = (int((len(a) - 1)/2) * 2) + 1
            smoothed_curves = savgol_filter(a, smoothing_width, smoothing_order)
            smoothing_width = temp
            smoothed_curves[smoothed_curves < 0] = 0
            diff = [b for b in array(a) - smoothed_curves]
            # abs_diff = [abs(b) for b in diff]
            normal_distribution_scale_factor = 1.4826
            self.signal_noise.append(mad(diff) * normal_distribution_scale_factor * 3)
            # self.signal_noise.append(rmse(array(a), smoothed_curves) * 3)
            # self.signal_noise.append(mean(abs_diff) * 3)
            # self.signal_noise.append(median(abs_diff) * 3)

            pass
        
        # After performing all of this filtration, aggregate all of the data
        #   in the remaining envelopes in to a 'condensed envelope'
        self.condense_envelopes()

    def _get_rt_list(self, rounded=False):
        if rounded:
            return [round(float(envelope.rt), 4) for envelope in self._envelopes]
        return [float(envelope.rt) for envelope in self._envelopes]

    def _get_peak_list(self, intensity_filter=0):
        peaks_list = list()
    
        # Get a 2D list of the intensities
        for i in range(self.n_isos):
            peaks_list.append(list())
    
        for envelope in self._envelopes:
            for i in range(self.n_isos):
                if envelope.get_m(i).ab < intensity_filter:
                    peaks_list[i].append(0)
                else:
                    peaks_list[i].append(envelope.get_m(i).ab)
        
        return peaks_list
    
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
        rt_median = median([float(envelope.rt) for envelope in self._envelopes])
        self.condensed_envelope = Envelope(
            peaks=peak_list,
            rt=rt_median,
            n_lookback=settings.peak_lookback,
            n_lookahead=settings.peak_lookahead
        )

        # TODO: mean or median?
        # NOTE: if we sum the signal then we need to sum the baseline
        self.condensed_envelope.baseline = \
            median([float(envelope.baseline) for envelope in self._envelopes])

    def divide_chromatography(self, should_plot=False):
        import scipy.ndimage as ndi
        import scipy.signal as sig
        import numpy as np

        sampling_rate = settings.sampling_rate
        smoothing_width = settings.smoothing_width
        smoothing_order = settings.smoothing_order
        intensity_filter = settings.intensity_filter

        peaks_list = self._get_peak_list(intensity_filter)
        combined_peaks = np.sum(np.array(peaks_list), axis=0)

        # REMOVE THESE LINES LATER, ONCE WE ACTUALLY NEED TO USE THIS!
        # print("REMOVE LINE 379-383 (or around there)")
        if max(combined_peaks) < intensity_filter:
            return
        combined_data = [combined_peaks[i] if combined_peaks[i] > intensity_filter else 0 for i in range(combined_peaks.size)]

        check_valid_list = [self._envelopes[i].is_valid and combined_data[i] > intensity_filter for i in
                            range(len(combined_data))]
        if True not in check_valid_list:
            return

        self.is_valid = True

        if settings.how_divided == "combined":
            division_data = combined_data
        elif settings.how_divided == "m0":
            division_data = peaks_list[0]
        elif settings.how_divided == "m1":
            division_data = peaks_list[1]
        elif settings.how_divided == "m2":
            division_data = peaks_list[2]
        else:
            division_data = combined_data

        # Gaussian smoothing:
        from scipy.signal import savgol_filter
        smoothed_curves = savgol_filter(division_data, smoothing_width, smoothing_order)
        smoothed_curves[smoothed_curves < 0] = 0
        gaussian_curves = ndi.gaussian_filter1d(smoothed_curves, sampling_rate)
        gaussian_min_cuts = sig.find_peaks([-x for x in gaussian_curves])[0]
        peak_finding = np.concatenate(
            (np.array([gaussian_curves[1]]), gaussian_curves, np.array([gaussian_curves[-2]])))

        max_cuts = sig.argrelmax(peak_finding)[0]
        max_cuts = [item - 1 for item in max_cuts]
        gaussian_max_cuts = max_cuts

        left_points = list()
        right_points = list()
        widths = list()
        heights = list()
        peaks = gaussian_max_cuts
        min_points = gaussian_min_cuts
        min_points = np.append(min_points, len(gaussian_curves))
        for i in range(len(peaks)):
            if i == 0:
                peak_position = peaks[i]
                # Left Sided Peak:
                left_peak = np.concatenate((gaussian_curves[:peaks[i] + 1],
                                            np.flip(gaussian_curves[:peaks[i]])))
                right_peak = np.concatenate((np.flip(gaussian_curves[peaks[i] + 1:min_points[0]]),
                                             gaussian_curves[peaks[i]:min_points[0]]))

                if int(len(left_peak)/2) != 0:
                    left_sides = sig.peak_widths(left_peak, [int(len(left_peak) / 2)], rel_height=settings.rel_height)
                else:
                    left_sides = (np.array([0]), left_peak, np.array([0]), np.array([0]))
                if int(len(right_peak) / 2) != 0:
                    right_sides = sig.peak_widths(right_peak, [int(len(right_peak) / 2)], rel_height=settings.rel_height)
                else:
                    right_sides = (np.array([0]), right_peak, np.array([0]), np.array([0]))
                left_points.append(left_sides[2][0])
                right_points.append(peak_position - int(len(right_peak) / 2) + right_sides[3][0])
                widths.append(right_sides[0][0] / 2 + left_sides[0][0] / 2)
                heights.append(max([right_sides[1][0], left_sides[1][0]]))
            else:
                point_correction = min_points[i - 1]
                left_peak = np.concatenate((gaussian_curves[min_points[i - 1]:peaks[i] + 1],
                                            np.flip(gaussian_curves[min_points[i - 1]:peaks[i]])))
                right_peak = np.concatenate((np.flip(gaussian_curves[peaks[i] + 1:min_points[i]]),
                                             gaussian_curves[peaks[i]:min_points[i]]))

                if int(len(left_peak) / 2) != 0:
                    left_sides = sig.peak_widths(left_peak, [int(len(left_peak) / 2)], rel_height=settings.rel_height)
                else:
                    left_sides = (np.array([0]), left_peak, np.array([0]), np.array([0]))
                if int(len(right_peak) / 2) != 0:
                    right_sides = sig.peak_widths(right_peak, [int(len(right_peak) / 2)], rel_height=settings.rel_height)
                else:
                    right_sides = (np.array([0]), right_peak, np.array([0]), np.array([0]))
                left_points.append(point_correction + left_sides[2][0])
                right_points.append(peaks[i] - int(len(right_peak) / 2) + right_sides[3][0])
                widths.append(right_sides[0][0] / 2 + left_sides[0][0] / 2)
                heights.append(max([right_sides[1][0], left_sides[1][0]]))
        full_peak_width_list = (np.array(widths), np.array(heights), np.array(left_points).astype("int32"),
                                         np.array(right_points).astype("int32"))

        if should_plot:
            import matplotlib.pyplot as plt
    
            plot_colors = ['tab:blue', 'tab:orange', 'tab:purple', 'olivedrab', 'chocolate']
            gaussian_colors = ['b', 'orangered', 'indigo', 'lawngreen', 'saddlebrown']
    
            plt.plot(combined_data, label="Combined", color=plot_colors[0])
            plt.plot(gaussian_curves, label='Combined Gaussian', zorder=2, color=gaussian_colors[0])
            plt.vlines(x=full_peak_width_list[2], ymin=0, ymax=max(peak_finding), label="Start of Window", color="r")
            plt.vlines(x=full_peak_width_list[3], ymin=0, ymax=max(peak_finding), label="End of Window", color="g")
            plt.vlines(x=peaks, ymin=0, ymax=max(peak_finding), label="Max Intensity in Window", color="lightgrey")
    
            plt.title(str(self.mz))
            plt.legend()
            # plt.xticks(range(0, len(self._get_rt_list()), 25), [str(i) for i in self._get_rt_list(True)[0::25]])
            plt.show()
        
        self._neutromer_peak_variance(should_plot)

        # rt_values = [envelope.rt for envelope in self._envelopes]
        # left_values = [rt_values[a] for a in full_peak_width_list[2]]
        # right_values = [rt_values[b] for b in full_peak_width_list[3]]
        # highest_value = [rt_values[c] for c in peaks]

        left_values = [a for a in full_peak_width_list[2]]
        right_values = [b for b in full_peak_width_list[3]]
        highest_value = [c for c in peaks]

        self.rt_windows = [(left_values[i], right_values[i]) for i in range(len(left_values))]
        self.rt_peak_index = highest_value

    def _neutromer_peak_variance(self, should_plot=False):
        import scipy.ndimage as ndi
        import scipy.signal as sig
        import numpy as np

        intensity_filter = settings.intensity_filter

        sampling_rate = settings.sampling_rate
        smoothing_width = settings.smoothing_width
        smoothing_order = settings.smoothing_order

        peaks_list = self._get_peak_list(intensity_filter)

        # if max(peaks_list) < intensity_filter:
        #     return

        neutromer_maximums = list()
        gaussian_curves = list()
        
        for neutromer_peak in range(len(peaks_list)):
            neutromer_peak_list = peaks_list[neutromer_peak]

            from scipy.signal import savgol_filter
            smoothed_curves = savgol_filter(neutromer_peak_list, smoothing_width, smoothing_order)
            smoothed_curves[smoothed_curves < 0] = 0
            gaussian_curve = ndi.gaussian_filter1d(smoothed_curves, sampling_rate)
            gaussian_min_cuts = sig.find_peaks([-x for x in gaussian_curve])[0]
            peak_finding = np.concatenate(
                (np.array([gaussian_curve[1]]), gaussian_curve, np.array([gaussian_curve[-2]])))
    
            max_cuts = sig.argrelmax(peak_finding)[0]
            max_cuts = [item - 1 for item in max_cuts]
            gaussian_curves.append(gaussian_curve)
            neutromer_maximums.append(max_cuts)
        
        if should_plot:
            import matplotlib.pyplot as plt
    
            plot_colors = ['tab:blue', 'tab:orange', 'tab:purple', 'olivedrab', 'chocolate']
            gaussian_colors = ['b', 'orangered', 'indigo', 'lawngreen', 'saddlebrown']
    
            for i in range(self.n_isos):
                plt.plot(peaks_list[i], label=f"M{i}", color=plot_colors[i])
                plt.plot(gaussian_curves[i], label=f'Gaussian M{i}', zorder=2, color=gaussian_colors[i])
    
            plt.title(str(self.mz))
            plt.legend()
            plt.show()
        
        self.neutromer_peak_maximums = neutromer_maximums

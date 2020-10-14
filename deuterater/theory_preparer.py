'''Theoretical Value Preparation

This purpose of this module is to calculate the expected theoretical values
of each identified analyte.

'''

import pandas as pd
import numpy as np  # noqa: 401

import traceback  # noqa: 401
from pathlib import Path
from functools import partial
import multiprocessing as mp

from tqdm import tqdm  # noqa: 401

import deuterater.settings as settings
import utils.data_chunker as nvc



class TheoryPreparer():
    def __init__(self, enrichment_path, out_path, settings_path):
        settings.load(settings_path)
        self.settings_path = settings_path
        self.enrichment_path = Path(enrichment_path)
        if self.enrichment_path.suffix == '.tsv':
            self._enrichment_df = pd.read_csv(
                filepath_or_buffer=str(self.enrichment_path),
                sep='\t'
            )
        elif self.enrichment_path.suffix == '.csv':
            self._enrichment_df = pd.read_csv(
                filepath_or_buffer=str(self.enrichment_path),
                sep=','
            )
        if settings.recognize_available_cores is True:
            self._n_partitions = mp.cpu_count()
        else:
            self._n_partitions = settings.n_partitions

        self._mp_pool = mp.Pool(self._n_partitions)
        self.out_path = out_path
        self.model = None

    def write(self):
        self.model.to_csv(
            path_or_buf=self.out_path,
            sep='\t',
            index=False
        )

    def prepare(self):
        if settings.debug_level == 0:
            args_list = self._enrichment_df.to_records(index=False).tolist()
            func = partial(TheoryPreparer._mp_prepare, self.settings_path)
            results = list(
                tqdm(
                    self._mp_pool.imap_unordered(func, args_list),
                    total=len(self._enrichment_df)
                )
            
            )
            
            
        if settings.debug_level >= 1:
            print('Beginning single-processor theory preparation.')
            results = []
            for row in tqdm(self._enrichment_df.itertuples(),
                            total=len(self._enrichment_df)):
                # TODO: how to handle functions. Default I would think
                df = pd.read_csv(filepath_or_buffer=row.file, sep='\t')
                df = TheoryPreparer._apply_filters(df)
                df['time'] = row.time
                df['enrichment'] = row.enrichment
                df["sample_group"]  = row.sample_group
                results.append(df)

        self.model = pd.concat(results)
        self.model = self.model.drop(columns=['drop'])

        if settings.use_empir_n_value:
            self.model = nvc.data_chunker(self.model)

        self._mp_pool.close()
        self._mp_pool.join()

    @staticmethod
    def _mp_prepare(settings_path, args):
        settings.load(settings_path)
        #file_path, time, enrichment = args
        file_path, time, enrichment, sample_group = args
        df = pd.read_csv(filepath_or_buffer=file_path, sep='\t')
        df = TheoryPreparer._apply_filters(df)
        df['time'] = time
        df['enrichment'] = enrichment
        df["sample_group"]  = sample_group
        return df

    # def _load(self):
    #     '''Pulls in the relevant data from the model.

    #     Parameters
    #     ----------
    #     model : :obj:`pandas.Dataframe`
    #         The unmodified data model

    #     Returns
    #     -------
    #     :obj:`pandas.Dataframe`
    #         A dataframe containing a copy of the relevant columns

    #     '''
    #     # TODO: Basic checks like whether the data looks right need done
    #     if not isinstance(self.model, pd.DataFrame):
    #         try:
    #             # TODO: csv/tsv flexibility?
    #             self.model = pd.read_csv(self.model, sep='\t')
    #         except Exception as e:
    #             # TODO: better exception logging
    #             print(e)
    #             traceback.print_tb(e.__traceback__)
    #             raise
    #     try:
    #         # TODO: csv/tsv flexibility?
    #         self._enrichment_df = pd.read_csv(self.enrichment_path, sep='\t')
    #     except Exception as e:
    #         # TODO: better exception logging
    #         print(e)
    #         traceback.print_tb(e.__traceback__)
    #         raise

    @staticmethod
    def _apply_filters(df):
        '''filters the internal dataframe

        This function does not modify the dataframe in place.

        Parameters
        ----------
        df : :obj:`pandas.Dataframe`
            The internal dataframe of the theory_value_prep function

        Returns
        -------
        :obj:`pandas.Dataframe`
            The filtered dataframe. Does not modify in place.
        '''
        # This is 'clean_up_data' in the old deuterater
        # This is a
        data = df.dropna(
            axis='index',
            subset=['mzs', 'abundances']
        ).copy()
        data['drop'] = False
        for row in data.itertuples():
            mask = ((data['mz'] - row.mz).abs() <
                    settings.mz_proximity_tolerance)
            data.loc[mask, 'drop'] = True
        data = data[~data['drop']]

        # TODO: Check to see if no data went through
        return data


def main():
    print('please use the main program interface')


if __name__ == '__main__':
    main()

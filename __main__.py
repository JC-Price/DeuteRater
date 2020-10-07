'''
__main__.py
====================================
The the root module of the Deuterater Project.

'''
import argparse
from pathlib import Path
import os
import multiprocessing as mp

from PyQt5 import QtWidgets

from deuterater.extractor import Extractor
from deuterater.theory_preparer import TheoryPreparer
# TODO: integrate n_value calculator
from deuterater.fraction_new_calculator import FractionNewCalculator
from deuterater.rate_calculator import RateCalculator
from main_gui import MainGuiObject

# TODO: add logging capabilities


def main():
    '''
    The function to run if this file is called.
    Parameters
    ----------
    '''
    parser = argparse.ArgumentParser(
        prog='Deuterater',
        description='Python tool for analyzing labeled analytes'
    )
    subparsers = parser.add_subparsers(
        title='subcommands',
        dest='subcommand',
        description='valid subcommands'
    )
    #hard coding 
    parser.add_argument(
        '--settings', action='store', type=str,
        default='./resources\settings.yaml',
        help='The settings file to use for the extraction'
    )
    parser.add_argument(
        '--filters', action='store', type=str,
        default='./resources/filters.yaml',
        help='The settings file to use for the extraction'
    )

    # Set up the extractor subparser ##########################################
    parser_extract = subparsers.add_parser(
        name='extract',
        help='Extracts isotopic envelopes from mzML data'
    )
    parser_extract.add_argument(
        '--id', action='store',
        type=str, required=True,
        help='The .tsv/ .csv file with identification data'
    )
    parser_extract.add_argument(
        '--mzml', action='store',
        type=str, required=True,
        help='The .mzml data gathered from the laboratory'
    )
    parser_extract.add_argument(
        '--output', action='store',
        type=str, default=r'./dummy_data./extract_out.tsv',
        help='The output destination'
    )

    # Set up the directory extractor subparser ################################
    parser_extract_dir = subparsers.add_parser(
        name='extract_dir',
        help='Extracts isotopic envelopes from mzML data directories'
    )
    parser_extract_dir.add_argument(
        '--id', action='store',
        type=str, required=True,
        help='The .tsv/ .csv file with identification data'
    )
    parser_extract_dir.add_argument(
        '--mzml_dir', action='store',
        type=str, required=True,
        help='The .mzml data directory'
    )
    parser_extract_dir.add_argument(
        '--output_dir', action='store',
        type=str, default=r'./dummy_data/extract_out/',
        help='The output destination'
    )

    # Set up the theory subparser #############################################
    # TODO: this might be better named something else like: 'EnrichmentManager'
    parser_theory = subparsers.add_parser(
        name='theory',
        help='Manages enrichment data for rate analysis'
    )
    parser_theory.add_argument(
        '--enrich', action='store',
        type=str, required=True,
        help='the table of enrichment data'
    )
    parser_theory.add_argument(
        '--output', action='store',
        type=str, default=r'./dummy_data./extract_out.tsv',
        help='The output destination'
    )

    # TODO: Set up the nvalue subparser #######################################

    # Set up the fraction new subparser #######################################
    parser_fracnew = subparsers.add_parser(
        name='fracnew',
        help='Calculates fraction new values for rate analysis'
    )
    parser_fracnew.add_argument(
        '--model', action='store',
        type=str, required=True,
        help='The extracted data prepared with the enrichment'
    )
    parser_fracnew.add_argument(
        '--output', action='store',
        type=str, default=r'./dummy_data./extract_out.tsv',
        help='The output destination'
    )

    # Set up the rate calculator subparser ####################################
    parser_ratecalc = subparsers.add_parser(
        name='rate',
        help='Performs rate analysis'
    )
    parser_ratecalc.add_argument(
        '--model', action='store',
        type=str, required=True,
        help='The extracted data prepared with the enrichment and frac_new'
    )
    parser_ratecalc.add_argument(
        '--output', action='store',
        type=str, default=r'./dummy_data./extract_out.tsv',
        help='The output destination'
    )
    parser_ratecalc.add_argument(
        '--graph', action='store',
        type=str, default=r'./dummy_data./graph/',
        help='The graph file destination'
    )
    
    parser_gui = subparsers.add_parser(
        name = 'gui',
        help = 'activate the gui'
        )
    # Parse Commands ##########################################################
    # TODO: Should we add a try-catch block around this so that
    #       we can direct them to the usage information?
    args = parser.parse_args()

    if args.subcommand == 'extract':
        subargs, other_args = parser_extract.parse_known_args()
        extractor = Extractor(
            id_path=subargs.id,
            mzml_path=subargs.mzml,
            out_path=subargs.output,
            settings_path=args.settings
        )
        extractor.load()
        extractor.run()
        extractor.write()

    elif args.subcommand == 'extract_dir':
        subargs, other_args = parser_extract_dir.parse_known_args()
        mzml_dir_p = Path(subargs.mzml_dir).resolve()
        mzMLs = [mzml_dir_p / f for f in os.listdir(subargs.mzml_dir)]
        mzMLs = [p for p in mzMLs if p.suffix in ['.mzML', '.mzml']]
        print('')
        for file in mzMLs:
            extractor = Extractor(
                id_path=subargs.id,
                mzml_path=str(file),
                out_path=str(Path(subargs.output_dir) / f'{file.stem}.tsv'),
                settings_path=args.settings
            )
            extractor.load()
            print(file)
            extractor.run()
            extractor.write()
            print('')

    elif args.subcommand == 'theory':
        subargs, other_args = parser_theory.parse_known_args()
        theorist = TheoryPreparer(
            enrichment_path=subargs.enrich,
            out_path=subargs.output,
            settings_path=args.settings
        )
        theorist.prepare()
        theorist.write()

    # elif args.subcommand == 'nval':
    #     subargs, other_args = parser_nval.parse_known_args()
    #     print(subargs)
    #     raise NotImplementedError(
    #         'N Value calculator has yet to be fully integrated'
    #     )

    elif args.subcommand == 'fracnew':
        subargs, other_args = parser_fracnew.parse_known_args()
        fnewcalc = FractionNewCalculator(
            model_path=subargs.model,
            out_path=subargs.output,
            settings_path=args.settings
        )
        fnewcalc.generate()
        fnewcalc.write()

    elif args.subcommand == 'rate':
        subargs, other_args = parser_ratecalc.parse_known_args()
        ratecalc = RateCalculator(
            model_path=subargs.model,
            out_path=subargs.output,
            graph_folder = subargs.graph,
            settings_path=args.settings
        )
        ratecalc.calculate()
        ratecalc.write()
    elif args.subcommand == 'gui':
        import sys
        mp.freeze_support()
        app = QtWidgets.QApplication(sys.argv)
        gui_object = MainGuiObject(None)
        gui_object.show()
        app.exec_()
    # if args.subcommand == 'full':
    #     subargs, other_args = parser_full.parse_known_args()
    #     print(subargs)
    #     raise NotImplementedError(
    #         'Full run has yet to be implemented'
    #     )


if __name__ == '__main__':
    import time
    start_time = time.time()
    main()
    print('--- %s seconds ---' % (time.time() - start_time))
    print('')

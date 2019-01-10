#! /usr/bin/env python3

# ------------------------------------------------
# built-ins
import sys
import unittest
import argparse

# pipped
import pandas as pd
# ------------------------------------------------

def add_postgap(suite, postgap):
    """
    Iterate through suite and construct a new suite where
    each test has postgap passed as a keyword arg to it's
    constructor.
    """
    new_suite = unittest.TestSuite()

    for item in suite:
        test_class = item.__class__
        if test_class == unittest.TestSuite:
            new_suite.addTest(add_postgap(item, postgap))
        else:
            test_name = item._testMethodName
            new_suite.addTest(test_class(test_name, postgap))

    return new_suite


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='run the postgap unit tests')
    parser.add_argument('filename', type=str, help='a postgap output file (gzipped tsv file)')
    parser.add_argument('--skip-data-checks', default=False, action='store_true', help='skip the data checks')
    parser.add_argument('--skip-health-checks', default=False, action='store_true', help='skip the health checks')

    args = parser.parse_args()

    loader = unittest.TestLoader()

    pattern = 'test*.py'
    if (args.skip_data_checks):
        suite = loader.discover('./health_checks')
    elif (args.skip_health_checks):
        suite = loader.discover('./data_checks')
    else:
        suite = loader.discover('.')

    postgap = pd.read_csv(args.filename, sep='\t', na_values=['None'])

    suite_with_postgap = add_postgap(suite, postgap)
    result = unittest.TextTestRunner(verbosity=2).run(suite_with_postgap)
    sys.exit(not result.wasSuccessful())

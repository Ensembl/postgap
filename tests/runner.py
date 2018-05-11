#! /usr/bin/env python3

# ------------------------------------------------
# built-ins
import sys
import unittest

# pipped
import pandas as pd
import papermill as pm
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
    loader = unittest.TestLoader()
    suite = loader.discover('.')
    
    filename = sys.argv[1]
    postgap = pd.read_csv(filename, sep='\t', na_values=['None'])

    suite_with_postgap = add_postgap(suite, postgap)
    result = unittest.TextTestRunner(verbosity=2).run(suite_with_postgap)
    sys.exit(not result.wasSuccessful())

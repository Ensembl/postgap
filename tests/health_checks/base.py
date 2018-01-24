import unittest

class TestPostgapBase(unittest.TestCase):
    """
    Base class for postgap tests. Provides common utility methods.
    """

    def __init__(self, test_name, postgap=None):
        super(TestPostgapBase, self).__init__(test_name)
        self.pg = postgap

    def assert_series_in_range(self, series, low, high):
        """
        Check if all values in a `pandas.Series` are in the range [low, high].

        If any values fail the condition, the first exception is printed.
        """
        between = series.between(low, high)
        all_between = between.all()
        self.assertTrue(all_between,
                        series[~between].head(1).to_string(index=False))

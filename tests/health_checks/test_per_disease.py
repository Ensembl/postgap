# ------------------------------------------------
# built-ins
import unittest

# local
from .base import TestPostgapBase
# ------------------------------------------------

class TestPostgapPerDisease(TestPostgapBase):

    def setUp(self):
        self.per_disease = self.pg.groupby('disease_efo_id')

    def test_disease_efo_and_name_combination_unique(self):
        self.assert_groupby_series_is_unique_per_group(
            self.per_disease.disease_name
        )


if __name__ == '__main__':
    unittest.main()

# Testing POSTGAP output
The `tests` folder contains three types of quality control utilities for POSTGAP output, namely **health checks**, **data checks** and **reports**.

## Health checks
These are unit tests that:
* check the data schema, such as:
  * consistent format for all values in a column
  * uniqueness of values in a column, when grouped by values in another
* can be run against either a partial or whole output file, such as:
  * an output file for a single EFO term
  * a large concatenated output file (across all EFO terms)

## Data checks (TODO)
These are unit tests that:
* check biological expectations, such as:
  * filtering of the MHC region
  * filtering of *trans* relations (ie. when genes and snps have different chromosomes)
* can be run against a whole output file only

## Reports (TODO)
These are summary or metadata files that:
* can be generated for either a partial or whole output file
* present summary statistics to allow comparison between POSTGAP output files

# Usage

### Installation requirements
**Assumption**: You have `python3` installed.

Ideally, use `virtualenv` as follows:
```
# (in tests folder)
virtualenv -p python3 venv
source venv/bin/activate
pip install -r requirements.txt
```

### Run tests against a file
For a file in POSTGAP TSV format, run:
```
python runner.py ./sample_data/postgap.20180108.asthma.tsv.gz
```

## Usage

### Install requirements
Ideally, use `virtualenv` as follows:
```
# (in tests folder)
virtualenv venv
pip install -r requirements.txt
```

### Run tests against a file
For a file `asthma.tsv.gz` in POSTGAP TSV format, run:
```
python runner.py ./sample_data/postgap.20180108.asthma.tsv.gz
```

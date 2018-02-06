# ------------------------------------------------
# built-ins
import sys
import datetime

# pipped
import papermill as pm
# ------------------------------------------------

if __name__ == '__main__':
    now = datetime.datetime.now()
    filename = sys.argv[1]
    filestem = filename.split('/')[-1]
    ipynb_filename = '{}.REPORT.{}.ipynb'.format(filestem, now.strftime('%Y%m%d%H%M%S'))
    pm.execute_notebook(
        './reports/template.ipynb',
        './__reports__/{}'.format(ipynb_filename),
        parameters = dict(filename=filename)
    )

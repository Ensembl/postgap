#!/bin/bash
date

output_dir='results'
if [ ! -f $output_dir ]; then
  mkdir -p $output_dir
fi
python POSTGAP.py --diseases bmi --summary_stats summary_bmi.txt --bayesian --TYPE ML --output $output_dir/res


date

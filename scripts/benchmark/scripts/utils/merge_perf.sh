#!/bin/bash

cat perf_ncores*.txt | awk 'BEGIN {printf("app\tncore\trep\ttime\tmem\n");} ! /^app/' > perf.txt

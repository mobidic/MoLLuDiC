#!/usr/bin/env python3

import sys
import pandas

pandas.set_option('display.max_rows', 500)
pandas.set_option('display.max_columns', 500)
pandas.set_option('display.width', 1000)

hsm = sys.argv[1]

df = pandas.read_csv(hsm, sep='\t')

selectdf = df.filter(items=['PCT_PF_UQ_READS','ON_BAIT_VS_SELECTED','PCT_TARGET_BASES_10X','PCT_TARGET_BASES_50X','AT_DROPOUT','GC_DROPOUT','MEAN_INSERT_SIZE'])


print(selectdf.to_csv(sep='\t', index=False))

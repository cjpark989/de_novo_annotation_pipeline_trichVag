# some_file.py
import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '/Users/chungjunepark/Documents/Lab')

from methods import *

#rnd-1_family-33, unknown, and length greater than 1000: 1000.bed
sixpack_primer("unknown_fam_33")

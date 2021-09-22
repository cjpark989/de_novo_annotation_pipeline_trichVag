# some_file.py
import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '/Users/chungjunepark/Documents/Lab')

from methods import *

#rnd-1_family-33, unknown, and length greater than 1000: 1000.bed

#counter a will represent sequence_(a).seq 
a = 0

#going from 0 ~ 56 for .seq files
with  open("longest_orf.fa","a") as l:
    while a < 57:
        current_file_name = "/Users/chungjunepark/Documents/Lab/RP210909/LTR/sequence_" + str(a)

        with open(current_file_name) as f:
            first_line = f.readline().rstrip()
        
        l.write(first_line)
        l.write(longest_orf(current_file_name)) 
        l.write("\n")

        #increment a by 1
        a+=1

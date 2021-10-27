import sys
import os
strain = "T_Stableri_BTPI"
sys.path.insert(1,"/Users/chungjunepark/Documents/"+strain+"/methods")
from methods import *

directory = "/Users/chungjunepark/Documents/T_stableri_BTPI/files_and_outputs/inputs/input_for_bedtools" 

def pipeline(strain_name,n_bp,type_name):
    l = categorize(strain_name,type_name)
    print("done with categorizing!")

    # iterates through the list of family names, and picks the valid ones
    for e in l:
        filter_strain(strain_name,n_bp,e,type_name)
    print("done with filter_strain!")

    for fam_name in os.listdir(directory): 
        fam_name = fam_name.replace("input_for_bedtools_"+strain_name+"_","")
        fam_name = fam_name.replace(".bed","")
        bedtools_runner(strain_name,fam_name,type_name)
        sixpack_primer(strain_name,fam_name,type_name)
        sixpack_runner(strain_name,fam_name,type_name)
        longest_orf(strain_name,fam_name,type_name)
        print("one family done")
        print()

    print("done!")
n_bp = int(input("Specify length: "))
strain_name = input("Specify strain: ")
type_name = input("Specify type: ")

pipeline(strain,n_bp,type_name)








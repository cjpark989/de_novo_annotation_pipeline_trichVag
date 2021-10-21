import sys
import os

sys.path.insert(1,"/Users/chungjunepark/Documents/"+strain+"/methods")
from methods import *

directory = "/Users/chungjunepark/Documents/T_stableri_BTPI/files_and_outputs/inputs/input_for_bedtools" 

def pipeline(strain_name,n_bp,type_name):
    l = categorize(strain_name,type_name)
    print("done with categorizing!")
    for e in l:
        filter_strain(strain_name,n_bp,e,type_name)

    print("done with filter_strain!")
    for filename in os.listdir(directory): 
        filename = filename.replace("input_for_bedtools_"+strain_name+"_","")
        filename = filename.replace(".bed","")
        bedtools_runner(strain_name,filename,type_name)
        sixpack_primer(strain_name,filename,type_name)
        sixpack_runner(strain_name,filename,type_name)
        #longest_orf(strain_name,filename,type_name)
n_bp = int(input("Specify length: "))
strain_name = input("Specify strain: ")
type_name = input("Specify type: ")

pipeline(strain,n_bp,type_name)


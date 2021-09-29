import sys
strain = "T_stableri_BTPI"
sys.path.insert(1,"/Users/chungjunepark/Documents/"+strain+"/methods")
from methods import *

def pipeline(strain_name,n_bp,fam_name,type_name):
    filter_strain(strain_name,n_bp,fam_name,type_name)
    bedtools_runner(strain_name,fam_name,type_name)
    sixpack_primer(strain_name,fam_name,type_name)
    sixpack(strain_name,fam_name,type_name)

def main():
    pipeline("T_stableri_BTPI",1000,"rnd-1_family-33","Unknown")

main()
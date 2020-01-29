import os

ROOT_FOLDER="/specific/netapp5/gaga/hagailevi/aneuploidy/" # "/media/hag007/Data/aneuploidy" # "/specific/netapp5/gaga/hagailevi/aneuploidy/"   #"D:\\aneuploidy"
CODE_FOLDER=os.path.join(ROOT_FOLDER,"code")
CACHE_FOLDER=os.path.join(ROOT_FOLDER,"cache")
OUTPUT_FOLDER=os.path.join(ROOT_FOLDER,"output")
DATASETS_FOLDER=os.path.join(ROOT_FOLDER,"datasets")

INTRA=0
INTER=1

CHR_INTERACTION_NAMES={INTRA:"intra",
                   INTER: "inter"}

#!/usr/bin/python

import re
import sys
import shutil

file=sys.argv[1]
gse=sys.argv[2]
print("Renaming " + file)

p = re.compile(".*(GSM\d{4,9}).*")
m = p.match(file)
dest = m.group(1)

path = "/fslhome/aguyer/fsl_groups/fslg_Bio465_BDAG/geneExpression/"+gse+"/unprocessed/"

source=path+file
destination=path+dest+".CEL"

print("Previous Name: " + source)
print("Destination Name: " + destination)
print(" ")

shutil.move(source, destination)

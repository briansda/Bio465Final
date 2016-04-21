import os, sys, glob
import utilities

inFilePathPattern = sys.argv[1]
outFilePath = sys.argv[2]

#print(inFilePathPattern)
#print glob.glob(os.path.expanduser(inFilePathPattern))
#exit()

uniqueFeatures = None
sampleIDs = []

for inFilePath in glob.glob(os.path.expanduser(inFilePathPattern)):
    print "Getting feature list from %s" % inFilePath
    sampleIDs.append(os.path.basename(inFilePath))

    inData = utilities.readMatrixFromFile(inFilePath)
    inData.pop(0)
    features = set([row[0] for row in inData])

    if uniqueFeatures == None:
        uniqueFeatures = features
    else:
        uniqueFeatures = uniqueFeatures & features

allDataDict = {}
#print(uniqueFeatures)

for inFilePath in glob.glob(os.path.expanduser(inFilePathPattern)):
    print "Getting data from %s" % inFilePath
    inData = utilities.readMatrixFromFile(inFilePath)
    inData.pop(0)

    inDataDict = {}
    for row in inData:
        inDataDict[row[0]] = row[1]

    allDataDict[os.path.basename(inFilePath)] = inDataDict

outFile = open(outFilePath, 'w')
outFile.write("\t".join([""] + sampleIDs) + "\n")
for feature in uniqueFeatures:
    print "Writting data for %s" % feature
    outFile.write("\t".join([feature] + [allDataDict[sampleID][feature] for sampleID in sampleIDs]) + "\n")
outFile.close()

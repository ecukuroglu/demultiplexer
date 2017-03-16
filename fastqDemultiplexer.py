#################
## Created by Engin Cukuroglu
## 10 March 2017
#################

import os, sys, time
import gzip, io
from itertools import combinations
from itertools import product

### write reads in the memory to the files
def writeBufferToFile(barcodeFileDict, stringBufferDict):
    for barcode in barcodeFileDict:
        barcodeFileDict[barcode].write('%s' %(stringBufferDict[barcode]))
        stringBufferDict[barcode] = ''
    return stringBufferDict


### read zipped fastq file
### checked the barcode in barcode dictionary
### write reads to assigned files
def demultiplexer(originalFastqFileDirectory, barcodeDict, barcodeFileNameDict, barcodeLength):
    barcodeFileDict = {}
    tempStringBufferDict = {}
    for barcode, tempFileDirectory in barcodeFileNameDict.items():
        barcodeFileDict[barcode] = io.BufferedWriter(gzip.open(tempFileDirectory, mode='wb'))
        tempStringBufferDict[barcode] = ''
    readCounter = 0
    with io.BufferedReader(gzip.open(originalFastqFileDirectory, 'rb')) as fin:
        for line in fin:
            line = line.strip()
            barcode = line[len(line)-barcodeLength:]
            if barcode in barcodeDict:
                barcode = barcodeDict[barcode]
            else:
                barcode = 'unknown'
            #print(barcode)
            tempLine = '%s\n' %(line)
            line = fin.readline()
            tempLine = '%s%s' %(tempLine, line)
            line = fin.readline()
            tempLine = '%s%s' %(tempLine, line)
            line = fin.readline()
            tempLine = '%s%s' %(tempLine, line)
            barcodeFileDict[barcode].write('%s' %(tempLine))
            #tempStringBufferDict[barcode] = '%s%s' %(tempStringBufferDict[barcode], tempLine)
            readCounter = readCounter + 1
            if readCounter == 50000:
                print(readCounter)
                #tempStringBufferDict = writeBufferToFile(barcodeFileDict, tempStringBufferDict)
                readCounter = 0
        #tempStringBufferDict = writeBufferToFile(barcodeFileDict, tempStringBufferDict)
    for barcode in barcodeFileDict.keys():
        print('barcode %s' %(barcode))
        barcodeFileDict[barcode].close()
            

### allow matching unread base indexes with actual barcodes
def barcodeCombinationGenerator_withUnreadBase(barcode, stationaryIndex, maxUnreadBase):
    changeableIndexes = []
    for i in list(range(len(barcode))):
        if not i in stationaryIndex: 
            changeableIndexes.append(i)
    barcodeDict = {}
    barcodeDict[barcode] = barcode
    if maxUnreadBase > 0:
        for i in list(range(1, maxUnreadBase+1)):
            for nIndexes in list(combinations(changeableIndexes, i)):
                tempString = list(barcode)
                for j in nIndexes:
                    tempString[j] = 'N'
                barcodeDict[''.join(tempString)] = barcode
    return barcodeDict
        
### allow mismatch in actual barcodes and generate a lookup dictionary
def barcodeCombinationGenerator_withMissingBase(barcode, stationaryIndex, missedMatchNumber):
    baseAlphabet = ['A', 'T', 'C', 'G', 'N']
    changeableIndexes = []
    for i in list(range(len(barcode))):
        if not i in stationaryIndex: 
            changeableIndexes.append(i)
    barcodeDict = {}
    barcodeDict[barcode] = barcode
    if missedMatchNumber > 0:
        for i in list(range(1,missedMatchNumber+1)):
            for nIndexes in list(combinations(changeableIndexes, i)):
                for tempBase in list(product(baseAlphabet, repeat=i)):
                    tempString = list(barcode)
                    for j in list(range(len(nIndexes))):
                        tempString[nIndexes[j]] = tempBase[j]
                    barcodeDict[''.join(tempString)] = barcode
    return barcodeDict


### Merge multiple barcode dictionaries and label the conflicts
def barcodeCombinationMerger(barcodeList,
                stationaryIndex_missed,
                stationaryIndex_unread,
                missedMatchNumber,
                maxUnreadBase):
    barcodeDictList = {}
    keyList = []
    for barcode in barcodeList:
        barcodeDictList[barcode] = barcodeCombinationGenerator_withMissingBase(barcode,
                                                                               stationaryIndex_missed,
                                                                               missedMatchNumber)
        barcodeDictList[barcode].update(barcodeCombinationGenerator_withUnreadBase(barcode,
                                                                                   stationaryIndex_unread,
                                                                                   maxUnreadBase))
        keyList = keyList + list(barcodeDictList[barcode].keys())
    dummySet = set()
    dummySet_add = dummySet.add
    duplicateElements = set(x for x in keyList if x in dummySet or dummySet_add(x))
    counter = 0
    barcodeDict = {}
    for barcode in barcodeDictList:
        if counter == 0:
            barcodeDict = barcodeDictList[barcode]
            counter = counter + 1
        else:
            barcodeDict.update(barcodeDictList[barcode])
    for tempElement in duplicateElements:
        barcodeDict[tempElement] = 'conflict'
    for barcode in barcodeList:
        barcodeDict[barcode] = barcode
    return barcodeDict
            
def __main__():
    maxUnreadBase = 5
    missedMatchNumber = 4
    stationaryIndex_missed = [6,7]
    stationaryIndex_unread = [6,7]
    barcodeLength = 8
   
    originalFastqFileDirectory = ''
    
    finalFastqFileDepository = 'finalFastq_Dummy_v3'
    if not os.path.exists(finalFastqFileDepository):
        os.makedirs(finalFastqFileDepository)
    barcodeFileNameDict = {}
    barcodeFileNameDict['unknown'] = '%s/undetermined_S50_L005_R1_001.fastq.gz' %(finalFastqFileDepository)
    barcodeFileNameDict['conflict'] = '%s/conflict_S50_L005_R1_001.fastq.gz' %(finalFastqFileDepository)
    barcodeFileNameDict['CGATGTCA'] = '%s/DHE117-CGATGTCA_S50_L005_R1_001.fastq.gz' %(finalFastqFileDepository)
    barcodeFileNameDict['ATCACGTG'] = '%s/DHE118-ATCACGTG_S50_L005_R1_001.fastq.gz' %(finalFastqFileDepository)
    barcodeList = ['CGATGTCA', 'ATCACGTG']
    barcodeDict = barcodeCombinationMerger(barcodeList,
                        stationaryIndex_missed,
                        stationaryIndex_unread,
                        missedMatchNumber,
                        maxUnreadBase)
    demultiplexer(originalFastqFileDirectory, barcodeDict, barcodeFileNameDict, barcodeLength)

t1 = time.time()
__main__()
t2 = time.time()
print(t2-t1)


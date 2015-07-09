# -*- coding: utf-8 -*-

import csv

import replicateHandling as rh
from DsfPlate import DsfPlate
from MeanWell import MeanWell

class DsfAnalysis:
    def __init__(self, analysisName):
        #initialisations
        self.name = analysisName
        self.plate = None
        self.meanWells = []
        self.contentsHash = {}
        return
        
    def loadCurves(self, dataFilePath, contentsMapFilePath):
        #create the DsfPlate object
        self.plate = DsfPlate(dataFilePath, contentsMapFilePath)
        return
    
    def analyseCurves(self):
        #perform analysis on the plate, note order here is important
        self.plate.computeOutliers()
        self.plate.computeSaturations()
        self.plate.computeMonotonicities()
        self.plate.computeInTheNoises()
        self.plate.computeTms()
        self.plate.computeComplexities()
        #create the mean wells of replicates on the plate
        self.__createMeanWells()
        #create grouped hash for plotting
        self.__createMeanContentsHash()
        return
    
    def __createMeanWells(self):
        seen = []
        #loop through each set of replicates
        for wellName in self.plate.repDict.keys():
            if wellName not in seen:
                reps = self.plate.repDict[wellName]
                seen += reps
                #get mean tm and tm error (sd of tms)
                tm, tmError = rh.meanSd([self.plate.wells[w].Tm for w in reps])
                complexMean = any([self.plate.wells[w].isComplex for w in reps])
                contents = self.plate.wells[wellName].contents
                #create a mean well and add it to list
                self.meanWells.append(MeanWell(tm, tmError, complexMean, reps, contents))
        return
    
    def __createMeanContentsHash(self):
        #loop through each mean well
        for well in self.meanWells:
            contents = well.contents
            #build the nested hashmaps as we go through the wells, start with cv1
            if contents.cv1 not in self.contentsHash.keys():
                self.contentsHash[contents.cv1] = {}
            #next its ph
            if contents.ph not in self.contentsHash[contents.cv1].keys():
                self.contentsHash[contents.cv1][contents.ph] = {}
            #then cv2, which maps to the mean well itself
            if contents.cv2 not in self.contentsHash[contents.cv1][contents.ph].keys():
                self.contentsHash[contents.cv1][contents.ph][contents.cv2] = well
        return
    
    def produceNormalisedOutput(self, filePath):
        #names of all the wells in sorted order
        sortedWellNames = sorted(self.plate.wells.keys())
        #list of temperatures, taken from first well since all have the same temperature list
        temperatures = self.plate.wells.values()[0].temperatures
        
        with open(filePath, 'w') as fp:
            fWriter = csv.writer(fp, delimiter='\t')
            fWriter.writerow(['Temperature'] + sortedWellNames)
            for i in range(len(temperatures)):
                #start each row with the temperature
                row = [temperatures[i]]
                #create each row as the value at that temperature on each well
                for wellName in sortedWellNames:
                    row.append(self.plate.wells[wellName].fluorescence[i])
                
                #write to the file
                fWriter.writerow(row)
        return
    
    def generateReport(self, outputFilePath):
        return
    
    
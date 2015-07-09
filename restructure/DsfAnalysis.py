# -*- coding: utf-8 -*-

import csv
import os
import pandas as pd

import replicateHandling as rh
from DsfPlate import DsfPlate, LYSOZYME, SIMILARITY_THRESHOLD
from MeanWell import MeanWell

#(mean, standard deviation) of lysozyme Tm over ~250 experiments
LYSOZYME_TM_THRESHOLD = (70.8720, 0.7339)

#the running location of this file
RUNNING_LOCATION = os.path.dirname(os.path.realpath(__file__))

class DsfAnalysis:
    def __init__(self, analysisName):
        #initialisations
        self.name = analysisName
        self.plate = None
        self.meanWells = []
        self.contentsHash = {}
        self.controlsHash = {}
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
        #check the controls on the plate
        self.__createControlsHash()
        return
    
    def __createMeanWells(self):
        seen = []
        #loop through each set of replicates
        for wellName in self.plate.repDict.keys():
            if wellName not in seen:
                reps = self.plate.repDict[wellName]
                seen += reps
                #get mean tm and tm error (sd of tms)
                tm, tmError = rh.meanSd([self.plate.wells[w].tm for w in reps])
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
    
    def __createControlsHash(self):
        #initialise the control results
        results = {"lysozyme": "Not Found",
                   "no dye": "Not Found",
                   "no protein": "Not Found"}
        
        #first check if lysozyme control is present on the plate
        if len(self.plate.lysozyme)>0:
            #get the mean well for the lysozyme control, it will have no ph, and no condition variable 2
            lysozymeMeanWell = self.contentsHash[LYSOZYME]['']['']
            
            #lysozyme Tm check, only uses mean lysozyme Tm, hence indexing ([0])
            if lysozymeMeanWell.tm > LYSOZYME_TM_THRESHOLD[0] - 2*LYSOZYME_TM_THRESHOLD[1] and\
            lysozymeMeanWell.tm < LYSOZYME_TM_THRESHOLD[0] + 2*LYSOZYME_TM_THRESHOLD[1]:
                results["lysozyme"] = "Passed"
            else:
                results["lysozyme"] = "Failed"
        
        #check if no dye control is present
        if len(self.plate.noDye)>0:
            #create a mean curve out of the replicates that are not outliers
            meanNoDyeCurve = [0 for x in self.plate.wells[self.plate.noDye[0]]]
            validCurvesInSum = 0
            for wellName in self.plate.noDye:
                well = self.plate.wells[wellName]
                #creates sum of curves being used
                if not well.isOutlier:
                    meanNoDyeCurve = [x+y for x,y in zip(well.fluorescence, meanNoDyeCurve)]
                    validCurvesInSum += 1
            #divide sum to give average curve
            meanNoDyeCurve = [x/validCurvesInSum for x in meanNoDyeCurve]
            
            #read in expected no dye control from file
            noDyeExpected = list(pd.Series.from_csv(RUNNING_LOCATION + "\\data\\noDyeControl.csv"))
            #if the curves are within required distance from one another, the control is passed
            if rh.sqrdiff(meanNoDyeCurve, noDyeExpected) < SIMILARITY_THRESHOLD:
                results["no dye"] = "Passed"
            else:
                results["no dye"] = "Failed"
        
        #check if no protein control is present
        if len(self.plate.noProtein)>0:
            #create a mean curve out of the replicates that are not outliers
            meanNoProteinCurve = [0 for x in self.plate.wells[self.plate.noDye[0]]]
            validCurvesInSum = 0
            for wellName in self.plate.noDye:
                well = self.plate.wells[wellName]
                #creates sum of curves being used
                if not well.isOutlier:
                    meanNoProteinCurve = [x+y for x,y in zip(well.fluorescence, meanNoProteinCurve)]
                    validCurvesInSum += 1
            #divide sum to give average curve
            meanNoProteinCurve = [x/validCurvesInSum for x in meanNoProteinCurve]
            
            #read in the expected curve for the no protein control
            noProteinExpected = list(pd.Series.from_csv(RUNNING_LOCATION + "\\data\\noProteinControl.csv"))
            #if the curves are within required distance from one another, the control is passed
            if rh.sqrdiff(meanNoProteinCurve, noProteinExpected) < SIMILARITY_THRESHOLD:
                results["no protein"] = "Passed"
            else:
                results["no protein"] = "Failed"
        
        #save the results hash
        self.controlsHash = results
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
    
    
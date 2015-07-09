# -*- coding: utf-8 -*-

import pandas as pd
from itertools import combinations

import replicateHandling as rh
from DsfWell import DsfWell
from Contents import Contents
from MeltdownException import MeltdownException

#control names that we recognise in contents map (lower cased)
LYSOZYME = 'lysozyme'
NO_DYE = 'no dye'
PROTEIN_AS_SUPPLIED = 'protein as supplied'
NO_PROTEIN = 'no protein'

#possible colours of plots
COLOURS = ["Blue","DarkOrange","Green","Magenta","Cyan","Red",
            "DarkSlateGray","Olive","LightSeaGreen","DarkMagenta","Gold","Navy",
            "DarkRed","Lime","Indigo","MediumSpringGreen","DeepPink","Salmon",
            "Teal","DeepSkyBlue","DarkOliveGreen","Maroon","GoldenRod","MediumVioletRed"]
            
#discarding bad replicates threshold. calculated as mean difference between any two of 168 normalised lysozyme curves
SIMILARITY_THRESHOLD = 1.72570084974
#gives threshold for monotonicity in a non normalised melt curve when multiplied by highest fluorescence value on the plate
PLATE_MONOTONICITY_THRESHOLD_FACTOR = 0.0005
#gives the 'in the noise' threshold when multiplied by the mean monotonicity threshold of the 'no protein' control wells
NOISE_THRESHOLD_FACTOR = 1.15


class DsfPlate:
    def __init__(self, dataFilePath, contentsMapFilePath):
        
        #initialise dict of well names to wells
        self.wells = {}
        #initialise lists of control names
        self.lysozyme = []
        self.noDye = []
        self.proteinAsSupplied = []
        self.noProtein = []
        #initialise mapping of condition variable 2's to colour
        self.cv2ColourDict = {}
        #initialised replicate dictionary, maps well name to name of all its replicate wells (including itself)
        self.repDict = {}
        
        #initial values for plate specific thresholds
        self.plateMonotonicThreshold = None
        self.noiseThreshold = None
        
        #==================read the data file as a pandas data frame
        try:
            data = pd.DataFrame.from_csv(dataFilePath, sep='\t', index_col='Temperature')
        except Exception as e:
            raise MeltdownException('There was a problem reading the data file\n' + e.message)
        
        #turn the column headers (well names) into strings, allows naming well names 1,2,3 etc
        data.columns = [str(colName) for colName in data.columns]
        #remove any columns that that are blank (default pcrd export includes empty columns sometimes)
        for column in data:
            if 'Unnamed' in column:
                data.pop(column)
        #remove any rows that are all blank (e.g. empty lines at the end of the file)
        data.dropna(how='all', inplace=True)
        #replace any empty cells (default value NaN) to be empty strings ('')
        data.fillna(value='', inplace=True)

        #==================read the in the contents map as a dataframe too
        try:
            contentsMap = pd.DataFrame.from_csv(contentsMapFilePath, sep='\t', index_col='Well')
        except Exception as e:
            raise MeltdownException('There was a problem reading the contents map file\n' + e.message)
        
        #remove any columns that that are blank (default pcrd export includes empty columns sometimes)
        for column in contentsMap:
            if 'Unnamed' in column:
                contentsMap.pop(column)
        #remove any rows that are all blank (e.g. empty lines at the end of the file)
        contentsMap.dropna(how='all', inplace=True)
        #replace any empty cells (default value NaN) to be empty strings ('')
        contentsMap.fillna(value='', inplace=True)
        #==================
        
        for wellName in data.columns:
            wellContents = self.__readContentsOfWell(contentsMap, wellName)
            #check if well is one fo the 4 supported controls, and add name to appropriate list if that is the case
            if wellContents.cv1.lower() == LYSOZYME:
                self.lysozyme.append(wellName)
                wellContents.control = 1
            elif wellContents.cv1.lower() == NO_DYE:
                self.noDye.append(wellName)
                wellContents.control = 1
            elif wellContents.cv1.lower() == PROTEIN_AS_SUPPLIED:
                self.proteinAsSupplied.append(wellName)
                wellContents.control = 1
            elif wellContents.cv1.lower() == NO_PROTEIN:
                self.noProtein.append(wellName)
                wellContents.control = 1
            
            #populate the list of wells
            self.__addWell(data[wellName], wellName, wellContents)
        
        #create a mapping of condition variable 2's to particular colours, to help with plotting
        self.__assignConditionVariable2Colours(contentsMap)
        #create a mapping of each well name to a list of wellnames that are replicates of itself (includeing itself)
        self.__createRepDict(contentsMap)
                
        return
    
    def __readContentsOfWell(self, contentsMap, wellName):
        #get the well's row from the dataframe
        contentsRow = contentsMap.xs(wellName)
        #get information from the row
        try:
            cv1 = contentsRow['Condition Variable 1']
        except Exception as e:
            raise MeltdownException('Could not read "Condition Variable 1" column for "' + wellName + '" from contents map\n' + e.message)
        
        try:
            cv2 = contentsRow['Condition Variable 2']
        except Exception as e:
            raise MeltdownException('Could not read "Condition Variable 2" column for "' + wellName + '" from contents map\n' + e.message)
        
        #as ph, dphdt, and control columns are not essential, if they are ommited, values take empty strings
        try:
            ph = contentsRow['pH']
        except KeyError:
            ph = ''
        except Exception as e:
            raise MeltdownException('Could not read "pH" column for "' + wellName + '" from contents map\n' + e.message)
        
        try:
            dphdt = contentsRow['d(pH)/dT']
        except KeyError:
            dphdt = ''
        except Exception as e:
            raise MeltdownException('Could not read "d(pH)/dT" column for "' + wellName + '" from contents map\n' + e.message)
        
        try:
            control = contentsRow['Control']
        except KeyError:
            control = ''
        except Exception as e:
            raise MeltdownException('Could not read "Control" column for "' + wellName + '" from contents map\n' + e.message)
        
        #create Contents object for the well
        contents = Contents(cv1, cv2, ph, dphdt, control)
        return contents
    
    def __addWell(self, fluorescenceSeries, name, contents):
        #create a dsf well object and add it to wells list
        well = DsfWell(list(fluorescenceSeries.values), list(fluorescenceSeries.index), name, contents)
        self.wells[name] = well
        return
    
    def __assignConditionVariable2Colours(self, contentsMap):
        colourIndex=0
        #for each unseen condition variable 2, map the next colour in the COLOUR list
        for cv2 in contentsMap['Condition Variable 2']:
            if cv2 not in self.cv2ColourDict.keys():
                #limit to how many condition variable 2's there can be
                if colourIndex == len(COLOURS):
                    raise MeltdownException('No more than '+ str(len(COLOURS)) +'different Condition Variable 2\'s are permitted')
                self.cv2ColourDict[cv2] = COLOURS[colourIndex]
                colourIndex += 1
        return
        
    def __createRepDict(self, contentsMap):
        #loop through each well name
        for wellName in contentsMap.index:
            self.repDict[wellName] = []
            #reloop through each well name, and add those with the same solution to the initial wells list of replicates
            for comparedWellName in contentsMap.index:
                well = contentsMap.xs(wellName)
                cmpWell = contentsMap.xs(comparedWellName)
                #replicate defined as having same condition variables 1 and 2 as well as same ph
                if well['Condition Variable 1'] == cmpWell['Condition Variable 1'] and\
                well['Condition Variable 2'] == cmpWell['Condition Variable 2'] and\
                well['pH'] == cmpWell['pH']:
                    #found a replicate
                    self.repDict[wellName].append(comparedWellName)
        return
    
    def computeOutliers(self):
        seen = []
        outlierWells = []
        for wellName in self.wells.keys():
            #carful not to loop over the same wells
            if wellName not in seen:
                reps = self.repDict[wellName]
                seen += reps
                #every possible pair of replicates
                pairs = combinations(reps, 2)
                #initialise the distance matrix, as described in replicate handling
                distMatrix = [[0 for x in range(len(reps))] for y in range(len(reps))]
                for pair in pairs:
                    #calculate the distance between the pairs
                    dist = rh.sqrdiff(self.wells[pair[0]].fluorescence, self.wells[pair[1]].fluorescence)
                    #add the appropriate distance values to the distance matrix
                    distMatrix[reps.index(pair[0])][reps.index(pair[1])] = dist
                    distMatrix[reps.index(pair[1])][reps.index(pair[0])] = dist
                #get list of replicates which are NOT outliers
                keep = rh.discardBad(reps, distMatrix, SIMILARITY_THRESHOLD)
                #add to the total list of outlier wells
                for rep in reps:
                    if rep not in keep:
                        outlierWells.append(wellName)
        #go thraough all the outlier wells and set their outlier and discarded flags to true
        for wellName in outlierWells:
            self.wells[wellName].isOutlier = True
            self.wells[wellName].isDiscarded = True
        return
    
    def computeSaturations(self):
        #let each stored well calculate if it is saturated
        for well in self.wells.values():
            well.computeSaturation()
        return
    
    def computeMonotonicities(self):
        self.__computePlateMonotonicThreshold()
        #make each well calculate whether it is monotonic, with the now calculated plate monotonicity threshold
        for well in self.wells.values():
            well.computeMonotonicity(self.plateMonotonicThreshold)
        return
    
    def computeInTheNoises(self):
        self.__computeNoiseThreshold()
        #make each well calculate whether it is in the noise, with the now calculated noise threshold
        for well in self.wells.values():
            well.computeInTheNoise(self.noiseThreshold)
        return
    
    def computeTms(self):
        #each well calculates its Tm on itself
        for well in self.wells.values():
            well.computeTm()
        return
    
    def computeComplexities(self):
        #each well calculates if it is complex on itself
        for well in self.wells.values():
            well.computeComplexity()
        return
    
    def __computePlateMonotonicThreshold(self):
        #get the highest fluorescence value from all wells before they were normalised
        overallMaxNonNormalised = 0
        for well in self.wells.values():
            if well.wellMax > overallMaxNonNormalised:
                overallMaxNonNormalised = well.wellMax
        #calculate the plates monotonic threshold used the constant factor
        self.plateMonotonicThreshold = PLATE_MONOTONICITY_THRESHOLD_FACTOR * overallMaxNonNormalised
        return
    
    def __computeNoiseThreshold(self):
        #if no no protein controls, leave the noise threshold as None, and let this be handled in DsfWell
        if len(self.noProtein)==0:
            self.noiseThreshold = None
            return
        #otherwise, calculate th noise threshold from the constant factor
        meanNoProteinMonotonicThreshold, sd = rh.meanSd([self.wells[wellName].wellMonotonicThreshold for wellName in self.noProtein])
        self.noiseThreshold = meanNoProteinMonotonicThreshold * NOISE_THRESHOLD_FACTOR
        return

















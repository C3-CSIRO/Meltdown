# -*- coding: utf-8 -*-

import pandas as pd

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

class DsfPlate:
    def __init__(self, dataFilePath, contentsMapFilePath):
        
        #initialise list of wells
        self.wells = []
        #initialise lists of control names
        self.lysozyme = []
        self.noDye = []
        self.proteinAsSupplied = []
        self.noProtein = []
        #initialise mapping of condition variable 2's to colour
        self.cv2ColourDict = {}
        #initialised replicate dictionary, maps well name to name of all its replicate wells (including itself)
        self.repDict = {}
        
        
        #read the data file as a pandas data frame
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
        
        
        
        #read the in the contents map as a dataframe too
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
        self.wells.append(well)
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
    
    def computeOutliers():
        return
    
    def computeSaturations():
        return
    
    def computeMonotonicities():
        return
    
    def computeInTheNoises():
        return
    
    def computeTms():
        return
    
    def computeComplexities():
        return
    
    def __computePlateMonotonicThreshold():
        return
    
    def __computeNoiseThreshold():
        return
    
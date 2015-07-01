# -*- coding: utf-8 -*-

import pandas as pd

class DsfPlate:
    def __init__(self, dataFilePath, contentsMapFilePath):
        
        #read the data file as a pandas data frame
        data = pd.DataFrame.from_csv(dataFilePath, sep='\t', index_col='Temperature')
        #turn the column headers (well names) into strings, allows naming well names 1,2,3 etc
        data.columns = [str(colName) for colName in data.columns]
        #remove any columns that that are blank (default pcrd export includes empty columns sometimes)
        for column in data:
            if 'Unnamed' in column:
                data.pop(column)
        #remove any rows that are all blank (e.g. empty lines at the end of the file)
        data.dropna(how='all', inplace=True)
        
        #read the in the contents map as a dataframe too
        contentsMap = pd.DataFrame.from_csv(contentsMapFilePath, sep='\t', index_col=None)
        #remove any columns that that are blank (default pcrd export includes empty columns sometimes)
        for column in contentsMap:
            if 'Unnamed' in column:
                contentsMap.pop(column)
        #remove any rows that are all blank (e.g. empty lines at the end of the file)
        contentsMap.dropna(how='all', inplace=True)
        
        for wellName in data.columns:
            wellContents = self.__readContentsOfWell(contentsMap, wellName)
        
                
        return
    
    def __readContentsOfWell(contentsMap, wellName):
        contents = 'tmp'
        return contents
    
    def __addWell(fluorescence, temperatures, name, contents):
        return
    
    def __assignConditionVariable2Colours(contentsMap):
        return
        
    def createRepDict(contentsMap):
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
    
# -*- coding: utf-8 -*-

import math

class DsfWell:
    def __init__(self,fluorescence,temperatures,name,contents):
        self.fluorescence = fluorescence
        self.temperatures = temperatures
        self.contents = contents
        self.name = name

        self.normalisationFactor = None
        self.wellMax = None
        self.wellMin = None
        self.wellNormalisedMax = None 
        self.wellNormalisedMin = None
        self.tm = None
        self.wellMonotonicThreshold = None

        self.isMonotonic = False
        self.isComplex =  False
        self.isOutlier = False
        self.isInTheNoise = False
        self.isSaturated = False
        self.isDiscarded = False
        
        #get min and max of the non normalised and normalised curves
        self.wellMin, self.wellMax = self.getMinAndMax()
        self.normalise()
        self.wellNormalisedMin, self.wellNormalisedMax = self.getMinAndMax()

        return
    
    def getMinAndMax(self):
        minimum = self.fluorescence[0]
        maximum = self.fluorescence[0]
        for value in self.fluorescence[1:]:
            if value > maximum:
                maximum = value
            if value < minimum:
                minimum = value
        return (minimum,maximum)
    
    def normalise(self):
        #throw error if there is only 1 fluorescence value?
        stepSize = math.fabs(self.temperatures[1] - self.temperatures[0])
        count = 0
        for height in self.fluorescence:
            count += height*stepSize
        self.fluorescence = [x / count for x in self.fluorescence]
        #used to calculate the monotenicity threshold
        self.normalisationFactor = count
        return
    
    def computeSaturation(self):
        if not self.isDsicarded:
            currentValue = self.fluorescence[0]
            maxInd = 0
            while currentValue != self.wellNormalisedMax:
                maxInd+=1
                currentValue = self.fluorescence[maxInd]

            diff = self.wellNormalisedMax - self.wellNormalisedMin

            # A boundry defining how much the points can fluctuate and still be considered flat
            lowFlatBoundry = self.wellNormalisedMax - 0.005*diff

            # Look each way to see how many temperature steps the curve stays flat for
            count = 0
            ind = maxInd - 1
            while ind>=0:
                if self.fluorescence[ind] > lowFlatBoundry:
                    count += 1
                    ind -= 1
                else:
                    break
            ind = maxInd+1
            while ind<len(self.fluorescence):
                if self.fluorescence[ind] > lowFlatBoundry:
                    count += 1 
                    ind += 1
                else:
                    break
            if count >= 10:
                self.isSaturated = True
                self.isDiscarded = True

        return
    
    def computeMonotonicity(plateMonotonicThreshold):
        #set self.wellMonotonicThreshold here
        return
    
    def computeInTheNoise(noiseThreshold):
        return
    
    def computeTm():
        return
    
    def computeComplexity():
        return
    
# -*- coding: utf-8 -*-


class DsfWell:
    def __init__(self,fluorescence,temperatures,name,cv1,cv2,ph,dphdt,isControl):
        self.fluorescence = fluorescence
        self.temperatures = temperatures
        self.contents = Contents(cv1,cv2,ph,dphdt,isControl)
        self.name = name

        self.normalisationFactor = None
        self.wellMax = None
        self.wellMin = None
        self.wellNormalisedMax = None 
        self.wellNormalisedMin = None
        self.tm = None

        self.isMonotonic = False
        self.isComplex =  False
        self.isOutlier = False
        self.isInTheNoise = False
        self.isSaturated = False
        self.isDiscarded = False

        return
    
    def getMinAndMax():
        minimum = fluorescence[0]
        maximum = fluorescence[0]
        for value in fluorescence[1:]:
            if value > maximum:
                maximum = value
            if value < minimum:
                minimum = value
        return (minimum,maximum)
    
    def normalise():
        count = 0
        for height in self.fluorescence:
            count += height*stepSize
        self.fluorescence = [x / count for x in self.fluorescence]
        #used to calculate the monotenicity threshold
        self.normalisationFactor = count
        return
    
    def computeSaturation():
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
        return
    
    def computeInTheNoise(noiseThreshold):
        return
    
    def computeTm():
        return
    
    def computeComplexity():
        return
    
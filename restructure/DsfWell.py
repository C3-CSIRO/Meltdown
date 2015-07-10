# -*- coding: utf-8 -*-

import math
import pandas as pd
import numpy as np
import Tkinter, tkMessageBox

#the max amount the flat saturated curves can fluctuate within the flat section
SATURATION_FLUCTUATION_THRESHOLD = 0.005
#how long a flat section on curve can be before it is considered saturated
LENGTH_OF_FLAT_CONSIDERED_SATURATED = 10
#the fraction of the total curve, at the end, that will not be iterated over for Tm calculation, as last section is unreliable
FRACTION_OF_CURVE_NOT_CHECKED_FOR_TM = 0.125
#number of times a monotonic curve can break its monotonicity but still be considered monotonic
MONOTONIC_CONTRADICTION_LIMIT = 5
#when finding sign changes in the derivative series, this is the forgiving threshold from 0
SIGN_CHANGE_THRESH = 0.000001
#the largest difference between the calculated Tm, and the mid point of the highest and lowest points on the curve before
#curve is considered complex
MAX_DIFFERENCE_BETWEEN_TMS_BEFORE_COMPLEX = 5

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
        if not self.isDiscarded:
            currentValue = self.fluorescence[0]
            maxInd = 0
            while currentValue != self.wellNormalisedMax:
                maxInd+=1
                currentValue = self.fluorescence[maxInd]

            diff = self.wellNormalisedMax - self.wellNormalisedMin

            # A boundry defining how much the points can fluctuate and still be considered flat
            lowFlatBoundry = self.wellNormalisedMax - SATURATION_FLUCTUATION_THRESHOLD*diff

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
            if count >= LENGTH_OF_FLAT_CONSIDERED_SATURATED:
                self.isSaturated = True
                self.isDiscarded = True

        return
    
    def computeMonotonicity(self, plateMonotonicThreshold):
        #calculate the well's individual monotonic threshold from the plate's one
        self.wellMonotonicThreshold = plateMonotonicThreshold / self.normalisationFactor
        
        #no need to calculate if curve is monotonic, if it is already tagged as discarded
        if self.isDiscarded:
            return
        
        #initially assume curve is decreasing monotonic
        decreasingMonotonic = True
        contradictions = 0
        prev=self.fluorescence[0]
        for point in self.fluorescence[1:]:
            if point > prev+self.wellMonotonicThreshold:
                #found contradiction to monotonicity, increase contradication counter
                contradictions+=1
                
            elif point < prev+self.wellMonotonicThreshold:
                #lower contradiction counter if points are decreasing, but don't take counter below 0
                if contradictions!=0:
                    contradictions-=1
                    
            #when point=previous do nothing with the contradiction counter
            
            #if enough consecutive contradictions are hit, clasify the curve as decreasing monotonic
            if contradictions == MONOTONIC_CONTRADICTION_LIMIT:
                decreasingMonotonic=False
                break
            prev=point
            
        #if curve is found to be monotonic, change appropriate variables
        if decreasingMonotonic:
            self.isMonotonic = True
            self.isDiscarded = True
        else:
            self.isMonotonic = False
        return
    
    def computeInTheNoise(self, noiseThreshold):
        #no need to calculate if curve already discarded
        if self.isDiscarded:
            return
        #curve is in the noise, and should be discarded, if its monotonic threshold is greater than the noise threshold given
        if self.wellMonotonicThreshold > noiseThreshold:
            self.isInTheNoise = True
            self.isDiscarded = True
        else:
            self.isInTheNoise = False
        return
    
    def computeTm(self):
        #if well is monotonic, saturated, in the noise, or an outlier, then don't try to find its Tm
        if self.isDiscarded:
            return
        x = self.temperatures
        y = self.fluorescence
        #get the derivative series (each point is the slope between successive points in the normalised curve)
        xdiff = np.diff(x)
        dydx = -np.diff(y)/xdiff
        #the derivative series has one less index since there is one fewer differences than points
        seriesDeriv = pd.Series(dydx, x[:-1])
        
        #now that we have the derivative series, we can find the Tm
        lowestPoint = 0
        lowestPointIndex = None
        #since the end of the melt curves is often very unpredictable, we only search for a Tm up to a point
        ignoreIndex = -int(len(seriesDeriv.index)*FRACTION_OF_CURVE_NOT_CHECKED_FOR_TM)
        #find the lowest point in the derivative series
        for ind in seriesDeriv.index[:ignoreIndex]:
            if seriesDeriv[ind]<lowestPoint:
                lowestPoint = seriesDeriv[ind]
                lowestPointIndex = ind
        
        #if lowest point was found to be at the start or end of the checked derivative series, then no curve fit is required
        if lowestPointIndex == seriesDeriv.index[0] or lowestPointIndex == seriesDeriv.index[ignoreIndex-1]:
            self.tm = lowestPointIndex
            return
        #if no lowest point could be found, leave the Tm as none
        if lowestPointIndex == None:
            self.tm = None
            #force curve to be complex if no Tm can be found, and it has not been discarded
            self.isComplex = True
            return
        
        #get the index's either side of the lowest point
        leftIndex = list(seriesDeriv.index)[list(seriesDeriv.index).index(lowestPointIndex) - 1]
        rightIndex = list(seriesDeriv.index)[list(seriesDeriv.index).index(lowestPointIndex) + 1]
        
        #matrices used to fit a parabola to the 3 points
        Y=[seriesDeriv[leftIndex],
           seriesDeriv[lowestPointIndex],
           seriesDeriv[rightIndex]]
           
        A=[[leftIndex**2,   leftIndex,   1],
           [lowestPointIndex**2, lowestPointIndex, 1],
           [rightIndex**2,  rightIndex,  1]]
           
        #solves for b, in the form Y=Ab
        (a,b,c) = np.linalg.solve(A,Y)
        
        #initialise tm to left most point of relevant curve
        tm=seriesDeriv[leftIndex]
        tmFluorescence=0
        #set tm to the lowest point on the fitted parabola rounded to nearest 0.01
        for x in np.arange(leftIndex,rightIndex,0.01):
            point = (a*(x**2) + b*x + c)
            if point < tmFluorescence:
                tmFluorescence = point
                tm = x
        #save the tm and exit
        self.tm = tm
        return
    
    def computeComplexity(self):
        #only calculate if curve is not discarded, and not already marked as complex
        if self.isDiscarded or self.isComplex:
            return
        x = self.temperatures
        y = self.fluorescence
        #get the derivative series (each point is the slope between successive points in the normalised curve)
        xdiff = np.diff(x)
        dydx = -np.diff(y)/xdiff
        #the derivative series has one less index since there is one fewer differences than points
        seriesDeriv = pd.Series(dydx, x[:-1])
        
        #checks from a derivative sign change between the highest and lowest points on the curve
        lowestPoint = 1
        lowestIndex = None
        highestPoint = 0
        highestIndex = None
        previous = None
        for i, value in enumerate(self.fluorescence[:-1]):
            if value > highestPoint:
                highestPoint = value
                highestIndex = i
        if highestIndex == 0 :
            highestPoint = 0
            highestIndex = None
            for i, value in enumerate(self.fluorescence[:-1]):
                if value<lowestPoint:
                    lowestPoint = value
                    lowestIndex = i
            for i, value in enumerate(self.fluorescence[:-1]):
                if i < lowestIndex:
                    continue
                if value > highestPoint:
                    highestPoint = value
                    highestIndex = i
        else:
            for i, value in enumerate(self.fluorescence[:-1]):
                if i > highestIndex:
                    break
                if value<lowestPoint:
                    lowestPoint = value
                    lowestIndex = i
        signChange = False
        for ind in seriesDeriv.index[lowestIndex+1:highestIndex]:
        
            if previous:
                if seriesDeriv[ind] + SIGN_CHANGE_THRESH < 0 and previous - SIGN_CHANGE_THRESH > 0:
                    signChange = True
                if seriesDeriv[ind] - SIGN_CHANGE_THRESH > 0 and previous + SIGN_CHANGE_THRESH < 0:
                    signChange = True
            previous = seriesDeriv[ind]
        #if we have a sign change, curve is complex
        if signChange:
            self.isComplex = True
            return
        
        #other check for complex curve, uses another estimate for Tm
        averagePoint = (lowestPoint +highestPoint) / 2
        i = lowestIndex
        while self.fluorescence[i]<averagePoint:
            i += 1;
        #if difference between previously calculated Tm, and new estimate is too large the curve is considered complex
        if math.fabs(self.temperatures[i] -self.tm) > MAX_DIFFERENCE_BETWEEN_TMS_BEFORE_COMPLEX:
            self.complex=True
        return
        
    def setAsOutlier(self):
        #since outlier computation is done outlside of DsfWell class, isDiscarded shouldn't need to be set elsewhere
        self.isOutlier = True
        self.isDiscarded = True
        return



def main():
    root = Tkinter.Tk()
    root.withdraw()
    tkMessageBox.showwarning("Inncorrect Usage", "Please run the 'RunMeltdown.bat' file from the same directory")
    return
    
    
if __name__ == "__main__":
    main()











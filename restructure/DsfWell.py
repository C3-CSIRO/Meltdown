# -*- coding: utf-8 -*-


class DsfWell:
    def __init__(self,fluorescence,temperatures,cv1,cv2,ph,dphdt,isControl):
        self.fluorescence = fluorescence
        self.temperatures = temperatures
        self.contents = Contents(cv1,cv2,ph,dphdt,isControl)
        return
    
    def getMinAndMax():
        return
    
    def normalise():
        return
    
    def computeSaturation():
        return
    
    def computeMonotonicity(plateMonotonicThreshold):
        return
    
    def computeInTheNoise(noiseThreshold):
        return
    
    def computeTm():
        return
    
    def computeComplexity():
        return
    
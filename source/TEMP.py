# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import scipy.stats

def lowestPoints():
    
    x = [0,1,2,3,4,5,6,7,8,9]
    y = [4,3,1,2,5,7,6,5,3,4]
    plt.plot(x,y)
    lowestPoint = 10
    lowestPointIndex = None
    lowestPoint2 = 10
    lowestPointIndex2 = None
    
    upAfterMin = False
    maxAfterMin = 0
    
    for i, point in enumerate(y):
        if lowestPointIndex != None and point > lowestPoint:
            upAfterMin = True
        
        if point <  lowestPoint:
            
            lowestPoint = point
            lowestPointIndex = i
        
        
        
        
    print "lowest point: ", lowestPoint, "at index: ", lowestPointIndex
    print "lowest point 2: ", lowestPoint2, "at index: ", lowestPointIndex2
    
    return


def lowestPoints2():
    x = [0,1,2,3,4,5,6,7,8,9]
    y = [4,3,1,2,5,7,6,5,3,4]
    plt.plot(x,y)
    
    mins = []
    lowest = 10
    indexOfLowest = None
    
    for i,point in zip(x,y):
        if point < lowest:
            lowest = point
            indexOfLowest = i
    
    prev = y[indexOfLowest]
    indexOfRightPeak = indexOfLowest
    for i,point in enumerate(y[indexOfLowest:]):
        if prev > point:
            break
        
        indexOfRightPeak = i+indexOfLowest
        prev = point
        
    prev = y[indexOfLowest]
    indexOfLeftPeak = indexOfLowest
    for i,point in enumerate(reversed(y[:indexOfLowest+1])):
        if prev > point:
            break
        
        indexOfLeftPeak = indexOfLowest-i
        prev = point
    
    
    print indexOfLeftPeak, indexOfLowest, indexOfRightPeak, lowest
    plt.show()
    
    return

def getMinAndPeaks(list_, startIndex, endIndex):
    
    
    return




def tmp():    
    a = [0.482490272,
    0.394849785,
    0.319526627,
    0.282608696,
    0.266055046,
    0.285714286,
    0.245742092,
    0.201729107,
    0.196825397,
    0.197058824,
    0.127586207,
    0.095419847]
    
    b = [0.40077821,
    0.354792561,
    0.297548605,
    0.251207729,
    0.266055046,
    0.232804233,
    0.253041363,
    0.181556196,
    0.136507937,
    0.179411765,
    0.137931034,
    0.091603053]
    
    ans = scipy.stats.ttest_ind(a, b, axis=0, equal_var=False)
    print ans
    return

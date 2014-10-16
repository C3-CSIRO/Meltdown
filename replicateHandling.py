"""
Authers: 
Nicholas Rosa and Marko Ristic

Date: 
*note github date

Synopsis: 
This .py file stores all the functions used in dealing with the duplicate 
wells from a given DSF experiment.

The aim is for these functions to be used to return a pandas dataframe of 
meltcurve data which contains no duplicates. the duplicates will be processed 
and averaged to produce an average melt curve, while discarding cases which are
calculated to be innacurate

simiarity is based on the sum of squares of differences at every point from a 
calculated 'ideal' curve. We have chosen the idea curve to be a normalised (from 0 to 1) average 
curve of Lysozme with the Tm moved to be at zero. Every consecutive curve being compared 
to this ideal, is also normalised and its Tm also moved to zero.

Tm of any curve was calculated as the lowest point on the derivative melt data.

"""

import Tkinter, tkMessageBox
import math
import numpy as np
from itertools import combinations


def meanSd(listOfNumbers):
    """
    Finds the mean and standard deviation of a group of numbers
    
    Input: Takes in a list of numbers
    
    Output: Returns a tuple that contains the mean and the standard deviation
    """
    etotal=0
    e2total=0
    N=len(listOfNumbers)
    for item in listOfNumbers:
        if item !=None:
            etotal+=item
            e2total+=math.pow(item,2)
        else:
            N -= 1
    if N == 0:
        return (None,None)
    moment1=etotal/N
    moment2=e2total/N
    variance = moment2 - math.pow(moment1,2)
    return (moment1, np.sqrt(variance))

#TODO get rid of overlap, decide on whether to divide by len(index)
def sqrdiff(series1,series2):
    indexoverlap=[]
    for x in series1.index:
        if x in series2.index:
            indexoverlap.append(x)
    sqrsum=0
    for ind in indexoverlap:
        sqrsum+=math.pow(series1[ind]-series2[ind],2)
    
    return sqrsum / len(indexoverlap)


def discardBad(wellNamesList,matrix,thresh):
    """using the dist matrix, groups are made of those proteins which are within the threshold
    of each other, and the most tightly packed grouped is returned. If all are unrelated, an 
    empty group is returned"""
    #keys are wells, values are a list of their group members
    wells = {}
    for name in wellNamesList:
        #the groups in which the wells are placed are the indexes of the wells in the input list
        wells[name]=[name]
    for i in range(len(matrix)):
        #loop over all the matrix except the diagonal (always 0)
        for j in [x for x in range(len(matrix[i])) if x!=i]:
            if matrix[i][j] > thresh:
                #remove edges that are greater than the threshold
                matrix[i][j] = -1
            else:
                #add the two wells within threshold to each others groupmembers' groups (including own)
                for wlj in wells[wellNamesList[j]]:
                    if wellNamesList[i] not in wells[wlj]:
                        wells[wlj].append(wellNamesList[i])
                for wli in wells[wellNamesList[i]]:
                    if wellNamesList[j] not in wells[wli]:
                        wells[wli].append(wellNamesList[j])
    #finding the beet group
    maxi=-1
    for key in wells.keys():
        if len(wells[key]) > maxi:
            maxi = len(wells[key])
            maxkey = key   #used as previous equal best if there is an equal max group
            choose = key   #used as the key to the group to be chosen
        elif len(wells[key]) == maxi:
            #if the length is the same because key and maxkey are in the same group, simply skip to next loop iteration
            if key in wells[maxkey]:
                continue
            #equal largest group, choose the one with the smallest 'largest difference' amoong the pairs
            largestValKey = -1
            largestValMaxkey = -1
            for pair in combinations(wells[key],2):
                #finds largest pair comparison value in current group
                if matrix[wellNamesList.index(pair[0])][wellNamesList.index(pair[1])] > largestValKey:
                    largestValKey = matrix[wellNamesList.index(pair[0])][wellNamesList.index(pair[1])]
                    choose = key
            for pair2 in combinations(wells[maxkey],2):
                #find largest pair comparison value in previous equal size group
                if matrix[wellNamesList.index(pair2[0])][wellNamesList.index(pair2[1])] > largestValMaxkey:
                    largestValMaxkey = matrix[wellNamesList.index(pair2[0])][wellNamesList.index(pair2[1])]
                    choose = maxkey
            #choose group with the smaller largest value
            if largestValKey > largestValMaxkey:
                choose = key
            elif largestValKey < largestValMaxkey:
                choose = maxkey
            else:
                #two values are equal, both groups are discarded
                #either both groups are length 1, or there's no way to choose the better one
                choose=""
    #appends the group members of the chosen group to be returned
    keep=[]
    if choose != "":
        for well in wells[choose]:
            keep.append(well)
        return keep
    else:
        return []


def main():
    root = Tkinter.Tk()
    root.withdraw()
    tkMessageBox.showwarning("Inncorrect Usage", "Please run the 'meltdown.py' file from the same directory")
    return
    
    
if __name__ == "__main__":
    main()

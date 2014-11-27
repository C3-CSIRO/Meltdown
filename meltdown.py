# -*- coding: UTF-8 -*-
#!/usr/bin/env python
#
#    *                                                       
#  (  `            (       )    (                            
#  )\))(      (    )\   ( /(    )\ )         (  (            
# ((_)()\    ))\  ((_)  )\())  (()/(    (    )\))(     (     
# (_()((_)  /((_)  _   (_))/    ((_))   )\  ((_)()\    )\ )  
# |  \/  | (_))   | |  | |_     _| |   ((_) _(()((_)  _(_/(  
# | |\/| | / -_)  | |  |  _|  / _` |  / _ \ \ V  V / | ' \)) 
# |_|  |_| \___|  |_|   \__|  \__,_|  \___/  \_/\_/  |_||_| 
#
#  meltdown.py
#
#  Authors:  Nicholas Rosa, Marko Ristic, Shane A. Seabrook, David Lovell, Del Lucent and Janet Newman
#
#  Website: TBD.github.com
#
#  License: XXX
#
# This software will generate automated reports from high-throughput thermofluor protein
# stability measurements.  Although originally designed to provide a report for the Buffer Screen 9
# experiment at the Collaborative Crystalisation Centre (C3), this can be extended to analyse any
# similar experiment in 96 well format.  Given two input xls files (one for fluorescence in
# each well as a function of temperature, and one describing the buffer formulation of each well in
# a standardized format) the following tasks are performed:
#
#  + data remediation including outlier removal of replicants
#  + Tm calculation and complex shape flag for curves
#  + report generation summarizing protein thermal stability as a function of buffer formulation,
#    Tm graph, formulation plots, control well checking, etc.
# 
# Additionally, given a suitable set of training data, the parameters used for outlier detection,
# Tm calculation, and curve classification can be recomputed.
#

import math
import pandas
import xlrd
import Tkinter, tkFileDialog, tkMessageBox
import cStringIO
import os
import sys, traceback
import zipfile as zf
import matplotlib.pyplot as plt
import numpy as np
from itertools import combinations
#reportlab needs to be installed separetly by anaconda, so a messagebox pops up alerting the user
#if there is trouble importing it
try:
    from reportlab.pdfgen import canvas
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.units import cm
    from reportlab.lib.utils import ImageReader
except:
    root = Tkinter.Tk()
    root.withdraw()
    tkMessageBox.showwarning("ReportLab not found", "You must use Anaconda to install reportlab before Meltdown can be run")

#files from directory
import replicateHandling as rh

##====GLOBAL VARIABLES====##
#discard bad threshold, mean difference between any two normalised lysozyme curves
#TODO LOG
#SIMILARITY_THRESHOLD = 8.10785059623e-05
#lysozyme tm over all of our files: (mean,std_dev)
LYSOZYME_TM_THRESHOLD = (70.87202380952381, 0.73394932964132509)
#Monotenicity threshold forgive value on non-normalised curves
DEFAULT_MONO_THRESH = 10
#the different colours of the saltconcentrations, in order of appearance
COLOURS = ["blue","darkorange","green","red","cyan","magenta"]
#NOTE first 4 are from custom input, second 4 are the oldnames  from default pcrd export
#individual order is also important:
CONTROL_WELL_NAMES = ["lysozyme","no dye","protein as supplied","no protein"]

##====DEBUGGING====##
#set this to false if you do not wish for the exported data files to be deleted after being analysed
deleteInputFiles = False



class DSFAnalysis:
    """
    Class used for performing the analysis of a DSF plate, methods will
    remove outliers, find mean curves, find tm's and check for monotenicity
    along with produce the final report
    """
    
    def __init__(self):
        """
        Initialises an empty list of deleted curves which is filled with unused curves
        when removing outliers from replicates
        """
        self.delCurves = []
        return
    
    def loadMeltCurves(self, fluorescenceXLS, labelsXLS):
        """
        Loads the two XLS or XLSX files and takes from them the information
        needed to start the analysis.
        """
        #populates the relevant instance variables
        self.name = fluorescenceXLS
        self.plate = DSFPlate(fluorescenceXLS, labelsXLS)
        self.wells = self.plate.wells
        self.originalPlate = DSFPlate(fluorescenceXLS, labelsXLS)
        self.removeOutliers()
        self.removeInsignificant()
        self.findMeanCurves()
        return
    
    def findMeanCurves(self):
        """
        Of the two plates created, we turn one (self.plate) into its own mean plate.
        Every group of retained replicates (remove outliers has already been called) is
        averaged out, and place in a well of the appropriate mean well name.
        well group A4, A5, A6 in a 3 replicate well would have mean well name A2
        """
        meanWells = {}
        visited = []
        numcount = 1
        letcount = 0
        ordcount = ord("A")
        names = sorted(self.plate.names,key=lambda well: int(well[1:]))
        for well in sorted(names, key=lambda well: well[0]):
            if well not in visited:
                reps = []
                reps += self.originalPlate.repDict[well]

                i = 0
                while i < len(reps):
                    visited.append(reps[i])
                    if reps[i] in self.delCurves:
                        del reps[i]
                        i -=1    
                    i+=1
                # if the first letter of the well name is not the same as ordcount
                if ord(well[0]) != ordcount:
                    ordcount = ord(well[0])
                    letcount += 1
                    numcount = 1
                if len(reps) > 0:
                    meanWells[chr(65+letcount) + str(numcount)] = (meanCurve([self.wells[well].fluorescence for well in reps]),well)
                else:
                    meanWells[chr(65+letcount) + str(numcount)] = (None,well)
                if well in self.plate.lysozyme:
                    self.plate.lysozyme = [chr(65+letcount) + str(numcount)]
                if well in self.plate.noDye:
                    self.plate.noDye = [chr(65+letcount) + str(numcount)]
                if well in self.plate.proteinAsSupplied:
                    self.plate.proteinAsSupplied = [chr(65+letcount) + str(numcount)]
                if well in self.plate.noProtein:
                    self.plate.noProtein = [chr(65+letcount) + str(numcount)]

                
                #creates a dictionary allowing you to find which wells created which mean well
                #e.g. 'A2': ["A4","A5","A6"]
                self.plate.meandict[chr(65+letcount) + str(numcount)] = self.originalPlate.repDict[well]
                
                #numcount increased for naming mean wells
                numcount += 1
                
        remove = []

        for well in self.plate.names:
            if well in meanWells.keys():
                if len(self.wells) == 0:
                    return
                if well not in self.plate.names:
                    new = self.plate.names[0]
                    self.wells[well] = self.wells[new]
                if meanWells[well][0]:
                    self.wells[well].fluorescence = meanWells[well][0] 
                else:
                    self.wells[well].fluorescence = None
                self.wells[well].contents = self.originalPlate.wells[meanWells[well][1]].contents
            elif well in self.plate.names:
                del self.wells[well]
                remove.append(well)
            else:
                remove.append(well)
        for well in remove:
            self.plate.names.remove(well)
        return


    
    def removeOutliers(self):
        """
        Removes curves from analysis by adding them to self.delCurves. The removal of
        a curve is decided by it similarity to its replicates. This is further described 
        in replicateHnadling.py
        """
        #With the DSFPlate object, we can just use self.wells.pop() to remove outliers
        visited = []
        discard = []
        for well in self.wells:
            if well not in visited:
                reps = []
                reps += self.originalPlate.repDict[well]
                pairs = combinations(reps,2)
                distMatrix = [[0 for x in range(len(reps))] for y in range(len(reps))]
                for pair in pairs:
                    dist = sqrDiffWellFluoro(self.wells[pair[0]].fluorescence,self.wells[pair[1]].fluorescence)
                    distMatrix[reps.index(pair[0])][reps.index(pair[1])] = dist
                    distMatrix[reps.index(pair[1])][reps.index(pair[0])] = dist
                keep = rh.discardBad(reps,distMatrix,SIMILARITY_THRESHOLD)
                for rep in reps:
                    visited.append(rep)
                    if rep not in keep:
                        discard.append(rep)
        for well in discard:
            self.wells[well].fluorescence = None
            self.delCurves.append(well)
        return

    def removeInsignificant(self):
        """
        Curves which were normalised with too large of a factor, are considered too 
        noise for any  sort of useful analysis. These carves are found and added to
        self.delCurves (removed from analysis)
        """
        thresholdm, i = rh.meanSd([self.originalPlate.wells[x].monoThresh for x in self.plate.noProtein])
        for well in self.wells:
            if well not in self.plate.lysozyme and well not in self.plate.noProtein and well not in self.plate.noDye:
                if self.wells[well].monoThresh > thresholdm/1.15:
                    #self.wells[well].fluorescence = None
                    self.delCurves.append(well)

            if self.wells[well].fluorescence:
                x = [x for x in self.wells[well].temperatures]
                y = [y for y in self.wells[well].fluorescence]
                xdiff = np.diff(x)
                dydx = -np.diff(y)/xdiff

                #the derivative series, has one less index since there is one fewer differences than points
                seriesDeriv = pandas.Series(dydx, x[:-1])
                mini = 0
                for ind in seriesDeriv.index[:-20]:
                    if seriesDeriv[ind]<mini:
                        mini = seriesDeriv[ind]

                if mini > -0.00005:
                    if well not in self.delCurves and well not in self.plate.lysozyme and well not in self.plate.noProtein and well not in self.plate.noDye:
                        self.delCurves.append(well)
        return
    
    def analyseCurves(self):
        """
        Calculate the Tms of all the curves in the plate, and also
        check that the controls on the plate are within the expected range
        """
        self.computeTms()
        self.returnControlCheck()
        return
        
    def returnControlCheck(self):
        """
        Returns a dictionary of the different controls, with values True or False
        depending on if the control passed or not. This is used when creating the report
        """
        #initialises all controls to fail
        result = {CONTROL_WELL_NAMES[0]:False,  #lysozyme Tm check
                  CONTROL_WELL_NAMES[1]:False,  #no dye similarity check
                  CONTROL_WELL_NAMES[3]:False}  #no protein check
               
        #only test the control if the control is present in the plate
        if len(self.plate.lysozyme)>0:
            #lysozyme Tm check
            if self.plate.wells[self.plate.lysozyme[0]].Tm > LYSOZYME_TM_THRESHOLD[0] - 2*LYSOZYME_TM_THRESHOLD[1] and\
               self.plate.wells[self.plate.lysozyme[0]].Tm < LYSOZYME_TM_THRESHOLD[0] + 2*LYSOZYME_TM_THRESHOLD[1]:
                    result[CONTROL_WELL_NAMES[0]] = True
                
        #only test the control if the control is present in the plate
        if len(self.plate.noDye)>0:
            #get the curves to compare as series
            noDyeControl = pandas.Series.from_csv("data/noDyeControl.csv")
            noDyeReal = pandas.Series(self.plate.wells[self.plate.noDye[0]].fluorescence, self.plate.wells[self.plate.noDye[0]].temperatures)
            #if the curves are within required distance from one another, the control is passed
            if rh.sqrdiff(noDyeReal, noDyeControl) < SIMILARITY_THRESHOLD:
                result[CONTROL_WELL_NAMES[1]] = True
                
        #only test the control if the control is present in the plate
        if len(self.plate.noProtein)>0:
            #get the curves to compare as series
            noProteinControl = pandas.Series.from_csv("data/noProteinControl.csv")
            noProteinReal = pandas.Series(self.plate.wells[self.plate.noProtein[0]].fluorescence, self.plate.wells[self.plate.noProtein[0]].temperatures)
            #if the curves are within required distance from one another, the control is passed
            if rh.sqrdiff(noProteinReal, noProteinControl) < SIMILARITY_THRESHOLD:
                result[CONTROL_WELL_NAMES[3]] = True
            
        return result
    
    def generateReport(self, outFile):
        """
        Generate a pdf report with the appropriate summary information.
        
        Input: Outfile is the name of the file to which the report is saved
        """
        pdf = canvas.Canvas(outFile,pagesize=A4) 
        pdf.setFont("Helvetica-Bold",16)
        pdf.drawString(cm,28*cm,"MELTDOWN")
        pdf.setFont("Helvetica",16)
        pdf.drawString(6*cm,28*cm,"Melt Curve Analysis")
        for i in range(1,len(self.name)):
            if self.name[-i] == '/' or self.name[-i] == '\\':
                if -i+40 < 0:
                    nameMatch = self.name[-i+1:-i+40] +"..."
                else:
                    nameMatch = self.name[-i+1:]
                break
        if nameMatch != None:
            pdf.drawString(cm,27*cm, nameMatch)
        else:
            pdf.drawString(cm,27*cm, "..."+self.name[-70:])
        pdf.setFont("Helvetica-Bold",12)
        pdf.drawImage("data/CSIRO_Grad_RGB_hr.jpg",17*cm,25.5*cm,3.5*cm,3.5*cm)
        best = ""
        maxi = 0
        mini = 100
        for well in self.plate.names:
            if self.wells[well].contents.name.lower() not in CONTROL_WELL_NAMES and self.wells[well].Tm != None and self.wells[well].mono == False and self.wells[well].Tm > maxi:
                maxi = self.wells[well].Tm
                best = well
        if best != "":
            if self.wells[best].TmError != None:
                pdf.drawString(3*cm,3.6*cm,"Highest Tm = " + str(round(self.wells[best].Tm,2)) + " +/- " + str(round(self.wells[best].TmError,2)))
            else:
                pdf.drawString(3*cm,3.6*cm,"Highest Tm = " + str(round(self.wells[best].Tm,2)))
            pdf.drawString(3*cm,3*cm,"("+self.wells[best].contents.name+" / "+self.wells[best].contents.salt+")")

        pdf.setFont("Helvetica",12)
        fig1 = plt.figure(num=1,figsize=(10,8))
        maxi = 0

        conditions = [self.wells[x].contents.name+" ("+str(self.wells[x].contents.pH)+")" for x in self.plate.names if self.wells[x].contents.name.lower() not in CONTROL_WELL_NAMES]
        labels = []
        for item in conditions:
            if item not in labels and item[:-7] not in labels:
                if "None" in item:
                    labels.append(item[:-7])
                else:
                    labels.append(item)
        labels = Contents.name
        tmHandles = []
        # plotting the lines for each salt concentration
        for i, saltConcentration in enumerate(Contents.salt):
            curves = []
            for well in self.plate.names:
                if self.wells[well].contents.salt == saltConcentration:
                    curves.append(well)
            tms = []#[self.wells[x].Tm for x in curves]
            badTms = []

            for well in curves:
                if self.wells[well].Tm != None and self.wells[well].TmError == None or self.wells[well].complex == True:
                    tms.append(None)
                    badTms.append(self.wells[well].Tm)
                else:
                    tms.append(self.wells[well].Tm)
                    badTms.append(None)

            for val in tms:
                if val != None:
                    if val > maxi:
                        maxi = val
                    if val < mini:
                        mini = val
            for val in badTms:
                if val != None:
                    if val > maxi:
                        maxi = val
                    if val < mini:
                        mini = val

            handle, = plt.plot([x for x in range(len(labels))],tms,color=COLOURS[i],marker="o",linestyle="None")
            plt.plot([x for x in range(len(labels))],badTms,color=COLOURS[i],marker="d",linestyle="None")
            if badTms:
                crosses = True
            tmHandles.append(handle)


        originalProteinMeanSd = rh.meanSd([self.originalPlate.wells[x].Tm for x in self.originalPlate.proteinAsSupplied])
        if originalProteinMeanSd[0]!= None:
            
            if(originalProteinMeanSd[0]-mini > maxi - originalProteinMeanSd[0]):
                plt.axis([-1,len(labels),mini-5,2*originalProteinMeanSd[0]-mini +10])
            else:
                plt.axis([-1,len(labels),2*originalProteinMeanSd[0] - maxi-5,maxi + 10])
        else:
            plt.axis([-1,len(labels),mini-5,maxi + 5])
        plt.gcf().subplots_adjust(bottom=0.35)
        plt.ylabel('Tm')

        plt.axhline(originalProteinMeanSd[0],0,1,linestyle="--",color="red")
        plt.xticks([x for x in range(len(labels))],labels,rotation="vertical")
        plt.legend(tmHandles,Contents.salt)

        imgdata = cStringIO.StringIO()
        fig1.savefig(imgdata, format='png',dpi=180)
        imgdata.seek(0)  # rewind the data
        Image = ImageReader(imgdata)
        pdf.drawImage(Image, cm, 4*cm, 16*cm, 11*cm)  
        plt.close()

        pdf.setFont("Helvetica",10)
        if crosses:
            pdf.drawString(7.9*cm, 14.2*cm, "Tms drawn in diamonds may be unreliable")
        crosses = False
        
        if originalProteinMeanSd[0]!=None:
            pdf.drawString(15.5*cm,10.4*cm,"Protein as supplied")

        pdf.setFillColor("black")


        controlChecks = self.returnControlCheck()

        #TODO Change position
        # Lysozyme Tm Control Check
        if controlChecks["lysozyme"] == False:
            pdf.setFillColor("red")
            pdf.drawString(3*cm,9.5*cm,"Lysozyme Control: Failed")
            pdf.setFillColor("black")

        # No dye control check 
        if controlChecks["no dye"] == False:
            pdf.setFillColor("red")
            pdf.drawString(3*cm,9*cm,"No Dye Control: Failed")
            pdf.setFillColor("black")
            
        # No Protein control check
        if controlChecks["no protein"] == False:
            pdf.setFillColor("red")
            pdf.drawString(3*cm,8.5*cm,"Dye Only Control: Failed")
            pdf.setFillColor("black")

        fig2 = plt.figure(num=1,figsize=(5,4))
        for well in self.originalPlate.proteinAsSupplied:
            if well in self.delCurves:
                plt.plot(self.originalPlate.wells[well].temperatures,self.originalPlate.wells[well].fluorescence\
                , 'grey')
            elif self.originalPlate.wells[well].mono == False:
                plt.plot(self.originalPlate.wells[well].temperatures,self.originalPlate.wells[well].fluorescence\
                , 'g')
            else:
               plt.plot(self.originalPlate.wells[well].temperatures,self.originalPlate.wells[well].fluorescence\
                , 'g',linestyle="--")

        plt.gca().axes.get_yaxis().set_visible(False)
        imgdata = cStringIO.StringIO()
        fig2.savefig(imgdata, format='png',dpi=140)
        imgdata.seek(0)  # rewind the data
        Image = ImageReader(imgdata)
        pdf.drawImage(Image, 0, 18*cm, 8*cm, 6*cm)
        plt.close()

        # Tm of the protein in the oringinal formulation
        try:
            originalProteinMeanSd = rh.meanSd([self.originalPlate.wells[x].Tm for x in self.originalPlate.proteinAsSupplied])
        except:
            originalProteinMeanSd = (None, None)
        if originalProteinMeanSd[0]!= None:
            pdf.drawString(cm,17.5*cm, "Protein as supplied: Tm = " +str(round(originalProteinMeanSd[0],2))+"(+/-"+str(round(originalProteinMeanSd[1],2))+")")
        else:
            pdf.drawString(cm,17.5*cm, "Protein as supplied: Tm = N/A")

        pdf.setFont("Helvetica-Bold",13)
        pdf.drawString(8*cm,22.75*cm,"Full interpretation of the results requires you to look ")
        pdf.drawString(8*cm,22.25*cm,"at the individual melt curves.")

        avTmError = 0
        count = 0
        tmCount = 0
        for well in self.originalPlate.names:
            if well not in self.originalPlate.lysozyme + self.originalPlate.noDye + self.originalPlate.noProtein + self.originalPlate.proteinAsSupplied:
                if self.originalPlate.wells[well].Tm != None and self.originalPlate.wells[well].mono == False and well not in self.delCurves:
                    tmCount += 1
                count += 1
        tmCount = int(round(tmCount/float(count),2)*100)
        pdf.drawString(8*cm,20.5*cm,str(tmCount)+"%")
        pdf.setFont("Helvetica",13)
        pdf.drawString(9.25*cm,20.5*cm,"of curves were used in Tm estimations")

        count = 0
        for well in self.plate.names:
            if self.wells[well].TmError != None:
                avTmError += self.wells[well].TmError
                count += 1
        pdf.drawString(8*cm,19.5*cm,"Average estimation of error is")
        pdf.setFont("Helvetica-Bold",13)
        avTmError = round(avTmError/float(count),1)
        pdf.drawString(14.1*cm,19.5*cm,str(avTmError)+" C")

        pdf.setFont("Helvetica",13)
        pdf.drawString(8*cm,18.5*cm, "Protein as supplied is")
        pdf.setFont("Helvetica-Bold",13)
        wellBehaved = True

        for well in self.originalPlate.proteinAsSupplied:
            if self.originalPlate.wells[well].Tm == None or self.originalPlate.wells[well].mono == True or well in self.delCurves:
                wellBehaved = False
        if originalProteinMeanSd[1] > 1.5:
            wellBehaved = False
        if wellBehaved:
            pdf.drawString(12.5*cm,18.5*cm,"well behaved")
        else:
            pdf.drawString(12.5*cm,18.5*cm,"not well behaved")

        if wellBehaved == False or avTmError >=1.5 or tmCount <= 50:
            pdf.drawString(8*cm,21.5*cm,"The summary graph appears to be unreliable")

        pdf.rect(7.75*cm,18.05*cm,12*cm,5.4*cm)

        # Saving the first page and moving onto the second page
        pdf.showPage()
        pdf.setFont("Helvetica",12)    
        
        pdf.setFont("Helvetica",9)
        pdf.drawString(cm, 1.3*cm,"Curves drawn with dashed lines are monotonic and excluded from Tm calculations")
        pdf.drawString(cm, 0.9*cm,"Curves with complex melt transitions are marked (^) and are drawn with a dotted line")
        pdf.drawString(cm, 0.5*cm,"Curves coloured grey are outliers, and are excluded from Tm calculations")
        pdf.setFont("Helvetica",12)
        fig3 = plt.figure(num=1,figsize=(5,4))

        xpos=2
        ypos = 3
        newpage = 1

        for sampleContents in Contents.name:
            curves = []
            for well in self.originalPlate.names:
                if self.originalPlate.wells[well].contents.name == sampleContents:
                    curves.append(well)
            complexDictionary = {}
            meanWellDictionary = {}
            for i, saltConcentration in enumerate(Contents.salt):
                for well in curves:
                    if self.originalPlate.wells[well].contents.salt == saltConcentration:
                        if self.originalPlate.wells[well].complex:
                            complexDictionary[i] = True
                        elif i not in complexDictionary:
                            complexDictionary[i] = False
                        if well in self.delCurves:
                            plt.plot(self.originalPlate.wells[well].temperatures,self.originalPlate.wells[well].fluorescence\
                            , 'grey')
                        # If the curve is not monotonic
                        elif self.originalPlate.wells[well].mono == False:
                            #if it is complex, plot it dotted
                            if self.originalPlate.wells[well].complex:
                                plt.plot(self.originalPlate.wells[well].temperatures,self.originalPlate.wells[well].fluorescence\
                                , COLOURS[i],linestyle=":")
                            #if not complex, plot it normally
                            else:
                                plt.plot(self.originalPlate.wells[well].temperatures,self.originalPlate.wells[well].fluorescence\
                                , COLOURS[i])
                        else:
                            # IF the curve is monotonic it is plotted with a 
                            # dashed line as it is not used to determin Tm
                            plt.plot(self.originalPlate.wells[well].temperatures,self.originalPlate.wells[well].fluorescence\
                            , COLOURS[i],linestyle="--")
                            
                        meanWellDictionary[i] = findKey(well,self.plate.meandict)
                        
            plt.axis([20,100,0.001,0.015])
            plt.gca().axes.get_yaxis().set_visible(False)
            imgdata = cStringIO.StringIO()
            fig3.savefig(imgdata, format='png',dpi=140)
            imgdata.seek(0)  # rewind the data
            Image = ImageReader(imgdata)
            pdf.drawImage(Image, cm+(xpos % 2)*9.5*cm,22.5*cm - (ypos % 3)*9*cm , 8*cm, 6*cm)
            pdf.setFillColor("black")
            pdf.setFont("Helvetica",12)
            pdf.drawString(cm+(xpos % 2)*9.5*cm,22*cm - (ypos % 3)*9*cm ,"Condition: " + sampleContents)
            pdf.drawString(cm+(xpos % 2)*9.5*cm,21.5*cm - (ypos % 3)*9*cm ,"Salt:")
            pdf.drawString(cm+(xpos % 2)*9.5*cm,21*cm - (ypos % 3)*9*cm ,"Tm: ")
            drawdpH = False
            for i in range(len(Contents.salt)):
                pdf.setFillColor(COLOURS[i])
                pdf.setFont("Helvetica",10)

                pdf.drawString(2*cm+((i+3)%3)*2.5*cm+(xpos % 2)*9.5*cm,21.5*cm -(i/3)*cm - (ypos % 3)*9*cm,Contents.salt[i])
                if complexDictionary[i]:
                    if self.wells[meanWellDictionary[i]].Tm != None:
                        if self.wells[meanWellDictionary[i]].TmError != None:
                            pdf.drawString(2*cm+((i+3)%3)*2.5*cm+(xpos % 2)*9.5*cm,21*cm - (ypos % 3)*9*cm,str(round(self.wells[meanWellDictionary[i]].Tm,2))+" (+/-"+str(round(self.wells[meanWellDictionary[i]].TmError,2))+")^")
                        else:
                            pdf.drawString(2*cm+((i+3)%3)*2.5*cm+(xpos % 2)*9.5*cm,21*cm - (ypos % 3)*9*cm,str(round(self.wells[meanWellDictionary[i]].Tm,2))+"^")
                    else:
                        pdf.drawString(2*cm+((i+3)%3)*2.5*cm+(xpos % 2)*9.5*cm,21*cm - (ypos % 3)*9*cm,"None")
                else:
                    if self.wells[meanWellDictionary[i]].Tm != None:
                        if self.wells[meanWellDictionary[i]].TmError != None:
                            pdf.drawString(2*cm+((i+3)%3)*2.5*cm+(xpos % 2)*9.5*cm,21*cm - (ypos % 3)*9*cm,str(round(self.wells[meanWellDictionary[i]].Tm,2))+" (+/-"+str(round(self.wells[meanWellDictionary[i]].TmError,2))+")")
                        else:
                            pdf.drawString(2*cm+((i+3)%3)*2.5*cm+(xpos % 2)*9.5*cm,21*cm - (ypos % 3)*9*cm,str(round(self.wells[meanWellDictionary[i]].Tm,2)))
                    else:
                        pdf.drawString(2*cm+((i+3)%3)*2.5*cm+(xpos % 2)*9.5*cm,21*cm - (ypos % 3)*9*cm,"None")
                if self.wells[meanWellDictionary[i]].contents.dpH != None and self.wells[meanWellDictionary[i]].contents.dpH != "" and self.wells[meanWellDictionary[i]].Tm != None and self.wells[meanWellDictionary[i]].contents.pH != None:
                    pdf.drawString(2*cm+((i+3)%3)*2.5*cm+(xpos % 2)*9.5*cm,20*cm - (ypos % 3)*9*cm,str(round(float(self.wells[meanWellDictionary[i]].contents.pH)+(self.wells[meanWellDictionary[i]].contents.dpH*(self.wells[meanWellDictionary[i]].Tm-20)),2)))
                    pdf.setFillColor("black")
                    if drawdpH ==False:
                        pdf.drawString(2*cm+(xpos % 2)*9.5*cm,20.5*cm - (ypos % 3)*9*cm,"Adjusted pH at Tm: "+self.wells[meanWellDictionary[i]].contents.pH+" at 20C")
                        drawdpH = True
            drawdpH = False
            xpos +=1
            if newpage % 2 == 0:
                ypos +=1
            if newpage % 6 == 0:
                pdf.showPage()
                pdf.setFont("Helvetica",9)
                pdf.drawString(cm, 1.3*cm,"Curves drawn with dashed lines are monotonic and excluded from Tm calculations")
                pdf.drawString(cm, 0.9*cm,"Curves with complex melt transitions are marked (^) and are drawn with a dotted line")
                pdf.drawString(cm, 0.5*cm,"Curves coloured grey are outliers, and are excluded from Tm calculations")
                pdf.setFont("Helvetica",12)
            newpage += 1 
            plt.close()
            fig3 = plt.figure(num=1,figsize=(5,4))

        plt.close()        


        pdf.save()

        return
    
    def computeTms(self):
        """
        To find the Tm of a mean curve we look at the curves that ocmprised it 
        (not including the discarded curves ofcourse) and find the Tm of each of
        those. We then find the mean and sd of these Tms which gives us a Tm estimate, 
        along with an error estimate
        
        When finding the Tm of a curve that is not a mean curve, the error is set to 
        None
        """
        #most calculating is done by getting the mean and sd of replicate Tms,
        #making this the useful part
        for well in self.originalPlate.names:
            #sets own mono instance variable to apropriate state
            self.originalPlate.wells[well].isMonotonic()
            if self.originalPlate.wells[well].mono == False:
                self.originalPlate.wells[well].computeTm()

        for well in self.plate.names:
                tms = [self.originalPlate.wells[x].Tm for x in self.plate.meandict[well]  if x not in self.delCurves and not self.originalPlate.wells[x].mono]
                complexs = [self.originalPlate.wells[x].complex for x in self.plate.meandict[well]  if x not in self.delCurves and not self.originalPlate.wells[x].mono]
                for data in complexs:
                    if data:
                        self.wells[well].complex = True
                self.wells[well].Tm , self.wells[well].TmError = rh.meanSd(tms)
                if len(tms) == 1:
                    self.wells[well].TmError = None
        return

    
class DSFPlate:
    
    def __init__(self, fluorescenceXLS, labelsXLS):
        """
        A container class that keeps all of the data for a DSF experiment
        organized in an ordered format.  This can be used to quickly
        pull out melt curves for wells, along with their contents.
        Lists of replicate can be found from dictionaries, with key 
        replicate well->mean well or key mean well->replicate wells
        
        +fluorescenceXLS is the .xls file that stores all the melt curve data
        
        +labelsXLS is the .xls file that stores the conditions, and pH if given a custom
        input file
        """
            
        #Open excel workbook and first sheet for both the RFU file and the summary/contents file
        wbData = xlrd.open_workbook(fluorescenceXLS)
        shData = wbData.sheet_by_index(0)
        wbContents = xlrd.open_workbook(labelsXLS)
        shContents = wbContents.sheet_by_index(0)
        
        #names are in order and used to iterate over, wells has each name as a key
        self.names = []
        self.wells = {}
        
        #creating dictionary corresponding to replicates on this plate
        self.repDict = {}
        #creating a dictionary that references every mean plate well, to those that it came from
        self.meandict = {}
        #Read rows containing names and conditions
        conditionWellNames = shContents.col_values(0, start_rowx=1, end_rowx=None)
        conditionNames = shContents.col_values(1, start_rowx=1, end_rowx=None)
        conditionSalts = shContents.col_values(2, start_rowx=1, end_rowx=None)
        conditionPhs = shContents.col_values(3, start_rowx=1, end_rowx=None)
        
        
        #checks whether summary for a d(pH)/dT column
        conditiondpHdT=[]
        try:
            for i in range(len(conditionWellNames)):
                conditiondpHdT.append(shContents.cell(i,4).value)
        except IndexError:
            for i in range(len(conditionWellNames)):
                conditiondpHdT.append(None)
        
        #saves the list of names in the Plate for future use
        self.names = conditionWellNames
        
        #control well locations, filled when creating replicate dictionary
        #uses names of controls from global list CONTROL_WELL_NAMES
        self.lysozyme = []
        self.noDye = []
        self.proteinAsSupplied = []
        self.noProtein = []
        #the names of the controls are gotten from the global list at the top of the file,
        #note that the order in that list is important

        #TODO get controls from contents map
        #control names when given a custom summary xls file
        for wellName in conditionWellNames:
            if wellName.lower() == CONTROL_WELL_NAMES[0]:#"lysozyme":
                self.lysozyme.append(wellName)
            elif wellName.lower() == CONTROL_WELL_NAMES[1]:#"no dye":
                self.noDye.append(wellName)
            elif wellName.lower() == CONTROL_WELL_NAMES[2]:#"protein as supplied":
                self.proteinAsSupplied.append(wellName)
            elif wellName.lower() == CONTROL_WELL_NAMES[3]:#"no protein":
                self.noProtein.append(wellName)
        
        #TODO FIX DAS ONE
        #dictionary keys are a single well, with value a list of its reps and itself, e.g. 'A2':['A1','A2','A3']
        self.repDict = {}
        for i, well in enumerate(conditionWellNames):
            #create and entry is one doesnt exist already
            if well not in self.repDict.keys():
                self.repDict[well] = []
            for j, conditionTuple in enumerate(zip(conditionNames, conditionSalts, conditionPhs)):
                #if we find another well with the same name, ph and salt, it is a replicate
                if conditionTuple == zip(conditionNames, conditionSalts, conditionPhs)[i]:
                    self.repDict[well].append(conditionWellNames[j])
        print self.repDict
        
        for i,name in enumerate(conditionWellNames):
            #creates a pandas series of each well, with index being temperature, and values fluorescence
            fluoroSeries = pandas.Series(shData.col_values(i+2,start_rowx=1,end_rowx=None),shData.col_values(1,start_rowx=1,end_rowx=None))
            #populate the well dictionary of labels and values
            self.wells[name] = DSFWell(fluoroSeries, conditionNames[i], conditionSalts[i], conditionPhs[i], conditiondpHdT[i])
            
        return
              

class DSFWell:
    
    def __init__(self, fluorescenceSeries, conditionName, conditionSalt, conditionPh, conditiondpHdT):
        """
        A well class that holds all of the data to characterize
        a well in a buffer plate.
        
        +fluorescenceSeries is a pandas Series describing the melt curve for a single well
        
        +dataLabels is a dictionary specifying formulation of the well.
        
        From the series a well is given a name, its list of temperatures, and its curve"""
        #name, temperatures, and curve from data
        self.name = fluorescenceSeries.name
        self.temperatures = fluorescenceSeries.index
        self.fluorescence = [x for x in fluorescenceSeries]
        #the curve is then normalised to have an area below the curve of 1
        count = 0
        for height in self.fluorescence:
            count += height
        self.fluorescence = [x / count for x in self.fluorescence]
        #the forgive monotonic threshold depends on the normalisation of the curve
        self.monoThresh = DEFAULT_MONO_THRESH / count
        #other attributes of the curve are set to false/none until later analysis of the curve
        self.complex = False
        self.mono = False
        #tm and tm error are calulated upon calling the computeTm() method
        self.Tm = None   
        self.TmError = None
        #the contents of the well is contained in an object of Contents inside well
        self.contents = Contents(conditionName, conditionSalt, conditionPh, conditiondpHdT)
        return  

    def isMonotonic(self):
        """
        Checks to see if the normalised curve, and hence original curve, is monotonic non-increasing
        
        This is done by looking at each consecutive point, forgiving an amount of monoThresh, and seeing if it
        it contradicts this, and is higher than the previous. 5 such consecutive points decide the curve 
        is not monotonic, with each decreasing pair of points inbetween, puts the consecutive contradiction 
        count back down by 1.
        """
        #assume monotonic non-increasing
        nonIncreasingMono = True
        contradictions = 0
        prev=self.fluorescence[0]
        for point in self.fluorescence[1:]:
            if point > prev+self.monoThresh:
                #found contradiction, even with forgive threshold, need (5) contradictions with less than (5)
                #non-contradicting points between them to be able to say it is not monotonically non-increasing
                contradictions+=1
                
            elif point < prev+self.monoThresh:
                #lower contradiction counter, if points are non-increasing as assumed, but not below 0
                if contradictions!=0:
                    contradictions-=1
                    
            #when point=previous do nothing with the contradiction counter
                    
            if contradictions == 5:
                #if 5 contradictions, it is definately not monotonic non-increasing
                nonIncreasingMono=False
                break
            prev=point
            
        #if it is non-increasing set mono to true and end now
        if nonIncreasingMono:
            self.mono = True
            return

        #if not non-increasing, then func is not monotonic
        self.mono = False
        return
        
    
    def computeTm(self):
        """
        Finds the temperature value corresponding to the Tm of a well curve. This is done 
        by finding the lowest point on the derivative of the well, fitting a parabola 
        to it and its two neighbouring indices, and taking the minimum (turning point) 
        of the parabola
        
        The well curve is also checked for having a 'complex' shape, in which case the report
        states that the Tm calculated might be unreliable. 
        This is calculated by finding the minimum and the maximum of of the curve, and seeing if there
        are any turning points between the two. If so, the curve shape is complex
        """
        #first step  is finding the derivative series of the well
        x = [x for x in self.temperatures]
        if self.fluorescence == None:
            self.Tm = None
            return
        y = [y for y in self.fluorescence]
    
        xdiff = np.diff(x)
        dydx = -np.diff(y)/xdiff
        #the derivative series, has one less index since there is one fewer differences than points
        seriesDeriv = pandas.Series(dydx, x[:-1])
        
        #now that we have the derivative series, we can find the Tm
        lowestPoint = 0
        lowestPointIndex = None
        
        #gets number of signchanges between max and min of the curve, used to determin if the curve
        #is complex or not
        lowestPoint2 = 1
        lowestIndex2 = None
        highestPoint = 0
        highestIndex = None
        signChangeCount = 0
        previous = None
        for i, value in enumerate(self.fluorescence[:-1]):
            if value > highestPoint:
                highestPoint = value
                highestIndex = i
        if highestIndex == 0 :
            highestPoint = 0
            highestIndex = None
            for i, value in enumerate(self.fluorescence[:-1]):
                if value<lowestPoint2:
                    lowestPoint2 = value
                    lowestIndex2 = i
            for i, value in enumerate(self.fluorescence[:-1]):
                if i < lowestIndex2:
                    continue
                if value > highestPoint:
                    highestPoint = value
                    highestIndex = i
        else:
            for i, value in enumerate(self.fluorescence[:-1]):
                if i > highestIndex:
                    break
                if value<lowestPoint2:
                    lowestPoint2 = value
                    lowestIndex2 = i
        for ind in seriesDeriv.index[lowestIndex2:highestIndex]:
            if previous:
                if seriesDeriv[ind] < 0 and previous > 0:
                    signChangeCount += 1
                if seriesDeriv[ind] > 0 and previous < 0:
                    signChangeCount += 1
                if seriesDeriv[ind] == 0:
                    signChangeCount += 1
            previous = seriesDeriv[ind]

            
        #finding the lowest point and its index on the derivative series
        #only search for Tm up to 90degrees, since last part is hard to predict
            #and often gives false positives
        for ind in seriesDeriv.index[:-20]:
            if seriesDeriv[ind]<lowestPoint:
                lowestPoint = seriesDeriv[ind]
                lowestPointIndex = ind

        #if lowest point is the first or last index, then no curve fit is required
        if lowestPointIndex in [seriesDeriv.index[0],seriesDeriv.index[-1]]:
            tm = lowestPointIndex
            self.Tm = tm
            
            #set complex to true if curve was complex
            if signChangeCount > 0:
                self.complex = True
            return
    
        #could not find any Tm
        if lowestPointIndex == None:
            self.Tm = None
            
            #if no tm, the curve cannot be complex, set to false (should already be false), before returning
            self.complex = False
            return     
        
        #the indices either side of the lowest index
        leftIndex = lowestPointIndex - 0.5
        rightIndex = lowestPointIndex + 0.5
        
        #matrices used to fit a parabola to the 3 points
        Y=[seriesDeriv[leftIndex],
           seriesDeriv[lowestPointIndex],
           seriesDeriv[rightIndex]]
           
        A=[[leftIndex**2,   leftIndex,   1],
           [lowestPointIndex**2, lowestPointIndex, 1],
           [rightIndex**2,  rightIndex,  1]]
           
        #solve for b, in the form Y=Ab
        (a,b,c) = np.linalg.solve(A,Y)
        
        #initialise tm to left most point of relevant curve
        tm=seriesDeriv[leftIndex]
        tmValue=0
        #make tm the lowest point on the fitted parabola rounded to nearest 0.01
        for x in np.arange(leftIndex,rightIndex,0.01):
            point = (a*(x**2) + b*x + c)
            if tmValue > point:
                tmValue = point
                tm = x
        self.Tm = tm
        
        #again check for complex shape before returning
        if signChangeCount > 0:
                self.complex = True
        averagePoint = (lowestPoint2 +highestPoint) / 2
        i = 0
        while self.fluorescence[i]<averagePoint:
            i += 1;

        # estimates tm by another method and if the difference is too large the curve is considred complex
        if (i/2.0+20 -self.Tm)**2 > 5**2:
            self.complex=True
        
        return

class Contents:
    """
    Keeps a static variable of all the different salt concentrtions and names in 
    the experiment as this is used for plotting in the report
    """
    salt = []
    name = []
    def __init__(self, conditionName, conditionSalt, conditionPh, conditiondpHdT):
        """
        Struct class for the contents of a well, stores the name, pH, salt, and
        dpH of each well.
        
        +sample name  is the contents as given by the contents map, or pcrd summary file
        
        +dpH is the dph as given in the contents file, otherwise its set to None
        """
        self.name = conditionName
        self.salt = conditionSalt
        self.pH = conditionPh
        self.dpH = conditiondpHdT
        
        #add names and salts to the static lists if not already present
        if self.name not in Contents.name:
            Contents.name.append(self.name)
        if self.salt not in Contents.salt:
            Contents.salt.append(self.salt)
        return

        

#==============Functions that can be used to replicate thresholds==============#
def determineOutlierThreshold(listOfLysozymeWellNames):
    """
    Used to reproduce the threshold when determining what curves
    are outliers from a group of replicates
    e.g. listOfLysozymeWellNames = ["A1","A2","A3"] (as on bufferscreen9)
    """
    lysozyme=[]
    results = []
    # TODO should not be hardcoded
    pathrfu = "../UropCrystallisation/data/bufferscreen9/rfuResults/xlsx"
    files = os.listdir(pathrfu)
    pathrfu = pathrfu + "/"
    for data in files:
        plate = DSFPlate(pathrfu+data,"../UropCrystallisation/data/Content_map.xlsx")
        for well in listOfLysozymeWellNames:
            lysozyme.append(plate.wells[well].fluorescence)
    for pair in combinations(lysozyme,2):
        results.append(sqrDiffWellFluoro(pair[0],pair[1]))
    total = 0
    for num in results:
        total+=num
    return total/len(results)
    
def lysozymeAllTm(pathder,wells):
    """
    Calculates the mean and standard deviation of the lysozyme Tm, used when checking controls
    in the experiment
    
    +pathder is the path to the folder containing all and only the rfu derivative
    results from a lot of experiments which contain lysozyme in some wells
    
    +wells is  a list of the wells on the plate that coorrespond to the location
    of lysozyme on the plate e.g. ["A1","A2","A3"]
    """
    Tms=[]
    pathder = pathder + "/"
    for data in os.listdir(pathder):
        df = pandas.DataFrame.from_csv(pathder+data,index_col=1)
        df.pop("Unnamed: 0")
        for well in wells:
            series = df[well]
            mini = 100
            Tm = -1
            for x in series.index:
                if series[x] < mini:
                    mini = series[x]
                    Tm = x
            Tms.append(Tm)
    return rh.meanSd(Tms)
#==============================================================================#


#useful small functions used throughout
#================================================#
def fixPcrdXlsxFile(xlsx):
    """
    Function that fixes an xlsx file error when exporting from a pcrd file. The error is 
    cause by a file inside the pcrd being missnamed, and calling this function on that file
    will fix this, along with keeping the original file saved with '.BAK' appended to its name
    """
    #rename original file, and create new one with same name of the first one
    #delete the backup if it exists already
    if os.path.exists(xlsx+".BAK"):
        os.remove(xlsx+".BAK")
    os.rename(xlsx, xlsx+".BAK")
    source = zf.ZipFile(xlsx+".BAK", "r")
    target = zf.ZipFile(xlsx, "w")
    
    #copy every file from the original xlsx to the new (target) xlsx
    for file in source.filelist:
        #found the wrongly named file
        if "sharedstrings.xml" in file.filename:
            #copy it with correct name
            target.writestr(file.filename[:-17]+"sharedStrings.xml", source.read(file))
        else:
            #otherwise copy file normally
            target.writestr(file.filename, source.read(file))
    source.close()
    target.close()
    
    return
    
    
def findKey(value,dict):
    """
    Takes a dictionary and a value and returns the first key it finds that corresponds to that value
    """
    for key, val in dict.iteritems():
        if value in val:
            return key


def meanCurve(curves):
    """
    Takes in a list of lists (list of individual curves) and returns a single list
    that is the mean curve of the ones supplied
    """
    mean = []
    for i in range(len(curves[0])):
        mean += [0]
    for curve in curves:
        for i in range(len(curve)):
            mean[i] += curve[i]
    return [x /len(curves) for x in mean]


def sqrDiffWellFluoro(fluoro1,fluoro2):
    """
    Gets the sum of squared differences between every pair of points in two fuorescence curves
    """
    #TODO LOG
    dist = 0
    count = 0
    x1 = [math.log(x) for x in fluoro1]
    x2 = [math.log(x) for x in fluoro2]
    while count < len(fluoro1):
        dist += math.pow(x1[count]-x2[count],2)
        count += 1
    return dist
#================================================#


def main():
    """
    parse arguments and call classes/functions needed for DSF analysis
    """
    #opens up windows for user to selec files
    root = Tkinter.Tk()
    root.withdraw()
    #DSF results file
    rfuFilepath = tkFileDialog.askopenfilename(title="Select the DSF results")
    # contents map, or default cfx manager summary file
    contentsMapFilepath = tkFileDialog.askopenfilename(title="Select the contents map")
    try:
        #if the files supplied are xlsx as opposed to xls file the pcrd naming error
        if ".xlsx" in rfuFilepath:
            fixPcrdXlsxFile(rfuFilepath)
        if ".xlsx" in contentsMapFilepath:
            fixPcrdXlsxFile(contentsMapFilepath)
        
        #the analysis
        mydsf = DSFAnalysis()
        mydsf.loadMeltCurves(rfuFilepath,contentsMapFilepath)
        mydsf.analyseCurves()
        
        # generates the report
        name = rfuFilepath.split(".")[0]
        mydsf.generateReport(name+".pdf")

        #also remove the exported xls/xlsx files after meltdown has been run on them
        #find the protein name, then all the files with that name in it, then delete them
        if deleteInputFiles:
            folder = rfuFilepath[:-len(rfuFilepath.split('/')[-1]) - 1]
            proteinName = rfuFilepath.split('/')[-1].split()[0]
            for fl in os.listdir(folder):
                if '.pdf' in fl:
                    continue
                if proteinName in fl:
                    os.remove(folder+'/'+fl)
            
        
    except:
        errors = open("error_log.txt",'w')
        etype, value, tb = sys.exc_info()
        errors.write(''.join(traceback.format_exception(etype, value, tb, None))) 
        root = Tkinter.Tk()
        root.withdraw()
        tkMessageBox.showwarning("Error", "Check error log")
        errors.close()
        
    finally:
        #delete temporary backup files
        if os.path.exists(rfuFilepath+".BAK"):
            os.remove(rfuFilepath+".BAK")
        if os.path.exists(contentsMapFilepath+".BAK"):
            os.remove(contentsMapFilepath+".BAK")
        
    return

#excecutes main() on file run
if __name__ == "__main__":
     main()

#TODO LOG
SIMILARITY_THRESHOLD = determineOutlierThreshold(["A1","A2","A3"])
print SIMILARITY_THRESHOLD;

"""
# Short piece of code for batch analysis of experiments
files = os.listdir("../UropCrystallisation/data/bufferscreen9/rfuResults/xlsx")
total = len(files)
for i, bsc9 in enumerate(files):
    mydsf = DSFAnalysis()
    filepath = "../UropCrystallisation/data/bufferscreen9/rfuResults/xlsx/" + bsc9
    mydsf.loadMeltCurves(filepath,"../UropCrystallisation/data/Content_map.xlsx")
    mydsf.analyseCurves()
    mydsf.generateReport("reports/"+bsc9+".pdf")
    print str(round(i/float(total)* 100,2))  +"%"
"""



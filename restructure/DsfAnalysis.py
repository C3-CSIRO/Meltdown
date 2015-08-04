# -*- coding: utf-8 -*-

import csv
import os
import math
import pandas as pd
import matplotlib.pyplot as plt
import sys
import Tkinter, tkMessageBox
import cStringIO

import replicateHandling as rh
from DsfPlate import DsfPlate, LYSOZYME, SIMILARITY_THRESHOLD
from MeanWell import MeanWell

#reportlab needs to be installed separetly by anaconda, so a messagebox pops up alerting the user if it can't import
try:
    from reportlab.pdfgen import canvas
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.units import cm
    from reportlab.lib.utils import ImageReader
except:
    root = Tkinter.Tk()
    root.withdraw()
    tkMessageBox.showerror("ReportLab not found", "You must use Anaconda to install reportlab before Meltdown can be run")
    sys.exit(1)

#(mean, standard deviation) of lysozyme Tm over ~250 experiments
LYSOZYME_TM_THRESHOLD = (70.8720, 0.7339)

#the running location of this file
RUNNING_LOCATION = os.path.dirname(os.path.realpath(__file__))

#largest tm error before the estimate is considered unreliable
MAX_TM_ERROR_BEFORE_UNRELIABLE = 1.5


class DsfAnalysis:
    def __init__(self, analysisName):
        #initialisations
        self.name = analysisName
        self.plate = None
        self.meanWells = []
        self.contentsHash = {}
        self.controlsHash = {}
        return
        
    def loadCurves(self, dataFilePath, contentsMapFilePath):
        #create the DsfPlate object
        self.plate = DsfPlate(dataFilePath, contentsMapFilePath)
        return
    
    def analyseCurves(self):
        #perform analysis on the plate, note order here is important
        self.plate.computeOutliers()
        self.plate.computeSaturations()
        self.plate.computeMonotonicities()
        self.plate.computeInTheNoises()
        self.plate.computeTms()
        self.plate.computeComplexities()
        #create the mean wells of replicates on the plate
        self.__createMeanWells()
        #create grouped hash for plotting
        self.__createMeanContentsHash()
        #check the controls on the plate
        self.__createControlsHash()
        return
    
    def __createMeanWells(self):
        seen = []
        #loop through each set of replicates
        for wellName in self.plate.repDict.keys():
            if wellName not in seen:
                reps = self.plate.repDict[wellName]
                seen += reps
                #get mean tm and tm error (sd of tms)
                tm, tmError = rh.meanSd([self.plate.wells[w].tm for w in reps if not self.plate.wells[w].isDiscarded])
                complexMean = any([self.plate.wells[w].isComplex for w in reps if not self.plate.wells[w].isDiscarded])
                repsNotDiscarded = sum([(not self.plate.wells[w].isDiscarded) for w in reps])
                contents = self.plate.wells[wellName].contents
                #create a mean well and add it to list
                self.meanWells.append(MeanWell(tm, tmError, complexMean, reps, repsNotDiscarded, contents))
        return
    
    def __createMeanContentsHash(self):
        #loop through each mean well and create a nested contents hash such that
        #{(cv1, ph): {cv2: meanWell}}
        for well in self.meanWells:
            contents = well.contents
            #build the nested hashmaps as we go through the wells, start with cv1,ph tuple
            if (contents.cv1, contents.ph) not in self.contentsHash.keys():
                self.contentsHash[(contents.cv1, contents.ph)] = {}
            #then cv2, which maps to the mean well itself
            if contents.cv2 not in self.contentsHash[(contents.cv1,contents.ph)].keys():
                self.contentsHash[(contents.cv1,contents.ph)][contents.cv2] = well
        return
    
    def __createControlsHash(self):
        #TODO no dye and no protein will always fail, now that we use aitchisons distance and dont divide by length of curve
        #NB they would always have passed before as threshold was set for not dividing by length and therefore massive (oops)
        #initialise the control results
        results = {"lysozyme": "Not Found",
                   "no dye": "Not Found",
                   "no protein": "Not Found"}
        
        #first check if lysozyme control is present on the plate
        if len(self.plate.lysozyme)>0:
            #get the mean well for the lysozyme control, it will have no ph, and no condition variable 2
            lysozymeMeanWell = self.contentsHash[(LYSOZYME,'')]['']
            
            #lysozyme Tm check, only uses mean lysozyme Tm, hence indexing ([0])
            if lysozymeMeanWell.tm > LYSOZYME_TM_THRESHOLD[0] - 2*LYSOZYME_TM_THRESHOLD[1] and\
            lysozymeMeanWell.tm < LYSOZYME_TM_THRESHOLD[0] + 2*LYSOZYME_TM_THRESHOLD[1]:
                results["lysozyme"] = "Passed"
            else:
                results["lysozyme"] = "Failed"
        
        #check if no dye control is present
        if len(self.plate.noDye)>0:
            #create a mean curve out of the replicates that are not outliers
            #initialise the mean curve sum to all zeros
            meanNoDyeCurve = [0 for x in self.plate.wells[self.plate.noDye[0]].temperatures]
            validCurvesInSum = 0
            for wellName in self.plate.noDye:
                well = self.plate.wells[wellName]
                #creates sum of curves being used
                if not well.isOutlier:
                    meanNoDyeCurve = [x+y for x,y in zip(well.fluorescence, meanNoDyeCurve)]
                    validCurvesInSum += 1
            #divide sum to give average curve
            meanNoDyeCurve = [x/validCurvesInSum for x in meanNoDyeCurve]
            
            #read in expected no dye control from file
            noDyeExpected = list(pd.Series.from_csv(RUNNING_LOCATION + "/data/noDyeControl.csv"))
            #if the curves are within required distance from one another, the control is passed
            if rh.aitchisonDistance(meanNoDyeCurve, noDyeExpected) < SIMILARITY_THRESHOLD:
                results["no dye"] = "Passed"
            else:
                results["no dye"] = "Failed"
        
        #check if no protein control is present
        if len(self.plate.noProtein)>0:
            #create a mean curve out of the replicates that are not outliers
            #initialise the mean curve sum to all zeros
            meanNoProteinCurve = [0 for x in self.plate.wells[self.plate.noProtein[0]].temperatures]
            validCurvesInSum = 0
            for wellName in self.plate.noDye:
                well = self.plate.wells[wellName]
                #creates sum of curves being used
                if not well.isOutlier:
                    meanNoProteinCurve = [x+y for x,y in zip(well.fluorescence, meanNoProteinCurve)]
                    validCurvesInSum += 1
            #divide sum to give average curve
            meanNoProteinCurve = [x/validCurvesInSum for x in meanNoProteinCurve]
            
            #read in the expected curve for the no protein control
            noProteinExpected = list(pd.Series.from_csv(RUNNING_LOCATION + "/data/noProteinControl.csv"))
            #if the curves are within required distance from one another, the control is passed
            if rh.aitchisonDistance(meanNoProteinCurve, noProteinExpected) < SIMILARITY_THRESHOLD:
                results["no protein"] = "Passed"
            else:
                results["no protein"] = "Failed"
        
        #save the results hash
        self.controlsHash = results
        return
    
    def produceNormalisedOutput(self, filePath):
        #names of all the wells in sorted order
        sortedWellNames = sorted(self.plate.wells.keys())
        #list of temperatures, taken from first well since all have the same temperature list
        temperatures = self.plate.wells.values()[0].temperatures
        
        with open(filePath, 'w') as fp:
            fWriter = csv.writer(fp, delimiter='\t')
            fWriter.writerow(['Temperature'] + sortedWellNames)
            for i in range(len(temperatures)):
                #start each row with the temperature
                row = [temperatures[i]]
                #create each row as the value at that temperature on each well
                for wellName in sortedWellNames:
                    row.append(self.plate.wells[wellName].fluorescence[i])
                
                #write to the file
                fWriter.writerow(row)
        return
    
    def generateReport(self, outputFilePath, version):
        #===================# headings and image #===================#
        #initialise the output pdf and print the heading and name of experiment
        pdf = canvas.Canvas(outputFilePath,pagesize=A4)
        pdf.setFont("Helvetica-Bold",16)
        pdf.drawString(cm,28*cm,"MELTDOWN v" + version)
        pdf.setFont("Helvetica",16)
        pdf.drawString(7*cm,28*cm,"Melt Curve Analysis")
        if len(self.name) < 40:
            pdf.drawString(cm,27*cm, self.name)
        else:
            pdf.drawString(cm,27*cm, self.name[:41] + '...')

        #put the csiro image in the top right
        pdf.drawImage(RUNNING_LOCATION + "/data/CSIRO_Grad_RGB_hr.jpg",17*cm,25.5*cm,3.5*cm,3.5*cm)
        
        
        #===================# protein as supplied graph and Tm #===================#
        #create a plot for the protein as supplied control, and plot the curves
        proteinAsSuppliedFigure = plt.figure(num=1,figsize=(5,4))
        for wellName in self.plate.proteinAsSupplied:
            well = self.plate.wells[wellName]
            if well.isDiscarded:
                #discarded curves are dotted
                plt.plot(well.temperatures, well.fluorescence, 'g', linestyle=":")
            elif well.isComplex:
                #complex curves are dashed
                plt.plot(well.temperatures, well.fluorescence, 'g', linestyle="--")
            else:
                #normal curves are full lines
                plt.plot(well.temperatures, well.fluorescence, 'g')
        #hide y axis, as RFU units are arbitrary
        plt.gca().axes.get_yaxis().set_visible(False)
        #put the image on the pdf
        imgdata = cStringIO.StringIO()
        proteinAsSuppliedFigure.savefig(imgdata, format='png',dpi=140)
        imgdata.seek(0)
        Image = ImageReader(imgdata)
        pdf.drawImage(Image, 0, 18*cm, 8*cm, 6*cm)
        plt.close()
        
        #print the tm of the protein as supplied below its graph, if the control was found
        pdf.setFont("Helvetica",10)
        if len(self.plate.proteinAsSupplied) > 0:
            #null strings for ph and condition variable 2 in the contents hash, as that's how controls are read
            meanSuppliedProtein = self.contentsHash[('protein as supplied', '')]['']
            suppliedProteinTm = meanSuppliedProtein.tm
            suppliedProteinTmError = meanSuppliedProtein.tmError
            if suppliedProteinTm != None and meanSuppliedProtein.replicatesNotDiscarded > 1:
                pdf.drawString(cm,17.5*cm, "Protein as supplied: Tm = " +str(round(suppliedProteinTm,2))+"(+/-"+str(round(suppliedProteinTmError,2))+")")
            elif suppliedProteinTm != None:
                pdf.drawString(cm,17.5*cm, "Protein as supplied: Tm = " +str(round(suppliedProteinTm,2)))
            else:
                pdf.drawString(cm,17.5*cm, "Protein as supplied: Tm = N/A")
        else:
            pdf.drawString(cm,17.5*cm, "Protein as supplied: Not Found")
        
        
        #===================# first page summary box (top right) #===================#
        #drawing the summary box to the right of the protein as supplied plot
        pdf.rect(7.75*cm,18.05*cm,12*cm,5.4*cm)
        pdf.setFont("Helvetica-Bold",13)
        pdf.drawString(8*cm,22.75*cm,"Full interpretation of the results requires you to look ")
        pdf.drawString(8*cm,22.25*cm,"at the individual melt curves.")
        
        #find and print percentage of non control wells that were used in Tm calculations
        numberOfNonControlWells = 0
        numberOfNonControlFoundTms = 0
        for well in self.plate.wells.values():
            if not well.contents.isControl:
                if well.tm != None:
                    numberOfNonControlFoundTms += 1
                numberOfNonControlWells += 1
        percentTmsFound = int(round(numberOfNonControlFoundTms/float(numberOfNonControlWells),2)*100)
        pdf.drawString(8*cm,20.5*cm,str(percentTmsFound)+"%")
        pdf.setFont("Helvetica",13)
        pdf.drawString(9.25*cm,20.5*cm,"of curves were used in Tm estimations")
        
        #find the average calculated tm error from all mean replicate tm errors, and print it
        tmErrorSum = 0.0
        numOfTmErrors = 0
        for well in self.meanWells:
            if well.tmError != None:
                tmErrorSum += well.tmError
                numOfTmErrors += 1
        pdf.drawString(8*cm,19.5*cm,"Average estimation of error is")
        if numOfTmErrors != 0:
            avgTmError = round(tmErrorSum/float(numOfTmErrors),1)
        pdf.setFont("Helvetica-Bold",13)
        pdf.drawString(14.1*cm,19.5*cm,str(avgTmError)+" C")

        #print whether the protein as supplied was well behaved
        pdf.setFont("Helvetica",13)
        pdf.drawString(8*cm,18.5*cm, "Protein as supplied is")
        pdf.setFont("Helvetica-Bold",13)
        proteinAsSuppliedIswellBehaved = True
        #not well behaved if any protein as supplied replicate has no Tm, or the tm error of the group is too high
        if len(self.plate.proteinAsSupplied) > 0:
            for wellName in self.contentsHash[('protein as supplied', '')][''].replicates:
                well = self.plate.wells[wellName]
                if well.tm == None:
                    proteinAsSuppliedIswellBehaved = False
                    break
            if suppliedProteinTmError >= MAX_TM_ERROR_BEFORE_UNRELIABLE:
                proteinAsSuppliedIswellBehaved = False
            #print outcome
            if proteinAsSuppliedIswellBehaved:
                pdf.drawString(12.5*cm,18.5*cm,"well behaved")
            else:
                pdf.drawString(12.5*cm,18.5*cm,"not well behaved")
        else:
            pdf.drawString(12.5*cm,18.5*cm,"not found")            
        #whether or not we are considering the summary graph to be unreliable,
        #depends on how protein as supplied behaved, and the average tm estimate error
        if not proteinAsSuppliedIswellBehaved or avgTmError >= MAX_TM_ERROR_BEFORE_UNRELIABLE:
            pdf.drawString(8*cm,21.5*cm,"The summary graph appears to be unreliable")

        
        #===================# controls #===================#
        #print out the results of the controls that are checked for
        pdf.setFillColor("blue")
        pdf.setFont("Helvetica",10)
        # lysozyme Tm control check
        pdf.drawString(1*cm,16.5*cm,"Lysozyme Control: " + self.controlsHash["lysozyme"])
        # no dye control check 
        pdf.drawString(1*cm,16*cm,"No Dye Control: " + self.controlsHash["no dye"])
        # no protein control check
        pdf.drawString(1*cm,15.5*cm,"No Protein Control: " + self.controlsHash["no protein"])
        
        
        #===================# first page summary graph #===================#
        #gets a set of all the condition var 2s in the experiment (without repeats), excluding ones only present in controls
        uniqueCv2s = set([cv2 for cv2Dict in self.contentsHash.values() for cv2 in [key for key in cv2Dict.keys() if not cv2Dict[key].contents.isControl]])
        #gets a sorted by ph list of (condition var 1, ph) tuples. these are unique, and do not include controls
        cv1PhPairs = sorted([key for key in self.contentsHash.keys() if any([not meanWell.contents.isControl for meanWell in self.contentsHash[key].values()])], key=lambda x: x[1])
        
        #turns the tuples into string names to display on the x axis
        xAxisConditionLabels = [pair[0]+"("+str(pair[1])+")" for pair in cv1PhPairs]
        
        #flag for if there were any complex curve tms found, so that the warning is displayed        
        foundUnreliable = False
        #list of plat handles, used in giving the legend the right colours
        legendHandles = []
        #y axis min and max initialisations, these are changed based on the highest and lowest Tms
        yAxisMin = yAxisMax = 0
        #save the meanwell which gives the highest Tm, and put this on the page
        highestTmMeanWell = None
        #creates the graph figure
        summaryGraphFigure = plt.figure(num=1,figsize=(10,8))
        
        for cv2 in uniqueCv2s:
            #the normal tms
            tms = []
            #the unreliable tms
            complexTms = []
            
            for cv1,ph in cv1PhPairs:
                #first check if have a well with the specified cv1, cv2, and ph
                try:
                    meanWell = self.contentsHash[(cv1, ph)][cv2]
                    conditionExists = True
                except KeyError:
                    #given (condition var 1, ph) pair does not have a well with current condition var 2
                    conditionExists = False
                
                #if we found a condition, add it's tm to the right list, and a None to the other
                if conditionExists:
                    newTm = meanWell.tm
                    if meanWell.isComplex:#TODO check when tms are unreliable (if estimate comes from only 1 ?)
                        tms.append(None)
                        complexTms.append(newTm)
                        foundUnreliable = True
                    else:
                        tms.append(newTm)
                        complexTms.append(None)
                #otherwise, add Nones to both lists
                else:
                    newTm = None
                    tms.append(None)
                    complexTms.append(None)
                
                #next we adjust the y axis min and max so that they fit the newly added tms
                if newTm:
                    #if this is the first Tm, change both the min and max to it
                    if yAxisMin == 0 and yAxisMax == 0:
                        yAxisMin = yAxisMax = newTm
                        highestTmMeanWell = meanWell
                    #otherwise update the y axis min or max according if required
                    elif newTm < yAxisMin:
                        yAxisMin = newTm
                    elif newTm > yAxisMax:
                        yAxisMax = newTm
                        highestTmMeanWell = meanWell
                
            #plot the tms and the complex tms, and add the non-complex ones to the legend handles
            handle, = plt.plot([x for x in range(len(xAxisConditionLabels))], tms, color=self.plate.cv2ColourDict[cv2], marker="o", linestyle="None")
            plt.plot([x for x in range(len(xAxisConditionLabels))], complexTms, color=self.plate.cv2ColourDict[cv2], marker="d", linestyle="None")
            legendHandles.append(handle)
        
        #set the min and max of the y axis, centre around protein as supplied Tm, if it's present
        if len(self.plate.proteinAsSupplied) > 0 and suppliedProteinTm != None:
            #draw a horizontal dashed red line for the protein as supplied Tm
            plt.axhline(suppliedProteinTm, 0, 1, linestyle="--", color="red")
            #centre around protein as supplied Tm
            distEitherSideOfSuppliedTm = max(math.fabs(suppliedProteinTm - yAxisMin), math.fabs(suppliedProteinTm - yAxisMax))
            plt.axis([-1, len(xAxisConditionLabels), suppliedProteinTm - distEitherSideOfSuppliedTm - 1, suppliedProteinTm + distEitherSideOfSuppliedTm + 1])
        else:
            #no protein as supplied Tm, just use calculated y axis min and max
            plt.axis([-1, len(xAxisConditionLabels), yAxisMin - 1, yAxisMax + 1])
            
        #label the axes
        plt.ylabel('Tm')
        plt.xticks([x for x in range(len(xAxisConditionLabels))], xAxisConditionLabels, rotation="vertical")
        
        #change the padding above the graph when legend gets bigger (i.e. there are more condition variable 2's)
        plt.gcf().subplots_adjust(bottom=0.35, top=0.85 - 0.035*(int(len(uniqueCv2s)/3)))
        #plot the legend
        plt.legend(legendHandles, uniqueCv2s, loc='lower center', bbox_to_anchor=(0.5, 1), ncol=3, fancybox=True, shadow=False, numpoints=1)
        
        #save the graph and print it on the pdf
        imgdata = cStringIO.StringIO()
        summaryGraphFigure.savefig(imgdata, format='png',dpi=180)
        imgdata.seek(0)
        Image = ImageReader(imgdata)
        pdf.drawImage(Image, cm, 4*cm, 16*cm, 11*cm)
        plt.close()

        #if there were any Tms computed as unreliable, print a warning above the graph
        pdf.setFillColor("black")
        if foundUnreliable:
            pdf.drawString(5.9*cm, 14.2*cm, "Tms drawn in diamonds may be unreliable")
        #if supplied protein has dashed line drawn for Tm, label it
        if len(self.plate.proteinAsSupplied) > 0 and suppliedProteinTm != None:
            pdf.drawString(15.5*cm,10.4*cm,"Protein as supplied")
        
        #if we found a highest Tm, print the condition that gave it, and it's Tm below the summary graph
        if highestTmMeanWell:
            pdf.setFont("Helvetica-Bold",12)
            if highestTmMeanWell.tmError != None:
                pdf.drawString(3*cm,2.6*cm,"Highest Tm = " + str(round(highestTmMeanWell.tm,2)) + " +/- " + str(round(highestTmMeanWell.tmError,2)))
            else:
                pdf.drawString(3*cm,2.6*cm,"Highest Tm = " + str(round(highestTmMeanWell.tm,2)))
            pdf.drawString(3*cm,2*cm,"("+highestTmMeanWell.contents.cv1+" / "+highestTmMeanWell.contents.cv2+")")
            pdf.setFont("Helvetica",12)
        
        
        #===================# individual condition graphs #===================#
        #start a new page
        pdf.showPage()
        pdf.setFont("Helvetica",10)
        #getting the y axis scale to be the same over all the boxes, finds appropriate min and max
        overallWellNormalisedMin = overallWellNormalisedMax = 0
        for well in self.plate.wells.values():
            #found first value, initialise min and max to this
            if overallWellNormalisedMin == 0 and overallWellNormalisedMax == 0:
                overallWellNormalisedMin = overallWellNormalisedMax = well.wellNormalisedMin
            #adjust the min and max if we have found a new smallest, or highest, value
            if well.wellNormalisedMin < overallWellNormalisedMin:
                overallWellNormalisedMin = well.wellNormalisedMin
            if well.wellNormalisedMax > overallWellNormalisedMax:
                overallWellNormalisedMax = well.wellNormalisedMax
        #save for use when plotting
        minYValue = overallWellNormalisedMin
        maxYValue = overallWellNormalisedMax
        #this is the added padding to the y axis when plotted
        paddingSize = (maxYValue - minYValue) * 0.05

        #variables used for knowing where to plot the next graph
        numberOfGraphsDrawn = 0
        xpos=2
        #graph positions will depend on the number of condition variable 2's
        if len(uniqueCv2s) < 6:
            #number of images to fit on a page
            maxGraphsPerPage = 6
            ySize = 9.2
            ypos = 3
            yNum = 3
        elif len(uniqueCv2s) < 13:
            maxGraphsPerPage = 4
            ySize = 13.8
            ypos = 2
            yNum = 2
        else:
            maxGraphsPerPage = 2
            ySize = 0
            ypos = 1
            yNum = 1
        
        #first we loop the condition variable 1 / pH pairs
        for cv1, ph in cv1PhPairs:
            #the plotting figure used for current cv1/ph image
            singleConditionFigure = plt.figure(num=1,figsize=(5,4))
            #start printing the tms at the top of the list, and assume no dph/dt is present for condition to begin with
            tmPrintOffset = 0
            hasDphdt = False
            #loop condition variable 2's present for the cv1/ph pair
            for cv2 in sorted(self.contentsHash[(cv1, ph)].keys()):
                #find the associated mean well
                meanWell = self.contentsHash[(cv1, ph)][cv2]
                
                #plot the curve for the well
                wells = [self.plate.wells[wellName] for wellName in meanWell.replicates]
                for well in wells:
                    #dotted line for discarded curves
                    if well.isDiscarded:
                        plt.plot(well.temperatures, well.fluorescence, self.plate.cv2ColourDict[cv2],linestyle=":")
                    #dashed line for complex curves
                    elif well.isComplex:
                        plt.plot(well.temperatures, well.fluorescence, self.plate.cv2ColourDict[cv2],linestyle="--")
                    #full line for normal curves
                    else:
                        plt.plot(well.temperatures, well.fluorescence, self.plate.cv2ColourDict[cv2])
                
                #print the tm calculated for the condition
                pdf.setFillColor(self.plate.cv2ColourDict[cv2])
                pdf.drawString(cm+(xpos % 2)*9.5*cm,22*cm - (ypos % yNum)*ySize*cm - tmPrintOffset*0.5*cm,cv2)
                #if condition was complex, append '^'
                if meanWell.isComplex:
                    #if Tm is calculabe, print it
                    if meanWell.tm != None:
                        #if Tm estimate is from more than 1 replicate, print the Tm error aswell
                        if meanWell.replicatesNotDiscarded > 1:
                            pdf.drawString(4.25*cm+(xpos % 2)*9.5*cm,22*cm - (ypos % yNum)*ySize*cm - tmPrintOffset*0.5*cm ,str(round(meanWell.tm,2))+" (+/-"+str(round(meanWell.tmError,2))+")^")
                        #estimate from only one replicate, do not print Tm error
                        else:
                            pdf.drawString(4.25*cm+(xpos % 2)*9.5*cm,22*cm - (ypos % yNum)*ySize*cm - tmPrintOffset*0.5*cm ,str(round(meanWell.tm,2))+"^")
                    #no calculabe Tm, print 'None' instead
                    else:
                        pdf.drawString(4.25*cm+(xpos % 2)*9.5*cm,22*cm - (ypos % yNum)*ySize*cm - tmPrintOffset*0.5*cm ,"None")
                #not complex, append nothing
                else:
                    #if Tm is calculabe, print it
                    if meanWell.tm != None:
                        #if Tm estimate is from more than 1 replicate, print the Tm error aswell
                        if meanWell.replicatesNotDiscarded > 1:
                            pdf.drawString(4.25*cm+(xpos % 2)*9.5*cm,22*cm - (ypos % yNum)*ySize*cm - tmPrintOffset*0.5*cm ,str(round(meanWell.tm,2))+" (+/-"+str(round(meanWell.tmError,2))+")")
                        #estimate from only one replicate, do not print Tm error
                        else:
                            pdf.drawString(4.25*cm+(xpos % 2)*9.5*cm,22*cm - (ypos % yNum)*ySize*cm - tmPrintOffset*0.5*cm ,str(round(meanWell.tm,2)))
                    #no calculabe Tm, print 'None' instead
                    else:
                        pdf.drawString(4.25*cm+(xpos % 2)*9.5*cm,22*cm - (ypos % yNum)*ySize*cm - tmPrintOffset*0.5*cm ,"None")
                
                #if any of the wells plotted on this graph will have a calculabe adjusted ph, calculate and print it
                if meanWell.contents.dphdt != '' and meanWell.contents.ph != '' and meanWell.tm != None:
                    adjustedPh = str(round(float(meanWell.contents.ph)+(meanWell.contents.dphdt*(meanWell.tm-20)),2))
                    pdf.drawString(7*cm+(xpos % 2)*9.5*cm,22*cm - (ypos % yNum)*ySize*cm - tmPrintOffset*0.5*cm, adjustedPh)
                    #set flag that at least one of the curves on this graph had its adjust ph calculated
                    hasDphdt = True
                #incrememnt the tm printing offset, for the next condition variable 2
                tmPrintOffset += 1
            
            #finalise the plot's axes
            plt.ylim(minYValue-paddingSize,maxYValue+paddingSize)
            plt.gca().axes.get_yaxis().set_visible(False)
            #save the graph figure, and print it to the pdf
            imgdata = cStringIO.StringIO()
            singleConditionFigure.savefig(imgdata, format='png',dpi=140)
            imgdata.seek(0)
            Image = ImageReader(imgdata)
            pdf.drawImage(Image, cm+(xpos % 2)*9.5*cm,23.5*cm - (ypos % yNum)*ySize*cm , 8*cm, 6*cm)
            plt.close()
            
            #print the condition name, and headings for calculated data
            pdf.setFillColor("black")
            pdf.setFont("Helvetica",12)
            pdf.drawString(cm+(xpos % 2)*9.5*cm,23*cm - (ypos % yNum)*ySize*cm ,cv1 + " (" + str(ph)+")")
            pdf.setFont("Helvetica",10)
            pdf.drawString(cm+(xpos % 2)*9.5*cm,22.5*cm - (ypos % yNum)*ySize*cm ,"Grouped by")
            pdf.drawString(4.25*cm+(xpos % 2)*9.5*cm,22.5*cm - (ypos % yNum)*ySize*cm ,"Tm")
            #only print the adjusted ph heading, if one of the conditions has had it calculated for the current graph
            if hasDphdt:
                pdf.setFillColor("black")
                pdf.drawString(7*cm+(xpos % 2)*9.5*cm,22.5*cm - (ypos % yNum)*ySize*cm ,"Adjusted pH at Tm")
                hasDphdt = False
            
            #udpate the plotting position variables accordingly
            xpos +=1
            if numberOfGraphsDrawn % 2 == 1:
                ypos +=1
            numberOfGraphsDrawn += 1 
            
            #if we have started a new page, print the plotting descriptions at the bottom of the page
            if numberOfGraphsDrawn % maxGraphsPerPage == 1:
                pdf.setFont("Helvetica",9)
                ##pdf.drawString(cm, 1.3*cm,"Curves drawn with dashed lines are unable to be analysed (monotonic, saturated, in the noise, and outliers)")
                ##pdf.drawString(cm, 0.9*cm,"and are excluded from Tm calculations")
                pdf.drawString(cm, 0.9*cm,"Monotonic, saturated, in the noise, and outlier curves are dotted, and excluded from Tm calculations")
                pdf.drawString(cm, 0.5*cm,"Curves drawn with dashed lines have unreliable Tm estimates")
                ##pdf.drawString(cm, 0.5*cm,"Curves drawn with dotted lines have unreliable estimates for Tms")
                pdf.setFont("Helvetica",10)

            #we have filled a page, move on to the next one
            if numberOfGraphsDrawn % maxGraphsPerPage == 0:
                #starts new page
                pdf.showPage()

        #save the pdf    
        pdf.save()
        return

def main():
    root = Tkinter.Tk()
    root.withdraw()
    tkMessageBox.showwarning("Inncorrect Usage", "Please read the instructions on how to run Meltdown")
    return
    
    
if __name__ == "__main__":
    main()



















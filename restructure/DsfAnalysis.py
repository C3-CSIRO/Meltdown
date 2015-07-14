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
                contents = self.plate.wells[wellName].contents
                #create a mean well and add it to list
                self.meanWells.append(MeanWell(tm, tmError, complexMean, reps, contents))
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
            noDyeExpected = list(pd.Series.from_csv(RUNNING_LOCATION + "\\data\\noDyeControl.csv"))
            #if the curves are within required distance from one another, the control is passed
            if rh.sqrdiff(meanNoDyeCurve, noDyeExpected) < SIMILARITY_THRESHOLD:
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
            noProteinExpected = list(pd.Series.from_csv(RUNNING_LOCATION + "\\data\\noProteinControl.csv"))
            #if the curves are within required distance from one another, the control is passed
            if rh.sqrdiff(meanNoProteinCurve, noProteinExpected) < SIMILARITY_THRESHOLD:
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
        pdf.drawImage(RUNNING_LOCATION + "\\data\\CSIRO_Grad_RGB_hr.jpg",17*cm,25.5*cm,3.5*cm,3.5*cm)
        
        
        
        #create a plot for the protein as supplied control, and plot the curves
        proteinAsSuppliedFigure = plt.figure(num=1,figsize=(5,4))
        for wellName in self.plate.proteinAsSupplied:
            well = self.plate.wells[wellName]
            if well.isDiscarded:
                #discarded curves are gray
                plt.plot(well.temperatures, well.fluorescence, 'grey')
            elif well.isComplex:
                #complex curves are dashed green
                plt.plot(well.temperatures, well.fluorescence, 'g', linestyle="--")
            else:
                #normal curves are just green
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
            suppliedProteinTm = self.contentsHash[('protein as supplied', '')][''].tm
            suppliedProteinTmError = self.contentsHash[('protein as supplied', '')][''].tmError
            if suppliedProteinTm != None and suppliedProteinTmError != None:
                pdf.drawString(cm,17.5*cm, "Protein as supplied: Tm = " +str(round(suppliedProteinTm,2))+"(+/-"+str(round(suppliedProteinTmError,2))+")")
            elif suppliedProteinTm != None:
                pdf.drawString(cm,17.5*cm, "Protein as supplied: Tm = " +str(round(suppliedProteinTm,2)))
            else:
                pdf.drawString(cm,17.5*cm, "Protein as supplied: Tm = N/A")
        else:
            pdf.drawString(cm,17.5*cm, "Protein as supplied: Not Found")
        
        
        
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

        
        
        #print out the results of the controls that are checked for
        pdf.setFillColor("blue")
        pdf.setFont("Helvetica",10)
        # lysozyme Tm control check
        pdf.drawString(1*cm,16.5*cm,"Lysozyme Control: " + self.controlsHash["lysozyme"])
        # no dye control check 
        pdf.drawString(1*cm,16*cm,"No Dye Control: " + self.controlsHash["no dye"])
        # no protein control check
        pdf.drawString(1*cm,15.5*cm,"No Protein Control: " + self.controlsHash["no protein"])
        
        
        
        
        
        
        
        
        #plotting the summary graph at the bottom of the first page
        #gets a set of all the condition var 2s in the experiment (without repeats)
        uniqueCv2s = set([cv2 for cv2Dict in self.contentsHash.values() for cv2 in cv2Dict.keys()])#TODO this should not include cv2's only found in controls
        #gets a sorted by ph list of (condition var 1, ph) tuples. these are unique
        cv1PhPairs = sorted(self.contentsHash.keys(), key=lambda x: x[1])#TODO this should not include controls
        
        
        #turns the tuples into string names to display on the x axis
        xAxisConditionLabels = [pair[0]+"("+str(pair[1])+")" for pair in cv1PhPairs]
        
        #flag for if there were any complex curve tms found, so that the warning is displayed        
        foundUnreliable = False
        #list of plat handles, used in giving the legend the right colours
        legendHandles = []
        #y axis min and max initialisations, these are changed based on the highest and lowest Tms
        yAxisMin = yAxisMax = 0
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
                #now make sure we aren't looking at a control well, since we don't want to plot those
                if conditionExists and meanWell.contents.isControl:
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
                    #otherwise update the y axis min or max according if required
                    elif newTm < yAxisMin:
                        yAxisMin = newTm
                    elif newTm > yAxisMax:
                        yAxisMax = newTm
                
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
            plt.axis([-1, len(xAxisConditionLabels), suppliedProteinTm - distEitherSideOfSuppliedTm - 5, suppliedProteinTm + distEitherSideOfSuppliedTm + 5])
        else:
            #no protein as supplied Tm, just use calculated y axis min and max
            plt.axis([-1, len(xAxisConditionLabels), yAxisMin - 5, yAxisMax + 5])
            
        #label the axes
        plt.ylabel('Tm')
        plt.xticks([x for x in range(len(xAxisConditionLabels))], xAxisConditionLabels, rotation="vertical")
        
        #change the padding above the graph when legend get bigger (i.e. there are more condition variable 2's)
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
            pdf.drawString(7.9*cm, 14.2*cm, "Tms drawn in diamonds may be unreliable")
            #if supplied protein has dashed line drawn for Tm, label it
        if len(self.plate.proteinAsSupplied) > 0 and suppliedProteinTm != None:
            pdf.drawString(15.5*cm,10.4*cm,"Protein as supplied") 
        
            
            
        
        
        
        
        
        #===================other pages plotting below
        # Variables used to keep track of where to draw the current graph
        xpos=2
        
        newpage = 1

        if len(Contents.salt) < 6:
            ySize = 9.2
            ypos = 3
            graphNum = 6
            yNum = 3
        elif len(Contents.salt) < 13:
            ySize = 13.8
            ypos = 2
            graphNum = 4
            yNum = 2
        else:
            ySize = 0
            ypos = 1
            graphNum = 2
            yNum = 1

        for sampleContentspH in Contents.name:
            if (newpage-1) % graphNum == 0:
                pdf.showPage()
                pdf.setFont("Helvetica",9)
                pdf.drawString(cm, 1.3*cm,"Curves drawn with dashed lines are unable to be analysed (monotonic, saturated, in the noise, and outliers)")
                pdf.drawString(cm, 0.9*cm,"and are excluded from Tm calculations")
                pdf.drawString(cm, 0.5*cm,"Curves drawn with dotted lines have unreliable estimates for Tms")
                pdf.setFont("Helvetica",12)
            sampleContents = sampleContentspH[0]
            curves = []
            for well in self.originalPlate.names:
                if self.originalPlate.wells[well].contents.name == sampleContents and self.originalPlate.wells[well].contents.pH == sampleContentspH[1]:
                    curves.append(well)
            complexDictionary = {}
            meanWellDictionary = {}
            for i, saltConcentration in enumerate(Contents.salt):
                meanWellDictionary[i] = None
                complexDictionary[i] = False
                for well in curves:
                    if self.originalPlate.wells[well].contents.salt == saltConcentration:
                        if self.originalPlate.wells[well].complex:
                            complexDictionary[i] = True
                            
                        #if curve is in delCurves, plot it dashed. Consists of monotonic curves, curves in the noise,
                        #saturated (flat) curves, and replicate outlier curves
                        if well in self.delCurves:
                            plt.plot(self.originalPlate.wells[well].temperatures,self.originalPlate.wells[well].fluorescence\
                            , COLOURS[i],linestyle="--")
                            
                        #If the curve complex (unreliable Tm), plot it dotted. 
                        #Consists of complex and Tms that are not steep enough
                        elif self.originalPlate.wells[well].complex == True:
                            plt.plot(self.originalPlate.wells[well].temperatures,self.originalPlate.wells[well].fluorescence\
                            , COLOURS[i],linestyle=":")
                            
                        #otherwise plot the graph normally
                        else:
                            plt.plot(self.originalPlate.wells[well].temperatures,self.originalPlate.wells[well].fluorescence\
                            , COLOURS[i])
                            
                        meanWellDictionary[i] = findKey(well,self.plate.meanDict)
                        
            plt.ylim(minYValue-paddingSize,maxYValue+paddingSize)
            plt.gca().axes.get_yaxis().set_visible(False)
            imgdata = cStringIO.StringIO()
            fig3.savefig(imgdata, format='png',dpi=140)
            imgdata.seek(0)  # rewind the data
            Image = ImageReader(imgdata)
            pdf.drawImage(Image, cm+(xpos % 2)*9.5*cm,23.5*cm - (ypos % yNum)*ySize*cm , 8*cm, 6*cm)
            pdf.setFillColor("black")
            pdf.setFont("Helvetica",12)
            pdf.drawString(cm+(xpos % 2)*9.5*cm,23*cm - (ypos % yNum)*ySize*cm ,sampleContents + " (" + str(sampleContentspH[1])+")")
            pdf.setFont("Helvetica",10)
            pdf.drawString(cm+(xpos % 2)*9.5*cm,22.5*cm - (ypos % yNum)*ySize*cm ,"Grouped by")
            pdf.drawString(4.25*cm+(xpos % 2)*9.5*cm,22.5*cm - (ypos % yNum)*ySize*cm ,"Tm")
            drawdpH = False
            offset=0
            for i in range(len(Contents.salt)):
                if Contents.salt[i]:
                    pdf.setFillColor(COLOURS[i])
                    pdf.drawString(cm+(xpos % 2)*9.5*cm,22*cm - (ypos % yNum)*ySize*cm - offset*0.5*cm,Contents.salt[i])
                    if complexDictionary[i]:
                        if meanWellDictionary[i] != None and self.wells[meanWellDictionary[i]].Tm != None:
                            if self.wells[meanWellDictionary[i]].TmError != None:
                                pdf.drawString(4.25*cm+(xpos % 2)*9.5*cm,22*cm - (ypos % yNum)*ySize*cm - offset*0.5*cm ,str(round(self.wells[meanWellDictionary[i]].Tm,2))+" (+/-"+str(round(self.wells[meanWellDictionary[i]].TmError,2))+")^")
                            else:
                                pdf.drawString(4.25*cm+(xpos % 2)*9.5*cm,22*cm - (ypos % yNum)*ySize*cm - offset*0.5*cm ,str(round(self.wells[meanWellDictionary[i]].Tm,2))+"^")
                        else:
                            pdf.drawString(4.25*cm+(xpos % 2)*9.5*cm,22*cm - (ypos % yNum)*ySize*cm - offset*0.5*cm ,"None")
                    else:
                        if meanWellDictionary[i] != None and self.wells[meanWellDictionary[i]].Tm != None:
                            if self.wells[meanWellDictionary[i]].TmError != None:
                                pdf.drawString(4.25*cm+(xpos % 2)*9.5*cm,22*cm - (ypos % yNum)*ySize*cm - offset*0.5*cm ,str(round(self.wells[meanWellDictionary[i]].Tm,2))+" (+/-"+str(round(self.wells[meanWellDictionary[i]].TmError,2))+")")
                            else:
                                pdf.drawString(4.25*cm+(xpos % 2)*9.5*cm,22*cm - (ypos % yNum)*ySize*cm - offset*0.5*cm ,str(round(self.wells[meanWellDictionary[i]].Tm,2)))
                        else:
                            pdf.drawString(4.25*cm+(xpos % 2)*9.5*cm,22*cm - (ypos % yNum)*ySize*cm - offset*0.5*cm ,"None")
                    if meanWellDictionary[i] != None and self.wells[meanWellDictionary[i]].contents.dpH != None and self.wells[meanWellDictionary[i]].contents.dpH != "" and self.wells[meanWellDictionary[i]].Tm != None and self.wells[meanWellDictionary[i]].contents.pH != None and self.wells[meanWellDictionary[i]].contents.pH != "":
                        pdf.drawString(7*cm+(xpos % 2)*9.5*cm,22*cm - (ypos % yNum)*ySize*cm - offset*0.5*cm ,str(round(float(self.wells[meanWellDictionary[i]].contents.pH)+(self.wells[meanWellDictionary[i]].contents.dpH*(self.wells[meanWellDictionary[i]].Tm-20)),2)))
                        pdf.setFillColor("black")
                        if drawdpH ==False:
                            pdf.drawString(7*cm+(xpos % 2)*9.5*cm,22.5*cm - (ypos % yNum)*ySize*cm ,"Adjusted pH at Tm")
                            drawdpH = True
                    offset += 1
            drawdpH = False
            xpos +=1
            if newpage % 2 == 0:
                ypos +=1
            
            newpage += 1 
            plt.close()
            fig3 = plt.figure(num=1,figsize=(5,4))

       
        plt.close()
            
            
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        #save the pdf    
        pdf.save()
        return

def main():
    root = Tkinter.Tk()
    root.withdraw()
    tkMessageBox.showwarning("Inncorrect Usage", "Please run the 'RunMeltdown.bat' file from the same directory")
    return
    
    
if __name__ == "__main__":
    main()



















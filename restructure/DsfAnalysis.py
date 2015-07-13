# -*- coding: utf-8 -*-

import csv
import os
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
            lysozymeMeanWell = self.contentsHash[LYSOZYME]['']['']
            
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
        pdf.setFillColor("black")
        if suppliedProteinTm != None:
            pdf.drawString(15.5*cm,10.4*cm,"Protein as supplied") 
        
        pdf.setFont("Helvetica",12)
        fig1 = plt.figure(num=1,figsize=(10,8))
        
        # Maximum and minimum for the y axis of the summary graph
        maxi = 0
        mini = 100
        # Labels for the summary graph
        names = sorted(Contents.name, key=lambda x: x[1])
        labels = [x[0]+"("+str(x[1])+")" for x in names]
        tmHandles = []
        unreliableDrawn = False
        for i, saltConcentration in enumerate(Contents.salt):
            # Tms to be drawn regularly 
            tms = []
            # Tms to be drawn as unreliable
            badTms = []

            for condition in names:
                found = False
                for well in self.plate.names:
                    if self.wells[well].contents.salt == saltConcentration and self.wells[well].contents.name == condition[0] and self.wells[well].contents.pH == condition[1]:
                        if (self.wells[well].Tm != None and len(self.plate.meanDict[well]) > 1 and self.wells[well].TmError == None) or self.wells[well].complex == True:
                            tms.append(None)
                            badTms.append(self.wells[well].Tm)
                        else:
                            tms.append(self.wells[well].Tm)
                            badTms.append(None)
                        found = True
                        break
                # If the particular combination of buffer, salt and pH is not found
                if not found:
                    tms.append(None)
                    badTms.append(None)

            # Figuring out the scale of the y axis for the summary graph
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
            # Handle for the legend
            try:
                handle, = plt.plot([x for x in range(len(labels))],tms,color=COLOURS[i],marker="o",linestyle="None")
            except IndexError:
                tkMessageBox.showerror("Error", "Only up to 24 types of each condition are supported.\n(there is no limit to the number of conditions)\n\ne.g. 6 different salt concentrations per buffer")
                sys.exit(1)

            plt.plot([x for x in range(len(labels))],badTms,color=COLOURS[i],marker="d",linestyle="None")
            unreliableDrawn = False
            if badTms:
                unreliableDrawn = True
            tmHandles.append(handle)

        originalProteinMeanSd = rh.meanSd([self.originalPlate.wells[x].Tm for x in self.originalPlate.proteinAsSupplied if x not in self.delCurves])
        
        
        # Setting the scale of the y axis
        if originalProteinMeanSd[0]!= None:
            if(originalProteinMeanSd[0]-mini > maxi - originalProteinMeanSd[0]):
                plt.axis([-1,len(labels),mini-5,2*originalProteinMeanSd[0]-mini +10])
            else:
                plt.axis([-1,len(labels),2*originalProteinMeanSd[0] - maxi-5,maxi + 10])
        else:
            plt.axis([-1,len(labels),mini-5,maxi + 5])
        #padding at the top changes depending on how many different condition cvariable 2's there are (legend gets bigger)
        plt.gcf().subplots_adjust(bottom=0.35, top=0.85 - 0.035*(int(len(Contents.salt)/3)))
        plt.ylabel('Tm')

        # Drawing the Tm of the protein as supplied as a horizontal line on the summary graph
        plt.axhline(originalProteinMeanSd[0],0,1,linestyle="--",color="red")
        # Setting x axis labels
        plt.xticks([x for x in range(len(labels))],labels,rotation="vertical")
        # Putting the legend for the summary graph at the top, and show only 1 dot instead of 2
        plt.legend(tmHandles,Contents.salt,loc='lower center', bbox_to_anchor=(0.5, 1),ncol=3, fancybox=True, shadow=False, numpoints=1)
        # Saving the summary graph as an image and drawing it on the page
        imgdata = cStringIO.StringIO()
        fig1.savefig(imgdata, format='png',dpi=180)
        imgdata.seek(0)  # rewind the data
        Image = ImageReader(imgdata)
        pdf.drawImage(Image, cm, 4*cm, 16*cm, 11*cm)  
        plt.close()

        pdf.setFont("Helvetica",10)
        if unreliableDrawn:
            pdf.drawString(7.9*cm, 14.2*cm, "Tms drawn in diamonds may be unreliable")
        unreliableDrawn = False
        
        
        
        
        
        
        
        
        
        
        
        
        
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



















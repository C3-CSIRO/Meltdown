# -*- coding: utf-8 -*-
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
#
#  Authors: Nicholas Rosa, Marko Ristic, Shane A. Seabrook, David Lovell, Del Lucent and Janet Newman
#
#  Website: TBD.github.com
#
#  License: XXX

import Tkinter, tkFileDialog, tkMessageBox
import os
import ConfigParser
import sys, traceback

from DsfAnalysis import DsfAnalysis
from MeltdownException import MeltdownException
import meltdownReleases

#current version of meltdown, displayed in error logs
VERSION = "v2.0.0"

#the running location of this file
RUNNING_LOCATION = os.path.dirname(os.path.realpath(__file__))

#open the settings.ini and set the appropriate flags
cfg = ConfigParser.ConfigParser()
cfg.readfp(open(RUNNING_LOCATION + '/../settings.ini'))
#get options
DELETE_INPUT_FILES = cfg.getboolean('Running Options', 'DeleteInputFiles')
CREATE_NORMALISED_DATA = cfg.getboolean('Extra Output', 'ProduceNormalisedData')
CHECK_FOR_NEW_VERSION = cfg.getboolean('Running Options', 'CheckForNewVersion')

def main():
    #opens up selection windows for user to use
    root = Tkinter.Tk()
    root.withdraw()
    
    if (CHECK_FOR_NEW_VERSION):
        try:
            if (meltdownReleases.checkIfLatestRelease(VERSION) == 0):
                tkMessageBox.showinfo("New version availiable!","There is a newer version of Meltdown avaliable for download at\nhttps://github.com/C3-CSIRO/Meltdown")
        except:
            #if user doesn't have internet, or github is down, then errors will occur, and we skip the check
            pass
    
    try:
        #choosing a dsf results data file
        rfuFilepath = tkFileDialog.askopenfilename(title="Select the DSF Experiment Results", filetypes=[("text files", ".txt")])
        #raise error if dialog is closed before selecting a file
        if rfuFilepath == '':
            raise MeltdownException("Data file not selected")
        
        #choosing a contents map for the data file
        contentsMapFilepath = tkFileDialog.askopenfilename(title="Select the Contents Map", filetypes=[("text files", ".txt")])
        #raise error if dialog is closed
        if contentsMapFilepath == '':
            raise MeltdownException("Contents map file not selected")

        #the analysis
        print 'reading in data ...'
        #name the analysis the name of the data file
        experiment = DsfAnalysis(rfuFilepath.split('/')[-1])
        experiment.loadCurves(rfuFilepath,contentsMapFilepath)
        print 'analysing ...'
        experiment.analyseCurves()
        
        # generating the report
        print 'generating report ...'
        name = rfuFilepath.split(".")[0]
        experiment.generateReport(name+".pdf", VERSION)

        #remove any exported files in the directory of the data file. These files are identified if they
        #have the same word at the start of their file name, this is assumed to be the protein name, and
        #all files with the same first word in the directory are deleted
        if DELETE_INPUT_FILES:
            print 'deleting input data files ...'
            folder = rfuFilepath[:-len(rfuFilepath.split('/')[-1]) - 1]
            proteinName = rfuFilepath.split('/')[-1].split()[0]
            for fl in os.listdir(folder):
                if '.pdf' in fl:
                    continue
                if proteinName in fl:
                    os.remove(folder+'/'+fl)
                    
        #generate a tab delimited .txt file of the normalised curves
        if CREATE_NORMALISED_DATA:
            #add -normalised to the end of the filename
            print 'creating normalised data ...'
            experiment.produceNormalisedOutput(rfuFilepath[:-4] + '-normalised.txt')
        
        print '*done*'
            
    #expected error, to do with reading input, will give descriptive messages
    except MeltdownException as e:
        tkMessageBox.showerror("Error", e.message)
        print '*error occured* ' + e.message
    #unexpected errors only direct you to the error log
    except Exception:
        tkMessageBox.showerror("Error", "Unexpected error, please check the error log for more information")
        #save to the error log before finishing
        errors = open(RUNNING_LOCATION + "/../error_log.txt",'w')
        etype, value, tb = sys.exc_info()
        errors.write("Version: "+VERSION+"\n")
        errors.write(''.join(traceback.format_exception(etype, value, tb, None))) 
        root = Tkinter.Tk()
        root.withdraw()
        errors.close()
        print '*error occured* check error log'
    return

#excecutes main() on file run
if __name__ == "__main__":
    main()
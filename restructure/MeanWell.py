# -*- coding: utf-8 -*-

import Tkinter, tkMessageBox

class MeanWell:
    def __init__(self, tm, tmError, isComplex, replicates, contents):
        #relevant info for a mean well
        self.tm = tm
        self.tmError = tmError
        self.isComplex = isComplex
        self.replicates = replicates
        self.contents = contents
        return


def main():
    root = Tkinter.Tk()
    root.withdraw()
    tkMessageBox.showwarning("Inncorrect Usage", "Please run the 'RunMeltdown.bat' file from the same directory")
    return
    
    
if __name__ == "__main__":
    main()
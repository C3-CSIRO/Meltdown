# -*- coding: utf-8 -*-

import Tkinter, tkMessageBox

class Contents:
    def __init__(self, cv1, cv2, ph, dphdt, isControl, order):
        #stores the contents of a well
        self.cv1 = cv1
        self.cv2 = cv2
        self.ph = ph
        self.dphdt = dphdt
        self.isControl = isControl
        self.order = order
        return


def main():
    root = Tkinter.Tk()
    root.withdraw()
    tkMessageBox.showwarning("Inncorrect Usage", "Please run the 'RunMeltdown.bat' file from the same directory")
    return
    
    
if __name__ == "__main__":
    main()
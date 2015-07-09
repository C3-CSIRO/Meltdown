"""
Place this files into a directory containing csv files you
wish to convert into tab delimited text files, and run this file
all the files will now have converted txt versions saved in the same
directory

"""


import csv
import glob as g

for csvfile in g.glob('*.csv'):
    fileIn =  open(csvfile, 'rb')
    fileOut = open(csvfile[:-4] + '.txt', 'wb')
    reader = csv.reader(fileIn)
    writer = csv.writer(fileOut, delimiter='\t')
    writer.writerows(reader)
    
    fileIn.close()
    fileOut.close()
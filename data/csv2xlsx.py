"""
Useful converter of csvs to xlsx

Since batch .pcrd exports put all data in .csv files, running this within the 
same directory as all the .csv files will create a copy of every file as an
.xlsx file. In this form it will work with meltdown.py

Although this converts to .xlsx instead of .xls it will still work with
Meltdown as the files are created differently to how .pcrd files export them

"""

import os
import glob
import csv
from xlsxwriter.workbook import Workbook

for csvfile in glob.glob(os.path.join('.', '*.csv')):
    workbook = Workbook(csvfile + '.xlsx')
    worksheet = workbook.add_worksheet()
    with open(csvfile, 'rb') as f:
        reader = csv.reader(f)
        for r, row in enumerate(reader):
            for c, col in enumerate(row):
                try:
                    worksheet.write(r, c, float(col))
                except ValueError:
                    worksheet.write(r, c, col)

    workbook.close()
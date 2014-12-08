# csv to tab deliminated text
import os, pandas



folder_loc = "../UropCrystallisation/data/jcsg/rfuResults"

# Short piece of code for batch analysis of experiments
files = os.listdir(folder_loc)
total = len(files)
for i, bsc9 in enumerate(files):
    if ".csv" in bsc9:
        filepath = folder_loc + "/" + bsc9
        print filepath
        dataTxt = pandas.DataFrame.from_csv(filepath, index_col='Temperature')
        dataTxt.to_csv(path_or_buf= folder_loc + "/text/"+bsc9+".txt",sep="\t")
        print str(round((i+1)/float(total)* 100,2))  +"%"

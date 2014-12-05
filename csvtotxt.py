# csv to tab deliminated text
import os, pandas

# Short piece of code for batch analysis of experiments
files = os.listdir("../UropCrystallisation/data/bufferscreen9/rfuResults")
total = len(files)
for i, bsc9 in enumerate(files):
    filepath = "../UropCrystallisation/data/bufferscreen9/rfuResults/" + bsc9
    print filepath
    dataTxt = pandas.DataFrame.from_csv(filepath, index_col='Temperature')
    dataTxt.to_csv(path_or_buf= "../UropCrystallisation/data/bufferscreen9/rfuResults/text/"+bsc9+".txt",sep="\t")
    print str(round(i/float(total)* 100,2))  +"%"

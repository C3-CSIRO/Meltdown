import pandas as pd
data = pd.DataFrame.from_csv('data/AS_30194_buffer_screen_9_081113/Content_map.txt', sep='\t')
data2 = pd.DataFrame.from_csv('data/AS_30194_buffer_screen_9_081113/AS_30194_buffer_screen_9_081113 -  Melt Curve RFU Results_FRET.txt', sep='\t', index_col='Temperature')
print data
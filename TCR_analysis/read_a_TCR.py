import pandas as pd
from os import listdir
from os.path import isfile, join

#get all of the names of the files from a folder into a big long list
mypath = 'C:\Research_Assistant\work\data\TCR_data'
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
print('onlyfiles',onlyfiles)

for file in onlyfiles:
    original_filepath = join(mypath, file)
    df = pd.read_csv(original_filepath,sep='\t')
    print(file[:-7],file[:-7])
    final_filepath = file[:-7] + '.csv'
    df.to_csv(final_filepath, index=False)


# df = pd.read_csv('C:\Research_Assistant\work\data\TCR_data\dcr_HVO_439679_day7_1_alpha.tsv.gz',sep='\t')
# print(df.head())
# df.to_csv('dcr_HVO_439679_day7_1_alpha.csv', index=False)
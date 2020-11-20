import pandas as pd
import glob


path = '/home/linc/c00219805/Research/Protien_Database/'
subfolder= '/extracted_new_samples/sample_wu_t1/'
fileType = "*.keys_keycombine2"  
fileType_triplets = "*.triplets_theta21_dist12" 
setting = 'keycombine2'
files=glob.glob(path+subfolder+'/theta21_dist12/'+fileType) 
files_triplets = glob.glob(path+subfolder+'/theta21_dist12/'+fileType_triplets)


df_all = pd.read_table(path+subfolder+'/theta21_dist12/tripletes_keycombine2.csv', sep = ',', header = 0)
print df_all
df_key_range = df_all[['key','fileName','+/-2_new_key']]
df_key_occurences_perFile = pd.DataFrame({'count' : df.groupby( [ 'key','fileName','+/-2_new_key'] ).size()}).reset_index()
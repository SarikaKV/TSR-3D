import copy
import os,os.path
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import normalize
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import cosine_similarity
from scipy.spatial.distance import pdist,squareform

class VectorizationFSPCA:

    def __init__(self, **kwargs):
        self.setting = kwargs["setting"] 
        self.sample_dict = kwargs["sample_dict"]
        self.f1_in = os.path.join(kwargs["outFolder"], 'localFeatureVect{}_before_fs.csv'.format(self.setting))
        self.feat_out = os.path.join(kwargs["outFolder"], 'selected_feature_names.csv')

        self.cosine = os.path.join(kwargs["outFolder"], 'cosine.csv')
        self.euclidean = os.path.join(kwargs["outFolder"], 'euclidean.csv')


    def vectorize(self):
        arrs = []   
        filenames = []   
        distance_filenames = []
        with open(self.f1_in) as fcsv:
            lines=fcsv.readlines()
            for idx,line in enumerate(lines):
                filenames.append(str(line.split(';')[0]).upper() )
                distance_filenames.append(str(self.sample_dict[str(line.split(';')[0]).upper()]) \
                        + '-' + str(line.split(';')[0]).upper())
                print(str(self.sample_dict[str(line.split(';')[0]).upper()]) \
                        + '-' + str(line.split(';')[0]).upper() )
                l = list(line.split(';')[1].split(','))
                l_arr = np.asarray(l[:-1]).astype(np.float) 
                arrs.append(l_arr)
        #print(data)
        x = np.array(arrs)
        print(x)
        #Scale
        # x = StandardScaler().fit_transform(x)
        # print(x)

        #Normalize
        x = normalize(x, axis=1, norm='l1')
        print(x)
        print(np.sum(x[0]))
        reduced_x = x

        pca = PCA(0.9999)
        pca.fit_transform(x)
        print(pca.n_components_)
        reduced_x = pca.transform(x)
        print(reduced_x)
        print(reduced_x.shape)

        df = pd.DataFrame(reduced_x, index = filenames)
        df.to_csv(self.feat_out)
        
        #print(squareform(pdist(df, metric='cosine')))
        pd.DataFrame(squareform(pdist(df.values, metric='cosine')), columns = distance_filenames, index = distance_filenames).to_csv(self.cosine)
        #pd.DataFrame(cosine_similarity(df)).tocsv()
            	

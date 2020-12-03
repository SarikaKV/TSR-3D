import copy
import os,os.path
import pandas as pd
import scipy.sparse
from scipy.sparse import csr_matrix

__author__ = "Venkata Sarika Kondra"

__version__ = "1.0.1"
__maintainer__ = "Venkata Sarika Kondra"
__email__ = "c00219805@louisiana.edu"


class VectorizationNoFS:

    def __init__(self, **kwargs):
        self.setting = kwargs["setting"] 
        self.fileList = kwargs["filesList"]
        self.f1_out = os.path.join(kwargs["outFolder"], 'feature_map_nofs{}.npz'.format(self.setting))
        self.feat_names = os.path.join(kwargs["outFolder"], 'feature_names{}.csv'.format(self.setting))

    def vectorize(self):
        keyset = {}
        rows = []
        columns = []
        data = []
        doc_count = 0        
        for i in self.fileList:
            for line in open(i,'r'):
                index = keyset.setdefault(line.split()[0].rstrip(), len(keyset))
                rows.append(doc_count)
                columns.append(index)
                data.append(int(line.split()[1].rstrip()))
            doc_count += 1
        #print(data)
        sparse_matrix = csr_matrix((data, (rows, columns)), dtype = int, shape=(len(self.fileList), len(keyset)))
        scipy.sparse.save_npz(self.f1_out, sparse_matrix)
        pd.DataFrame(keyset.items(), columns = ['key_no','index']).to_csv(self.feat_names)
            
            	

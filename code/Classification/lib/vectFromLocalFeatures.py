import copy
import os,os.path

__author__ = "Venkata Sarika Kondra"

__version__ = "1.0.1"
__maintainer__ = "Venkata Sarika Kondra"
__email__ = "c00219805@louisiana.edu"


class Vectorization:

    def __init__(self, **kwargs):
        self.setting = kwargs["setting"] 
        self.outFolder = kwargs['outFolder']
        self.fileList = kwargs["filesList"]
        self.f1_in = open(kwargs["outFolder"] + '//localFeatureSelection' + self.setting + '.txt', 'r')
        self.f1_out = open(kwargs["outFolder"] + '//localFeatureVect' + self.setting + '.csv', 'w')
        self.f2_out = open(kwargs["outFolder"] + '//feature_names' + self.setting + '.csv', 'w')        
        self.add_header = kwargs["add_header"]
        self.details_dict = kwargs["details_dict"]
        if self.add_header:
            self.fm_out = open(kwargs["outFolder"] + '//feature_map_with_header.csv', 'w')

    def vectorize(self):
        keyDict = {}
        for i in self.f1_in:
            i = i.split()
            keyDict[i[0].rstrip()] = 0
            
        self.f1_in.close()
        
        file_count = 0
        for i in self.fileList:
            filename = open(self.outFolder + '//{}.feature_names'.format(i.split('.')[0]), 'w')
            f2_in = open(i,'r')
            keyDict1 = copy.deepcopy(keyDict)
            for j in f2_in:
                j = j.split()
                if j[0].rstrip() in keyDict1:
                    keyDict1[j[0].rstrip()] = j[1].rstrip()
            f2_in.close()
            # print('dict: ', len(keyDict.keys()))
            # print('dict1: ', len(keyDict1.keys()))
            
            self.f1_out.writelines([str(i.split('.')[0]), ';'])
            if self.add_header:
                if file_count == 0:
                    self.fm_out.writelines(",".join(['group', 'protein'] + list(keyDict1.keys())))
                    self.fm_out.writelines("\n")
                self.fm_out.writelines([str(self.details_dict[(i.split('.')[0]).upper()]),',', str(i.split('.')[0]), ','])
            for k in keyDict1:
                filename.writelines([k, '\n'])
                self.f1_out.writelines([str(keyDict1[k]), ','])
                if self.add_header:
                    self.fm_out.writelines([str(keyDict1[k]), ','])
            self.f1_out.writelines(['0', '\n'])
            if self.add_header:
                self.fm_out.writelines(['0', '\n'])
            file_count += 1
            filename.close()
        self.f1_out.close()
        if self.add_header:
            self.fm_out.close()
        
            	

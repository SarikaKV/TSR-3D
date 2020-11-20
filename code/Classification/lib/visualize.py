#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import scipy,csv, argparse, os, time
import scipy.cluster.hierarchy as hac
from scipy.cluster.hierarchy import ward, average,dendrogram, linkage, fcluster, cut_tree
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics import homogeneity_score, completeness_score
import scipy.spatial as sp
from sklearn.manifold import MDS
from collections import Counter

from mpl_toolkits.mplot3d import Axes3D
#import matplotlib as mpl
#mpl.use('Agg')
import seaborn as sns
sns.set_palette("Greens")
import matplotlib.pyplot as plt

from os.path import expanduser

class Visualize:

    def __init__(self,**kwargs):
        self.setting = kwargs["setting"]
        self.folder = kwargs["outFolder"]
        self.writer = pd.ExcelWriter(os.path.join(self.folder, 'similarity_values.xlsx'), \
                engine='xlsxwriter')
        self.has_columns =     kwargs["has_columns"]
        self.details_dict = kwargs["details_dict"]

    def modify_file(self, file):
        """Replace Feng's output contents of semicolon to comma."""
        text = open(file, "r")
        text = ''.join([i for i in text]) \
                .replace(";", ",")
        f = open(file, "w")
        f.writelines(text)
        f.close()

    def fancy_dendrogram(self, *args, **kwargs):
        max_d = kwargs.pop('max_d', None)
        if max_d and 'color_threshold' not in kwargs:
            kwargs['color_threshold'] = max_d
        annotate_above = kwargs.pop('annotate_above', 0)

        ddata = dendrogram(*args, **kwargs)

        if not kwargs.get('no_plot', False):
            plt.title('Hierarchical Clustering Dendrogram (truncated)')
            plt.xlabel('sample index or (cluster size)')
            plt.ylabel('distance')
            for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
                x = 0.5 * sum(i[1:3])
                y = d[1]
                if y > annotate_above:
                    plt.plot(x, y, 'o', c=c)
                    plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                                 textcoords='offset points',
                                 va='top', ha='center')
            if max_d:
                plt.axhline(y=max_d, c='k')
        return ddata
    def get_all_figures(self, type, cluster_method, similarity_x, x, fileList, min_cutoff, max_cutoff, step_size):
        """Generate heatmap, clustermap and dendrogram for the similarity matrix."""
        # clus2_dict = {'chymotrypsin' : 11, 'elastase' : 12, 'plasmin' : 13,'prothrombin' : 14, 'subtilisin' : 15, 'trypsin' : 16, 
        # 'src' : 17, 'cdk' : 18, 'SRC' : 17, 'CDK' : 18, 'tyr_phosphatase' : 19, 'Tyr Phosphatase' : 19, 'Tyr_Phosphatase' : 19, 'ser_thr_phosphatase' : 20, 
        # 'Ser/Thr Phosphatase' : 20, 'Ser_Thr_Phosphatase' : 20, 'elastase_1' : 21, 'elastase_2' : 22, 'digestive' : 23, 'kinase' : 24, 'regulatory' : 25, 
        # 'phosphatase' : 26, 'Phosphatase' : 26, 'protease' : 27, 'kinase_phosphatase' : 28, 'STAT' : 29, 'stat' : 29, 'CBP KIX' : 30, 'cbp kix' : 30}


        clus2_dict = {v: k for k,v in enumerate(set(self.details_dict.values()))}


        df_evid = pd.DataFrame(list(clus2_dict.items()),columns = ['group_name','group_no'])

        plt.gcf().clear()


        df_hm = pd.DataFrame(similarity_x)
        df_hm = df_hm[sorted(list(df_hm))]
        df_hm = df_hm.sort_index()
        ax_heatmap = sns.heatmap(df_hm,  
                square = True, cbar_kws={ "shrink": 0.5}) #Plot the correlation as heat map
        fig = ax_heatmap.get_figure()
        fig.savefig(os.path.join(self.folder, 'heatmap_{}.png'.format(type)))
            
        df = pd.DataFrame(similarity_x)
        #print('df',(100-df).lt(0).sum().sum())
        linkage = hac.linkage(sp.distance.squareform(100 - df, checks = False), method=cluster_method)
        ax_clustermap = sns.clustermap(100 - df, row_linkage=linkage, col_linkage=linkage)
        ax_clustermap.savefig(os.path.join(self.folder, 'clustermap_{}.png'.format(type)))
        pd.DataFrame(list(df.index[[ax_clustermap.dendrogram_row.reordered_ind]]))\
                .to_csv(os.path.join(self.folder,"clustermap_reordered_names_{}.csv".format(type)))
            
        corr_condensed = hac.distance.squareform(x.values, checks = False)
        linkage_matrix = hac.linkage(corr_condensed, method=cluster_method)
        fig, ax = plt.subplots(figsize=(15, 20)) # set size
        ax = dendrogram(linkage_matrix, orientation="left", labels=fileList) 
         
        
        plt.tick_params(\
                axis= 'x',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                labelsize=11,
                )
        plt.tick_params(\
                axis= 'y',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                #labelsize=11,
                )
        plt.title("{} Jaccard Dendogram".format(type.capitalize()) )
        

        best_cutoff, best_ari = self.get_cutoff_plots(type, linkage_matrix, df_evid, fileList, min_cutoff, max_cutoff, step_size, cluster_method)
        #best_cutoff = 13
        plt.axvline(c = 'g', linestyle='--', x=best_cutoff)
        plt.savefig(os.path.join(self.folder,'dendo_{}_{}_{}.png'.format(cluster_method, type, best_cutoff)))
        #plt.show()

        plt.gcf().clear()
        

        ####################Change the cluster numbers according to dataset##################
        #Using cluster numbers
        cluster_no = len(set(self.details_dict.values())) #
        clusters = fcluster(linkage_matrix, cluster_no, criterion='maxclust')

        #Using distance cutoff
        # cluster_no = best_cutoff #(it is max_d)
        # clusters = fcluster(linkage_matrix, cluster_no, criterion='distance')

        print("Unique clusters: {}".format(len(set(clusters))))
        
        orig_groups = [x.split('-')[0] for x in fileList]
        df_clust = pd.DataFrame(list(zip(fileList, clusters, orig_groups)), columns = ['protein', 'group', 'group_name'])    
        #print(df_clust)    
        
        df_results = df_clust.join(df_evid.set_index('group_name'), on='group_name')
        #print(df_results)
        print("ARI for {} clusters of {} in data: {}".format(len(set(self.details_dict.values())), 
            set(self.details_dict.values()), 
            adjusted_rand_score(df_results['group'].values, df_results['group_no'].values)))
        print("Homogenity : {}".format(
            homogeneity_score(df_results['group'].values, df_results['group_no'].values)))
        print("Completeness: {}".format( 
            completeness_score(df_results['group'].values, df_results['group_no'].values)))
        df_results.to_csv(os.path.join(self.folder,"dendro_clusters_{}_{}_{}.csv".format(cluster_method, type, cluster_no)))
        #print(type, silhouette_score(corr_condensed, df_results['group'].values))
        print(cluster_method)
        

    def get_cutoff_plots(self, type, linkage_matrix, df_evid, fileList, min_cutoff, max_cutoff, step_size, cluster_method):
        cutoffs = []
        aris = []
        #print(np.arange(float(min_cutoff), float(max_cutoff), step_size))
        for i in np.arange(float(min_cutoff), float(max_cutoff), step_size):
            cluster_no = i #(it is max_d)
            clusters = fcluster(linkage_matrix, cluster_no, criterion='distance')
            orig_groups = [x.split('-')[0] for x in fileList]
            df_clust = pd.DataFrame(list(zip(fileList, clusters, orig_groups)), columns = ['protein', 'group', 'group_name'])
            df_results = df_clust.join(df_evid.set_index('group_name'), on='group_name')
            ari = round(adjusted_rand_score(df_results['group'].values, df_results['group_no'].values), 2)

            cutoffs.append(round(i,2))
            aris.append(ari)
        #     #print(i, ari)
        # fig = plt.figure(figsize=(12, 12))
        # ax = fig.add_subplot(111)
        # plt.plot(cutoffs, aris, 'b')
        # prev_y = 0
        # for x, y in zip(cutoffs, aris):
        #     if prev_y != y:
        #         ax.annotate(str(y), xy=(x, y), textcoords='data', fontsize = 8)
        #     prev_y = y
        # plt.xlabel("Cut-off")
        # plt.ylabel("ARI")
        # plt.title("{}: Cutoffs Vs. ARIs".format(type.capitalize()) )
        # #plt.show()
        # plt.savefig(os.path.join(self.folder,'cutoff_plot_{}_{}.png'.format(cluster_method, type)))
        #plt.gcf().clear()
        d = dict(zip( aris, cutoffs))
        pd.DataFrame(dict(zip(cutoffs, aris)).items(), columns = ['Cutoff', 'ARIs']).sort_values(['Cutoff']).to_csv(
            os.path.join(self.folder,'cutoff_data_{}_{}.csv'.format(cluster_method, type)))
        print(type, 'Best cutoff: ', d[np.max(aris)], 'Best ARI: ', np.max(aris))
        return d[np.max(aris)], np.max(aris)

    def visualize(self, sim_type, cluster_method, min_cutoff, max_cutoff, step_size):
        if self.has_columns:
            # print(os.path.join(self.folder,\
            #     "{}.csv".format(sim_type)))
            x = pd.read_csv(os.path.join(self.folder,\
                "{}.csv".format(sim_type)),\
                header = 0,index_col=0)    
            x = x.dropna()
            column_names =      x.columns
            column_names = [str.replace(protien, 'beta-barrel', 'beta_barrel') for protien in column_names]
            column_names = [str.replace(protien, 'cdk8-cycc', 'cdk8_cycc') for protien in column_names]
            column_names = ["{}-{}".format(str(self.details_dict[protien.split('-')[1]]),protien.split('-')[1]) for protien in column_names]
            x.columns = column_names
        else:
            #Replace semicolon in Feng's output to comma        
            self.modify_file(os.path.join(self.folder,"{}.csv".format(sim_type)))

            x = pd.read_csv(os.path.join(self.folder,"{}.csv".format(sim_type)),\
                 header = None, index_col = 0)
            x = x.dropna()
            #print(x.shape)
            column_names =     x.index
            column_names = ["{}-{}".format(str(self.details_dict[protien]),protien) for protien in column_names]
            print(column_names)
            x.columns = column_names
        x.index = column_names
        #print('x-negative',x.lt(0).sum().sum())
        max_val = x.max().max()
        #print('before', x.head(5))
        #print('afer', (x/max_val).head(5))
        if max_val > 1:
            x = x/max_val
        x = x.astype(float)
        similarity_x = (1-x) *100
        similarity_x[similarity_x < 0] = 0
        #print('sim-negative',similarity_x.lt(0).sum().sum(), \
        #'(np)sim-negative',np.sum((similarity_x < 0).values.ravel()))
        similarity_x.to_excel(self.writer,sheet_name=sim_type)
        
        self.get_all_figures(sim_type, cluster_method, similarity_x, x, x.columns, min_cutoff, max_cutoff , step_size)
        return x

    def visualize_all(self):
        start_time=time.time()                

        # All Maps
        # Comment any of the below if you have a Memory Error (Sarika)
        #average: {generalised: 0.5, ce_rmsd: 1.42, tm_rmsd: 1.43}
        #ward: {generalised: 0.61, ce_rmsd: 2.5, tm_rmsd: 3.3}

        # self.visualize('normal', 'ward',0, 3, 0.05)
        # self.visualize('wu', 'ward',0, 3, 0.05)
        # self.visualize('sarika', 'ward',0, 1.5, 0.05)
        # self.visualize('generalised', 'ward',0, 4, 0.1)

        # self.visualize('ce_rmsd', 'ward',0, 40, 0.125)
        # self.visualize('tm_rmsd', 'ward',0,40, 0.125)

        # self.visualize('ce_zscore_updated', 'ward',0, 40, 0.01)
        # self.visualize('tm_zscore_updated', 'ward',0, 4, 0.01)
        # self.visualize('dali_zscore_updated', 'ward', 0, 300, 0.5)
        # self.visualize('sequence_rmsd_updated', 'ward', 0, 14, 0.1)

        #average
        # self.visualize('normal', 'average',0, 1, 0.01)
        # self.visualize('wu', 'average',0, 1, 0.01)
        # self.visualize('sarika', 'average',0, 1, 0.001)
        # self.visualize('generalised', 'average',0, 0.7, 0.005)

        # self.visualize('ce_rmsd', 'average',0, 7, 0.01)
        # self.visualize('tm_rmsd', 'average',0,7, 0.01)

        # self.visualize('ce_zscore_updated', 'average',0, 40, 0.005)
        # self.visualize('tm_zscore_updated', 'average',0, 4, 0.001)
        # self.visualize('dali_zscore_updated', 'average', 0, 300, 0.5)
        # self.visualize('sequence_rmsd_updated', 'average', 0, 3, 0.01)
        ########################################################################

        # self.visualize('normal', 'ward',0, 1, 0.01)
        # self.visualize('wu', 'ward',0, 1, 0.01)
        # self.visualize('sarika', 'ward',0, 1, 0.01)
        self.visualize('generalised', 'ward',0, 1, 0.01)
        
        # self.visualize('tm_rmsd', 'ward',0,1, 0.01)
        # self.visualize('ce_rmsd', 'ward',0, 1, 0.01)        
        
        # self.visualize('tm_zscore_updated', 'ward',0, 1, 0.01)
        # self.visualize('ce_zscore_updated', 'ward',0, 1, 0.01)
        # self.visualize('dali_zscore_updated', 'ward', 0, 1, 0.01)
        # self.visualize('sequence_rmsd_updated', 'ward', 0, 1, 0.01)
        print('*'*200)

        #average
        # self.visualize('normal', 'average',0, 1, 0.01)
        # self.visualize('wu', 'average',0, 1, 0.01)
        # self.visualize('sarika', 'average',0, 1, 0.001)
        self.visualize('generalised', 'average',0, 1, 0.01)
        
        # self.visualize('tm_rmsd', 'average',0, 1, 0.01)
        # self.visualize('ce_rmsd', 'average',0, 1, 0.01)
        
        # self.visualize('tm_zscore_updated', 'average',0, 1, 0.001)
        # self.visualize('ce_zscore_updated', 'average',0, 1, 0.01)        
        # self.visualize('dali_zscore_updated', 'average',0, 1, 0.01)
        # self.visualize('sequence_rmsd_updated', 'average', 0, 1, 0.01)


        self.writer.close()
        end_time=time.time()
        total_time=((end_time)-(start_time))
        print("Time taken for visualization: {}".format(total_time))

        

#!/usr/bin/env python
#all_angle_lengths.py

""" Calculates all angles, lengths and representative angles and lengths.

    For a given triple of amino acids, this calculates all the angles , lengths 
    and representative angle and length. The angles are the angles formed by the 
    median of corresponding vertex and opposite edge.
    

    Attributes:
        path: Location where input required resources are stored.
        subfolder: The sample for which calculations are required.
        aminoAcidCode: Lexicographic file required for rule-based assignment.
        
"""

import math,glob,os,time
from collections import Counter

import pandas as pd
from joblib import Parallel, delayed, cpu_count


__author__ = "Sumi Singh, Venkata Sarika Kondra"

__version__ = "1.0.2"
__maintainer__ = "Venkata Sarika Kondra"
__email__ = "c00219805@louisiana.edu"


round_off_to = 2
total_samples = 12
setting = 'corrected'
def calcDist(indexLabel1,indexLabel2):
    """Calculate Distance between two points in 3D space."""
    x1=xCord[indexLabel1]
    x2=xCord[indexLabel2]
    y1=yCord[indexLabel1]
    y2=yCord[indexLabel2]
    z1=zCord[indexLabel1]
    z2=zCord[indexLabel2]
    distance=(((x1-x2)**2+(y2-y1)**2+(z2-z1)**2)**0.5)
    return distance
def indexFind(index_of_2,i1,j1,k1):
    if index_of_2==i1:
        indexOf0=j1
        indexOf1=k1
    elif index_of_2==j1:
        indexOf0=i1
        indexOf1=k1
    elif index_of_2==k1:
        indexOf0=i1
        indexOf1=j1
    return indexOf0, indexOf1
def processFiles(fileName):
    """Calculates all angles, all lengths, representative angle and maxDist after performing rule-based labelling.

    Arguments:
        fileName: The protein file in PDB/ENT format.

    Returns:
        all_angleList: A Counter having all angles formed by their medians on opposite edges of the 
        non-collinear triangle formed by the three amino acids at i, j and k 
                    and their frequencies of occurences in this protein file rounded to next significant digit.
        rep_angleList: A Counter having representative angle and its frequency
        all_lengthsList: A counter having lengths of all edges of the non-collinear triangle formed.
        maxDist: Maximum length among all lengths calculated above.

    """
    print fileName
    count_t1 = 0
    inFile=open(fileName,'r')
    all_angleList = Counter()
    rep_angleList = Counter()
    all_lengthsList = Counter()
    maxDist_List = Counter()
    global xCord, yCord, zCord
    aminoAcidName={}
    xCord={}
    yCord={}
    zCord={}
    seq_number={}
    counter=0
    for i in inFile:
        if (i[0:6].rstrip()=="NUMMDL"):
            numOfModels=i[10:14].rstrip()
        if ((i[0:6].rstrip()=="ENDMDL")or (i[0:6].rstrip()=='TER')):
            break
        if (i[0:6].rstrip()=="MODEL" and int(i[10:14].rstrip())>1):
            break
            
        if(i[0:4].rstrip())=="ATOM" and(i[13:15].rstrip())=="CA" and(i[16]=='A'or i[16]==' ')and i[17:20]!= "UNK" :
            aminoAcidName[counter]=int(aminoAcidLabel[i[17:20]])
            xCord[counter]=(float(i[30:38]))
            yCord[counter]=(float(i[38:46]))
            zCord[counter]=(float(i[46:54]))
            seq_number[counter]=str(i[22:27])
            counter+=1

    protLen=len(yCord)
    initialLabel=[]
    sortedLabel=[]
    sortedIndex=[]
    outDist={}
    for m in range(0,3):
        initialLabel.append(0)
        sortedLabel.append(0)
        sortedIndex.append(0)

    for i in range(0,protLen-2):
        for j in range(i+1,protLen-1):
            for k in range(j+1, protLen):
                global i1,j1,k1
                i1=i
                j1=j
                k1=k
                keepLabelIndex={}
                keepLabelIndex[aminoAcidName[i]]=i
                keepLabelIndex[aminoAcidName[j]]=j
                keepLabelIndex[aminoAcidName[k]]=k
                initialLabel[0]=aminoAcidName[i]
                initialLabel[1]=aminoAcidName[j]
                initialLabel[2]=aminoAcidName[k]
                sortedLabel=list(initialLabel)
                sortedLabel.sort(reverse=True)

                #Perform Rule- based labelling

                if (sortedLabel[0]==sortedLabel[1])and(sortedLabel[1]==sortedLabel[2]):
                    dist1_2Temp=calcDist(i,j)
                    dist1_3Temp=calcDist(i,k)
                    dist2_3Temp=calcDist(j,k)
                    if dist1_2Temp>=(max(dist1_2Temp,dist1_3Temp,dist2_3Temp)):
                        indexOf0=i
                        indexOf1=j
                        indexOf2=k
                    elif dist1_3Temp>=(max(dist1_2Temp,dist1_3Temp,dist2_3Temp)):
                        indexOf0=i
                        indexOf1=k
                        indexOf2=j
                    else:
                        indexOf0=j
                        indexOf1=k
                        indexOf2=i
                elif(aminoAcidName[i]!=aminoAcidName[j])and(aminoAcidName[i]!=aminoAcidName[k]) and(aminoAcidName[j]!=aminoAcidName[k]): 
                    for index_ in range(0,3):
                        sortedIndex[index_]=keepLabelIndex[sortedLabel[index_]]
                    indexOf0=sortedIndex[0]
                    indexOf1=sortedIndex[1]
                    indexOf2=sortedIndex[2]
                elif(sortedLabel[0]==sortedLabel[1])and(sortedLabel[1]!=sortedLabel[2]):
                    indexOf2=keepLabelIndex[sortedLabel[2]]
                    indices=indexFind(indexOf2,i,j,k)
                    a=indexOf2
                    b=indices[0]
                    c=indices[1]
                    dist1_3Temp=calcDist(b,a)
                    dist2_3Temp=calcDist(c,a)
                    if dist1_3Temp>=dist2_3Temp:
                        indexOf0=indices[0]
                        indexOf1=indices[1] 
                    else:
                        indexOf0=indices[1]
                        indexOf1=indices[0]
                elif(sortedLabel[0]!=sortedLabel[1])and(sortedLabel[1]==sortedLabel[2]):
                    indexOf0=keepLabelIndex[sortedLabel[0]]
                    indices=indexFind(indexOf0,i,j,k)
                    if calcDist(indexOf0,indices[0])>= calcDist(indexOf0,indices[1]):
                        indexOf1=indices[0]
                        indexOf2=indices[1] 
                    else:
                        indexOf2=indices[0]
                        indexOf1=indices[1]
                dist01=calcDist(indexOf0,indexOf1)
                s2=dist01/2
                dist02=calcDist(indexOf0,indexOf2)
                s1=dist02
                dist12=dist01
                dist03=calcDist(indexOf1,indexOf2)

                # All lengths calculation 
                all_lengthsList[round(dist01,round_off_to)] += 1
                all_lengthsList[round(dist02,round_off_to)] += 1
                all_lengthsList[round(dist03,round_off_to)] += 1

                maxDist_List[round(max(dist01,dist02,dist03),round_off_to)] +=1

                s3=(((xCord[indexOf0]+xCord[indexOf1])/2-xCord[indexOf2])**2
                    +((yCord[indexOf0]+yCord[indexOf1])/2-yCord[indexOf2])**2
                    +((zCord[indexOf0]+zCord[indexOf1])/2-zCord[indexOf2])**2)**0.5
                
                
                Theta1=180*(math.acos((s1**2-s2**2-s3**2)/(2*s2*s3)))/3.14
                if Theta1<=90:
                   all_angleList[round(Theta1,round_off_to)] +=1
                   rep_angleList[round(Theta1,round_off_to)] +=1
                else:
                    all_angleList[round(abs(180-Theta1),round_off_to)] +=1
                    rep_angleList[round(abs(180-Theta1),round_off_to)] +=1
                    
                #if Theta1>90:   
                 #   Theta1=abs(180-Theta1)
                #print 'Second Theta1, ',Theta1
                #Theta 2
                dist02=calcDist(indexOf1,indexOf0)
                s1=dist02
                dist01=calcDist(indexOf1,indexOf2)
                s2=dist01/2
                s3=(((xCord[indexOf1]+xCord[indexOf2])/2-xCord[indexOf0])**2
                    +((yCord[indexOf1]+yCord[indexOf2])/2-yCord[indexOf0])**2
                    +((zCord[indexOf1]+zCord[indexOf2])/2-zCord[indexOf0])**2)**0.5
                         
                Theta2=180*(math.acos((s1**2-s2**2-s3**2)/(2*s2*s3)))/3.14  
                #if Theta2 > 90:
                 #    Theta2 = abs(180-Theta2)
                if Theta2<=90:
                   all_angleList[round(Theta2,round_off_to)] +=1
                else:
                    all_angleList[round(abs(180-Theta2),round_off_to)] +=1

                #Theta 3
                dist02=calcDist(indexOf2,indexOf1)
                s1=dist02
                dist01=calcDist(indexOf2,indexOf0)
                s2=dist01/2
                s3=(((xCord[indexOf2]+xCord[indexOf0])/2-xCord[indexOf1])**2+
                    ((yCord[indexOf2]+yCord[indexOf0])/2-yCord[indexOf1])**2+
                    ((zCord[indexOf2]+zCord[indexOf0])/2-zCord[indexOf1])**2)**0.5
                
                Theta3=180*(math.acos((s1**2-s2**2-s3**2)/(2*s2*s3)))/3.14 
                #if Theta3 > 90:
                 #    Theta3 = abs(180-Theta3)
                if Theta3<=90:
                   all_angleList[round(Theta3,round_off_to)] +=1
                else:
                    all_angleList[round(abs(180-Theta3),round_off_to)] +=1
                # Either writting output to a file or using dictionary or 
                # counter will save you from memory exceptions in this case.
                #all_angleList[round(Theta1,round_off_to)] +=1
                #all_angleList[round(Theta2,round_off_to)] +=1
                #all_angleList[round(Theta3,round_off_to)] +=1

                #rep_angleList[round(Theta1,round_off_to)] +=1

                count_t1 = count_t1+1

    print 'count_t1:',count_t1

    return [all_angleList,rep_angleList,all_lengthsList,maxDist_List]

for i in range(1,total_samples +1):
    #Path where the sample PDB files present.
    path = '/home/linc/c00219805/Research/Protien_Database/extracted_new_samples/'

    subfolder="sample"+str(i) +"//"   ## Change sample no here for generatig results for different samples 
    #subfolder="sample_4t1//"
    print subfolder
    os.chdir(path)

    files=glob.glob(path+subfolder+"*.pdb")#Change file extension here if not pdb. Others may be .ent

    aminoAcidCode=open('/home/linc/c00219805/Research/Protien_Database/'+"aminoAcidCode_lexicographic _new.txt","r")
    #Start the timer
    start_time=time.time()
    aminoAcidLabel={}
    for amino in aminoAcidCode:
        amino=amino.split()
        aminoAcidLabel[amino[0]]=int(amino[1])
    aminoAcidCode.close()
        
    alltheta = Counter()
    all_repAngle = Counter()
    all_length = Counter()
    maxLength = Counter()
    #Parallel processing to use the power of CPU. This uses 2 cores less than CPU cores for parallelization
    a = Parallel(n_jobs=cpu_count() - 2, verbose=10, backend="multiprocessing", batch_size="auto")(delayed(processFiles)(fileName) for fileName in files)
    for t in a:
        alltheta += t[0]
        all_repAngle += t[1]
        all_length += t[2]
        maxLength += t[3]


    #Uncomment this for sequential processing and add similar steps for representative angle and lengths
    # for fileName in files:
    #     alltheta +=processFiles(fileName) 

    #For All thetas
    df_counts = pd.DataFrame.from_dict(alltheta, orient='index').reset_index()
    df_counts = df_counts.rename(columns={'index':'theta', 0:'freq'})
    print df_counts  
    df_counts.to_csv(path+subfolder+'//all_angles'+str(round_off_to)+setting+'.csv', sep=',')  

    #For representative theta
    df_counts = pd.DataFrame.from_dict(all_repAngle, orient='index').reset_index()
    df_counts = df_counts.rename(columns={'index':'theta', 0:'freq'})
    df_counts.to_csv(path+subfolder+'//theta'+str(round_off_to)+setting+'.csv', sep=',')  

    #For All Lengths
    df_counts = pd.DataFrame.from_dict(all_length, orient='index').reset_index()
    df_counts = df_counts.rename(columns={'index':'dist123', 0:'freq'})
    print df_counts  
    df_counts.to_csv(path+subfolder+'//all_length'+str(round_off_to)+setting+'.csv', sep=',')  

    #For max Dist
    df_counts = pd.DataFrame.from_dict(maxLength, orient='index').reset_index()
    df_counts = df_counts.rename(columns={'index':'maxDist', 0:'freq'})
    df_counts.to_csv(path+subfolder+'//maxDist'+str(round_off_to)+setting+'.csv', sep=',')  

    #End timer and calculate total time taken
    end_time=time.time()
    total_time=((end_time)-(start_time))
    print ("Code End Angle & Length calculation.")
    print ("TOTAL TIME IN MIN=",round(total_time/60,0))

                
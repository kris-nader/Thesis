# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 17:50:49 2021

@author: Kristen Michelle Nader
Hill Climbing trial ~~ i dont see the stochasticity here so??
"""

aa=["A","R","N","D","C","E","Q","G","H",
    "I","L","K","M","F","P","S","T","W","Y","V"]

### START W HSP70
## Starts the same create the same individual list
## Starts with the alanine seq

import pandas as pd
## mutate 1st position
def create_indlist(X,mutate_position):  
    with open('individual_list.txt','w') as myfile:
        for i in aa:
            myfile.write("A"+X+str(mutate_position)+i+";\n")
    myfile.close()
   
def readSummary(location_file,mutate_position):
    #data=pd.read_csv(location_file, sep="\t")
    data=pd.read_csv(location_file, sep="\t")
    data.columns=["Pdb", "Group1", "Group2","IntraclashesGroup1","IntraclashesGroup2","Interaction_Energy","StabilityGroup1","StabilityGroup2"]
  #  individual_list=pd.read_csv("C:\\Users\\user\\Desktop\\hill_climbing\\position2\\individual_list.txt", sep="\t",header=None)
    temp=[]
    for i in range(len(data)) : 
      temp.append(int(data.iloc[i, 0].split("_")[-1:][0].split(".")[0]))
    data["split"]=temp
    data=data.sort_values(by='split')
    data=data[['Pdb','Interaction_Energy']]
    #data["aa_mod"]=list(individual_list[0])
    return (data)
    
def find_best_binders(data):
    keep_next_round=[]
    index=[]
    file=data.Pdb
    interaction_energy=data.Interaction_Energy
    fitness_base=max(data.Interaction_Energy)
    for i in range(len(file)):
        if interaction_energy[i]<fitness_base:
            index.append(int(file[i].split("_")[-1:][0].split(".")[0]))
            keep_next_round.append(file[i])
    temp_df=pd.DataFrame({'Pdb': keep_next_round, 'Index_made': index })
    temp_df=temp_df.sort_values(by='Index_made')
    
    temp_df = temp_df.reset_index()
    temp_df=temp_df[['Pdb','Index_made']]
    
    
    return(temp_df)

def create_indlist2(X,mutate_position,data):
    with open('individual_list3.txt','w') as myfile:
        for j in find_best_binders(-6.30,data)[1]:
            for i in aa:
                myfile.write(j[:-1]+","+"A"+X+str(mutate_position)+i+";\n")
    myfile.close()
    
def generate_FOLDX_mutate_runfiles(data,position):
    for i in range(0,len(data)):
        temp_file="config_mutate_pos"+str(position)+"_pdb"+str(data.Index_made[i])+"_hsp70.cfg"
        with open(temp_file,'w',newline='\n') as configfile:
            print(data.iloc[i,0][2:])
            text_list_config = ["command=BuildModel\n","pdb="+data.iloc[i,0][2:]+"\n",
                "pdb-dir=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position-1)+"\n",
                "rotabaseLocation=/home/krinad/thesisfiles/hsp70/rotabase.txt\n",
                "mutant-file=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position)+"/individual_list.txt\n",
                "output-dir=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position)+"\n"]
            configfile.writelines(text_list_config)
            configfile.close()
        
    for i in range(0,len(data)):
        temp_file_q="FOLDX_mutate_pos"+str(position)+"_pdb"+str(data.Index_made[i])+"_hsp70.q"
        with open(temp_file_q,'w',newline='\n') as qfile:
            text_list_q = ["#!/bin/bash\n",
                           "#$ -N FOLDX_4po2AC_pos"+str(position)+"_hc\n",
                           "#$ -cwd\n",
                           "#$ -V\n" ,
                           "#$ -q all.q\n" ,
                           "#$ -l h_vmem=2G\n",
                           "\n",
                           "/switchlab/group/tools/2021_FoldX5_LoopX/foldx -f /home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position)+"/config_mutate_pos"+str(position)+"_pdb"+str(data.Index_made[i])+"_hsp70.cfg\n"]
            qfile.writelines(text_list_q)
            qfile.close()

def generate_FOLDX_analyse_runfiles(data,position):
            ## create config file for analyze complex
            for i in range(0,len(data)):
                temp_file="config_analyse_pos"+str(position)+"_pdb"+str(data.Index_made[i])+"_hsp70.cfg"
                with open(temp_file,'w',newline='\n') as configfile:
                    print(data.iloc[i,0][2:])
                    text_list_config = ["command=AnalyseComplex\n","pdb-list/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position)+"/pdb-list.txt\n",
                        "rotabaseLocation=/home/krinad/thesisfiles/hsp70/rotabase.txt\n",
                        "analyseComplexChains=A,C",
                        "output-dir=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position)+"/analyse\n"]
                    configfile.writelines(text_list_config)
                    configfile.close()
            # create foldx job for analyze complex
            for i in range(0,len(data)):
                    temp_file_q="FOLDX_analyse_pos"+str(position)+"_pdb"+str(data.Index_made[i])+"_hsp70.q"
                    with open(temp_file_q,'w',newline='\n') as qfile:
                        text_list_q = ["#!/bin/bash\n",
                                       "#$ -N FOLDX_4po2AC__ana_pos"+str(position)+"_hc\n",
                                       "#$ -cwd\n",
                                       "#$ -V\n" ,
                                       "#$ -q all.q\n" ,
                                       "#$ -l h_vmem=2G\n",
                                       "\n",
                                       "/switchlab/group/tools/2021_FoldX5_LoopX/foldx -f /home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position)+"/config_analyse_pos"+str(position)+"_pdb"+str(data.Index_made[i])+"_hsp70.cfg\n"]
                        qfile.writelines(text_list_q)
                        qfile.close()
                        
def generate_FOLDX_analyse_runfiles_mod(data,position):
    temp_file="config_analyse_pos"+str(position)+"_pdb"+"_hsp70.cfg"
    with open(temp_file,'w',newline='\n') as configfile:
        text_list_config = ["command=AnalyseComplex\n","pdb-list=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position)+"/pdb-list.txt\n",
            "rotabaseLocation=/home/krinad/thesisfiles/hsp70/rotabase.txt\n",
            "analyseComplexChains=A,C\n",
            "output-dir=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position)+"/analyse\n"]
        configfile.writelines(text_list_config)
        configfile.close()
    temp_file_q="FOLDX_analyse_pos"+str(position)+"_pdb"+"_hsp70.q"
    with open(temp_file_q,'w',newline='\n') as qfile:
        text_list_q = ["#!/bin/bash\n",
                       "#$ -N FOLDX_4po2AC__ana_pos"+str(position)+"_hc\n",
                       "#$ -cwd\n",
                       "#$ -V\n" ,
                       "#$ -q all.q\n" ,
                       "#$ -l h_vmem=2G\n",
                       "\n",
                       "/switchlab/group/tools/2021_FoldX5_LoopX/foldx -f /home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position)+"/config_analyse_pos"+str(position)+"_pdb"+"_hsp70.cfg\n"]
        qfile.writelines(text_list_q)
        qfile.close()
    temp_file="pdb-list.txt"
    with open(temp_file,'w',newline='\n') as pdbfile:
        for i in data.Pdb:
            for j in range(1,21):
                pdbfile.write("/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position)+"/"+i[2:len(i)-4]+"_"+str(j)+".pdb;\n")
    pdbfile.close()
    
        
    
         #"/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position)+"/"+   

    
        

#!/bin/bash
#$ -N FOLDX_4po2AC_pos1_hc
#$ -cwd
#$ -V
#$ -q all.q
#$ -l h_vmem=2G

# /switchlab/group/tools/2021_FoldX5_LoopX/foldx -f /home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos2/config_mutate_pos2_hsp70.cfg


#################

### ROUND 1

#create_indlist("C", 2)
data=readSummary("C:\\Users\\user\\Desktop\\hill_climbing\\position2\\Summary_pos2_hc_file.txt",2)
find_next=find_best_binders(data)

## ROUND 2
#create_indlist2("C", 3, data) # this is if you want to put them all in the same job
#create_indlist("C", 3)
#generate_FOLDX_mutate_runfiles(find_next, 3)
#generate_FOLDX_analyse_runfiles(find_next,3)
#generate_FOLDX_analyse_runfiles_mod(find_next,3)

data3=readSummary("C:\\Users\\user\\Desktop\\hill_climbing\\position3\\Summary_analyse_pos3.txt",2)


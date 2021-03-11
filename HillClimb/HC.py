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
import subprocess
## mutate 1st position
def create_indlist(X,mutate_position):  
    with open('individual_list.txt','w',newline='\n') as myfile:
        for i in aa:
            myfile.write("A"+X+str(mutate_position)+i+";\n")
    myfile.close()
   
def readSummary(location_file):
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

def create_indlist2(X,mutate_position,best_binders):
    with open('individual_list.txt','w',newline='\n') as myfile:
            for i in aa:
                myfile.write("A"+X+str(mutate_position)+i+";\n")
    myfile.close()
    
def generate_FOLDX_mutate_runfiles_pos2(originalPDB,position):
    
    temp_file="config_mutate_pos"+str(position)+"_pdb"+"_original"+"_hsp70.cfg"
    with open(temp_file,'w',newline='\n') as configfile:
                text_list_config = ["command=BuildModel\n","pdb="+originalPDB+"\n",
                    "pdb-dir=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/\n",
                    "rotabaseLocation=/home/krinad/thesisfiles/hsp70/rotabase.txt\n",
                    "mutant-file=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position)+"/individual_list.txt\n",
                    "output-dir=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position)+"\n"]
                configfile.writelines(text_list_config)
                configfile.close()
      
    temp_file_q="FOLDX_mutate_pos"+str(position)+"_pdb"+"_original"+"_hsp70.q"
    with open(temp_file_q,'w',newline='\n') as qfile:
        text_list_q = ["#!/bin/bash\n",
                       "#$ -N FOLDX_4po2AC_pos"+str(position)+"_hc\n",
                       "#$ -cwd\n",
                       "#$ -V\n" ,
                       "#$ -q all.q\n" ,
                       "#$ -l h_vmem=2G\n",
                       "\n",
                       "/switchlab/group/tools/2021_FoldX5_LoopX/foldx -f /home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position)+"/"+temp_file+"\n"]
        qfile.writelines(text_list_q)
        qfile.close()

def generate_FOLDX_analyse_runfiles_pos2(originalPDB, position):
    temp_file="config_analyse_pos"+str(position)+"_pdb_"+"original"+"_hsp70.cfg"
    with open(temp_file,'w',newline='\n') as configfile:
        text_list_config = ["command=AnalyseComplex\n","pdb-list=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position)+"/pdb-list.txt\n",
            "rotabaseLocation=/home/krinad/thesisfiles/hsp70/rotabase.txt\n",
            "analyseComplexChains=A,C\n",
            "output-dir=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position)+"/analyse/\n"]
        configfile.writelines(text_list_config)
        configfile.close()

    temp_file_q="FOLDX_analyse_pos"+str(position)+"_pdb_"+"original"+"_hsp70.q"
    with open(temp_file_q,'w',newline='\n') as qfile:
        text_list_q = ["#!/bin/bash\n",
                       "#$ -N FOLDX_4po2AC_pos"+str(position)+"_ana_hc\n",
                       "#$ -cwd\n",
                       "#$ -V\n" ,
                       "#$ -q all.q\n" ,
                       "#$ -l h_vmem=2G\n",
                       "\n",
                       "/switchlab/group/tools/2021_FoldX5_LoopX/foldx -f /home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position)+"/config_analyse_pos"+str(position)+"_pdb_"+"original"+"_hsp70.cfg\n"]
        qfile.writelines(text_list_q)
        qfile.close()
    temp_file="pdb-list.txt"
    with open(temp_file,'w',newline='\n') as pdbfile:
            for j in range(1,21):
                pdbfile.write("/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position)+"/"+"4po2AC_Repair_"+str(j)+".pdb;\n")
            pdbfile.close()
    
        
    
def generate_FOLDX_mutate_runfiles(data,position):
    for i in range(0,len(data)):
            #temp_file="config_mutate_pos"+str(position)+"_pdb_"+str(data.Pdb[i].split("_")[-1].split(".")[0])+"_hsp70.cfg"
            temp_file="config_mutate_pos"+str(position)+"_pdb_"+str(data.Pdb[i][16:len(data.Pdb[i])-4])+"_hsp70.cfg"
            print(temp_file)
            with open(temp_file,'w',newline='\n') as configfile:         
               
                text_list_config = ["command=BuildModel\n","pdb="+data.iloc[i,0][2:]+"\n",
                    "pdb-dir=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position-1)+"\n",
                    "rotabaseLocation=/home/krinad/thesisfiles/hsp70/rotabase.txt\n",
                    "mutant-file=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position)+"/individual_list.txt\n",
                    "output-dir=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position)+"\n"]
                configfile.writelines(text_list_config)
                configfile.close()
                 
    for i in range(0,len(data)):
        temp_file_q="FOLDX_mutate_pos"+str(position)+"_pdb_"+str(data.Pdb[i][16:len(data.Pdb[i])-4])+"_hsp70.q"
        with open(temp_file_q,'w',newline='\n') as qfile:
            text_list_q = ["#!/bin/bash\n",
                           "#$ -N FOLDX_4po2AC_pos"+str(position)+"_hc\n",
                           "#$ -cwd\n",
                           "#$ -V\n" ,
                           "#$ -q all.q\n" ,
                           "#$ -l h_vmem=2G\n",
                           "\n",
                           "/switchlab/group/tools/2021_FoldX5_LoopX/foldx -f /home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position)+"/config_mutate_pos"+str(position)+"_pdb_"+str(data.Pdb[i][16:len(data.Pdb[i])-4])+"_hsp70.cfg"+"\n"]
            qfile.writelines(text_list_q)
            qfile.close()
    with open("commandfile1.txt","w",newline="\n") as commandfile:
        for i in range(0,len(data)):
            temp_file_q="FOLDX_mutate_pos"+str(position)+"_pdb_"+str(data.Pdb[i][16:len(data.Pdb[i])-4])+"_hsp70.q"
            commandfile.write("qsub " +temp_file_q+"\n")
    commandfile.close()    
                

    

def generate_FOLDX_analyse_runfiles(data,position):
            ## create config file for analyze complex
            for i in range(0,len(data)):
                temp_file="config_analyse_pos"+str(position)+"_pdb_"+str(data.Index_made[i])+"_hsp70.cfg"
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
                    temp_file_q="FOLDX_analyse_pos"+str(position)+"_pdb_"+str(data.Index_made[i])+"_hsp70.q"
                    with open(temp_file_q,'w',newline='\n') as qfile:
                        text_list_q = ["#!/bin/bash\n",
                                       "#$ -N FOLDX_4po2AC__ana_pos"+str(position)+"_hc\n",
                                       "#$ -cwd\n",
                                       "#$ -V\n" ,
                                       "#$ -q all.q\n" ,
                                       "#$ -l h_vmem=2G\n",
                                       "\n",
                                       "/switchlab/group/tools/2021_FoldX5_LoopX/foldx -f /home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position)+"/config_analyse_pos"+str(position)+"_pdb_"+str(data.Index_made[i])+"_hsp70.cfg\n"]
                        qfile.writelines(text_list_q)
                        qfile.close()
                        
def generate_FOLDX_analyse_runfiles_mod(data,position):
    temp_file="config_analyse_pos"+str(position)+"_pdb_all_hsp70.cfg"
    with open(temp_file,'w',newline='\n') as configfile:
        text_list_config = ["command=AnalyseComplex\n","pdb-list=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position)+"/pdb-list.txt\n",
            "rotabaseLocation=/home/krinad/thesisfiles/hsp70/rotabase.txt\n",
            "analyseComplexChains=A,C\n",
            "output-dir=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position)+"/analyse\n"]
        configfile.writelines(text_list_config)
        configfile.close()
    temp_file_q="FOLDX_analyse_pos"+str(position)+"_pdb_all_hsp70.q"
    with open(temp_file_q,'w',newline='\n') as qfile:
        text_list_q = ["#!/bin/bash\n",
                       "#$ -N FOLDX_4po2AC__ana_pos"+str(position)+"_hc\n",
                       "#$ -cwd\n",
                       "#$ -V\n" ,
                       "#$ -q all.q\n" ,
                       "#$ -l h_vmem=2G\n",
                       "\n",
                       "/switchlab/group/tools/2021_FoldX5_LoopX/foldx -f /home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position)+"/config_analyse_pos"+str(position)+"_pdb_all_hsp70.cfg\n"]
        qfile.writelines(text_list_q)
        qfile.close()
    temp_file="pdb-list.txt"
    with open(temp_file,'w',newline='\n') as pdbfile:
        for i in data.Pdb:
            for j in range(1,21):
                pdbfile.write("/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/hill_climb/pos"+str(position)+"/"+i[2:len(i)-4]+"_"+str(j)+".pdb;\n")
    pdbfile.close()
    

#################

### ROUND 2-- position 2
# create individual list for position 2
### originalPDB : 4po2AC_Repair_1_0.pdb--> rename: 4po2AC_Repair.pdb

#create_indlist("C", 2)
#generate_FOLDX_mutate_runfiles_pos2("4po2AC_Repair.pdb",2)
## BASH SCRIPT
#mkdir analyse
#generate_FOLDX_analyse_runfiles_pos2("4po2AC_Repair.pdb",2)
### BASH SCRIPT
# cd analyse
# mkdir Summary
# mv Summary_* ./Summary
# cd Summary
# tail -n2 Summary_homologymodel_4po2AC_hsc70_Repair_1_0_1_1_AC.fxout | head -n1  > Summary_file_pos2.txt
# tail -n1 -q *.fxout >> Summary_file_pos2.txt
data=readSummary("C:\\Users\\user\\Desktop\\hill_climbing\\pos2\\Summary_file_pos2.txt")
find_next=find_best_binders(data)

## ROUND 3 -- position 3
#create_indlist2("C", 3, data) # this is if you want to put them all in the same job
#create_indlist("C", 3)
#generate_FOLDX_mutate_runfiles(find_next, 3)
#generate_FOLDX_analyse_runfiles_mod(find_next,3)

data3=readSummary("C:\\Users\\user\\Desktop\\hill_climbing\\pos3\\Summary_file_pos3.txt")
find_next=find_best_binders(data3)

## ROUND 4 -- position 4
#create_indlist("C", 4)
#generate_FOLDX_mutate_runfiles(find_next, 4)
#generate_FOLDX_analyse_runfiles_mod(find_next,4)

data4=readSummary("C:\\Users\\user\\Desktop\\hill_climbing\\pos4\\Summary_file_pos4.txt")
find_next=find_best_binders(data4)

## ROUND 5 -- position 5
#create_indlist("C", 5)
#generate_FOLDX_mutate_runfiles(find_next, 5)
#generate_FOLDX_analyse_runfiles_mod(find_next,5)
data5=readSummary("C:\\Users\\user\\Desktop\\hill_climbing\\pos5\\Summary_file_pos5.txt")
find_next=find_best_binders(data5)

## ROUND 6 -- position 6
create_indlist("C", 6)
generate_FOLDX_mutate_runfiles(find_next, 6)



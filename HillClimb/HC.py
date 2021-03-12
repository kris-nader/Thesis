# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 17:50:49 2021

@author: Kristen Michelle Nader
Hill Climbing trial ~~ i dont see the stochasticity here so??
"""

aa=["A","R","N","D","C","E","Q","G","H",
    "I","L","K","M","F","P","S","T","W","Y","V"]

dic_aa={"ALA":'A',"ARG":'R',"ASN":'N',"ASP":'D',"CYS":'C',"GLU":'E',"GLN":'Q',
        "GLY":'G',"HIS":'H',"ILE":'I',"LEU":'L',"LYS":'K',"MET":'M',"PHE":'F',
        "PRO":'P',"SER":'S',"THR":'T',"TRP":"W","TYR":'Y',"VAL":'V'}

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

def find_top_X(data,X):
    data=data.sort_values(by='Interaction_Energy')
    data = data.reset_index()
    data=data[['Pdb','Interaction_Energy']]
    return(data.iloc[0:int(X),])
    
            
    
    
    

# def create_indlist2(X,mutate_position,best_binders):
#     with open('individual_list.txt','w',newline='\n') as myfile:
#             for i in aa:
#                 myfile.write("A"+X+str(mutate_position)+i+";\n")
#     myfile.close()
    
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

def find_substate_read_file(location_file):
    pdb_data=pd.read_csv(location_file, sep="\s+", header=None, skiprows=1778)
    pdb_data.columns=["ATOM", "Position", "atom_type","aa","chain","substrate_position","one","two","three","four","five"]
    substrate=[]
    count=1
    for i in range(0,len(pdb_data.aa)):
        if(count!=int(pdb_data.substrate_position[i])):
            substrate.append(dic_aa[pdb_data.aa[i]])
            count=count+1
    return("".join(substrate))
    
def find_substate_name(file,length_of_substrate):
    substrate=['A']*length_of_substrate
    temp=file[:-4].split("_")[2:]
    for i in range(0,len(temp)):
        substrate[i]=aa[int(temp[i])-1]
    return("".join(substrate))
        
def return_top_X_peptides(top_X,length_of_substrate):
    peptides=[]
    for i in range(0,len(top_X)):
        peptides.append(find_substate_name(top_X.Pdb[i],6))
    return (peptides)    
    

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
top_10_2=find_top_X(data,10)
peptides_data2=return_top_X_peptides(top_10_2,6)

## ROUND 3 -- position 3
#create_indlist2("C", 3, data) # this is if you want to put them all in the same job
#create_indlist("C", 3)
#generate_FOLDX_mutate_runfiles(find_next, 3)
#generate_FOLDX_analyse_runfiles_mod(find_next,3)

data3=readSummary("C:\\Users\\user\\Desktop\\hill_climbing\\pos3\\Summary_file_pos3.txt")
find_next=find_best_binders(data3)
top_10_3=find_top_X(data3,10)
peptides_data3=return_top_X_peptides(top_10_3,6)

## ROUND 4 -- position 4
#create_indlist("C", 4)
#generate_FOLDX_mutate_runfiles(find_next, 4)
#generate_FOLDX_analyse_runfiles_mod(find_next,4)

data4=readSummary("C:\\Users\\user\\Desktop\\hill_climbing\\pos4\\Summary_file_pos4.txt")
find_next=find_best_binders(data4)
top_10_4=find_top_X(data4,10)
peptides_data4=return_top_X_peptides(top_10_4,6)

## ROUND 5 -- position 5
#create_indlist("C", 5)
#generate_FOLDX_mutate_runfiles(find_next, 5)
#generate_FOLDX_analyse_runfiles_mod(find_next,5)
data5=readSummary("C:\\Users\\user\\Desktop\\hill_climbing\\pos5\\Summary_file_pos5.txt")
find_next=find_best_binders(data5)
top_10_5=find_top_X(data5,10)
peptides_data5=return_top_X_peptides(top_10_5,6)


## ROUND 6 -- position 6
#create_indlist("C", 6)
#generate_FOLDX_mutate_runfiles(find_next, 6)
#generate_FOLDX_analyse_runfiles_mod(find_next,6)

data6=readSummary("C:\\Users\\user\\Desktop\\hill_climbing\\pos6\\Summary_file_pos6.txt")
find_next=find_best_binders(data6)


data6.loc[data6['Interaction_Energy'] == min(data6.Interaction_Energy)]
#print(find_substate_name("./4po2AC_Repair_19_4_1_13.pdb ",6))
top_10_6=find_top_X(data6,10)
peptides_data6=return_top_X_peptides(top_10_6,6)



peptides_6rounds=pd.DataFrame(list(zip(peptides_data2, peptides_data3,
                                       peptides_data4,peptides_data5,
                                       peptides_data6)),
                              columns =['pos2', 'po3','pos4','pos5','pos6'])



## ROUND 7 -- position 7
#create_indlist("C", 7)
#generate_FOLDX_mutate_runfiles(find_next, 7)



#find_sub=find_substate("C:\\Users\\user\\Desktop\\hill_climbing\\pos3\\4po2AC_Repair_19_1.pdb")


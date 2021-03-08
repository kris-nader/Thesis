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
    data=pd.read_csv("C:\\Users\\user\\Desktop\\hill_climbing\\position2\\Summary_pos2_hc_file.txt", sep="\t")
    data.columns=["Pdb", "Group1", "Group2","IntraclashesGroup1","IntraclashesGroup2","Interaction_Energy","StabilityGroup1","StabilityGroup2"]
    individual_list=pd.read_csv("C:\\Users\\user\\Desktop\\hill_climbing\\position2\\individual_list.txt", sep="\t",header=None)
    temp=[]
    for i in range(len(data)) : 
      temp.append(int(data.iloc[i, 0].split("_")[4].split(".")[0]))
    data["split"]=temp
    data=data.sort_values(by='split')
    data=data[['Pdb','Interaction_Energy']]
    data["aa_mod"]=list(individual_list[0])
    return (data)
    
def find_best_binders(data):
    keep_next_round=[]
    which_aa=[]
    file=data.Pdb
    interaction_energy=data.Interaction_Energy
    aa_find=data.aa_mod
    fitness_base=max(data.Interaction_Energy)
    for i in range(len(file)):
        if interaction_energy[i]<fitness_base:
            keep_next_round.append(file[i])
            which_aa.append(aa_find[i])
    
    return(keep_next_round,which_aa)

def create_indlist2(X,mutate_position,data):
    with open('individual_list3.txt','w') as myfile:
        for j in find_best_binders(-6.30,data)[1]:
            for i in aa:
                myfile.write(j[:-1]+","+"A"+X+str(mutate_position)+i+";\n")
    myfile.close()

             
#################

### ROUND 1

#create_indlist("C", 2)
data=readSummary("C:\\Users\\user\\Desktop\\hill_climbing\\position2\\Summary_pos2_hc_file.txt",2)
find_next=find_best_binders(data)

## ROUND 2
create_indlist2("C", 3, data)



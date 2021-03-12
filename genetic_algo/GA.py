# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 10:45:39 2021

@author: Kristen Michelle Nader
"""
import random
import pandas as pd
import math

aa=["A","R","N","D","C","E","Q","G","H",
    "I","L","K","M","F","P","S","T","W","Y","V"]

hydrophobicaa=["A","V","I","L","M","F","Y","W"]

## Start with n(100 in this case) randomized peptide sequences

def create_score(temp_string):
    # computes basic score: number of hydrophobic amino acids in the seq
    # returns sum for the sequence
    temp_dict={"A":0, "V":0,"I":0, "L":0,"M":0, "F":0, "Y":0, "W":0}
    for i in list(temp_string):
        #print(i)
        if i in temp_dict:
            temp_dict[i]=temp_dict[i]+1
    return sum(temp_dict.values())


def create_seq(n,k):
    # n is the number of sequences to create
    # k is the length of the sequence
    # returns a dictionary of the score and the sequence
    randomlist = {}
    for j in range(0,n): 
        temp="".join(random.sample(aa, k))
        if temp not in randomlist:
            randomlist[temp]=create_score(temp)
    return randomlist

def create_indlist(randomlist,X,start_position,end_position):  
     with open('individual_list.txt','w',newline='\n') as myfile:
            for i in randomlist:
                count=2
                temp_line=[]
                for k in list(i):
                    #print (k)
                    if (count==7):
                        temp_line.append("A"+X+str(count)+k+";\n")
                    else:
                        temp_line.append("A"+X+str(count)+k+",")
                    count=count+1
                myfile.write("".join(temp_line))
                        
            myfile.close()
                

### so i cant just take the best interactions and 
###combine them bc i might get stuck in a local optimum
### eliteism: suppose 60%  of best chrom and pass them untouched and then the rest 
### 100__> 20 passed and 100 cross over

def compute_fitness(location_file):
    # based on the ineraction energy
    data=pd.read_csv(location_file, sep="\t")
    temp=[]
    for i in range(len(data)) : 
      temp.append(int(data.iloc[i, 0].split("_")[4].split(".")[0]))
    data["split"]=temp
    data=data.sort_values(by='split')
    data["peptideseq"]=list(randomlist.keys())
    data=data[['peptideseq','Interaction_Energy']]
    data=data.sort_values(by='Interaction_Energy')
    return (data)

def selection(peptides,interaction_energy,n):
    # peptides: the peptides and the pdb names
    # fitness: interaciton energy
    # n is the number of peptides to be selected
    # perform a roulette selection
    nselected=0
    percentages=[f/sum(abs(interaction_energy)) for f in abs(interaction_energy)] ##??
    selectedIndividuals = []
    while nselected!=n:
        for i in range(0, len(interaction_energy)):
            if nselected==n:
                break
            x=random.random()
            if (x<=percentages[i]):
                selectedIndividuals.append(list(peptides)[i])
                nselected=nselected+1
    return selectedIndividuals


## this one will split in a difffernt way
### pep1: 'R-AT-PFV'  lets say CuttingPoint=1  child1: 'R-NP-PFV'  
### pep2: 'S-NP-ELY'           CuttingPoint2=3 child2: 'S-AT-ELY'

## problem: what if the cutt1> cutt 2?? --- LOOK AT THIS
def crossover(peptides):
    children = []
    #data.peptideseq
    for i in range(0,len(peptides),2):
        prob = random.random()
        if prob<=crossoverProbability:
            cuttingPoint=random.randint(0, nvar-1)
            cuttingPoint2=random.randint(0, nvar-1)
            child1= peptides[i][0:cuttingPoint] + peptides[i+1][cuttingPoint:cuttingPoint2] + peptides[i][cuttingPoint2:]
            child2= peptides[i+1][0:cuttingPoint] + peptides[i][cuttingPoint:cuttingPoint2] + peptides[i+1][cuttingPoint2:]
            if(len(child1)>6 or len(child2)>6):
                while( len(child1)>6 or len(child2)>6):
                    cuttingPoint=random.randint(0, nvar-1)
                    cuttingPoint2=random.randint(0, nvar-1)
                    child1= peptides[i][0:cuttingPoint] + peptides[i+1][cuttingPoint:cuttingPoint2] + peptides[i][cuttingPoint2:]
                    child2= peptides[i+1][0:cuttingPoint] + peptides[i][cuttingPoint:cuttingPoint2] + peptides[i+1][cuttingPoint2:]
                children.append(child1)
                children.append(child2)
            else:
                children.append(child1)
                children.append(child2)
                    
        else:
            print("enetered")
            print(peptides[i])
            children.append(peptides[i])
            children.append(peptides[i+1])
    return(children)



#METHOD1:
## this one will split in a difffernt way
### pep1: 'R-AT-PFV'  lets say CuttingPoint=1  child1: 'R-NP-PFV'  
### pep2: 'S-NP-ELY'           CuttingPoint2=3 child2: 'S-AT-ELY'

#METHOD2:
## this one will split in half so
### pep1: 'RAT-PFV'   child1: 'RAT-ELY'
### pep2: 'SNP-ELY'   child2: 'SNP-PFV'

def crossover_half(peptides):
    children = []
    for i in range(0,len(peptides),2):
        prob = random.random()
        if prob<=crossoverProbability:
            cuttingPoint=int(len(peptides[i])/2)
            child1=peptides[i][0:cuttingPoint+1]+peptides[i+1][cuttingPoint+1:]
            child2=peptides[i+1][0:cuttingPoint+1]+peptides[i][cuttingPoint+1:]
            children.append(child1)
            children.append(child2)           
        else:
            children.append(peptides[i])
            children.append(peptides[i+1])
    return (children)


def mutate(peptides):
    mutatedPeptides = []
    local_mutate=random.random()
    for peptide in peptides:
         local_mutate=random.random()
         print(local_mutate)
         if local_mutate<=mutationRate:
             position_mutate=random.randint(0,5)
             print(position_mutate)
             aa_mutate=random.sample(aa,1)[0]
             temp_peptide = list(peptide)
             temp_peptide[position_mutate]= aa_mutate
             #computeFitness("".join(list(temp_peptide)
             if (get_interaction_energy("").join(temp_peptide))< maxFitness):
                 peptide = "".join(temp_peptide)
                 mutatedPeptides.append(peptide)
    return (peptides)
    
def get_interaction_energy(peptide):
    generate_single_individual_list(peptide)
    generate_FOLDX_mutate_runfiles("4po2AC_Repair.pdb",0)
    generate_FOLDX_analyse_runfiles("4po2AC_Repair.pdb",0)
    temp_ie=pd.read_csv("C:\\Users\\user\\Desktop\\GA\\Summary_4po2AC_Repair_1_AC.fxout", skiprows=8, sep="\t")
    temp_ie.columns=["Pdb", "Group1", "Group2","IntraclashesGroup1","IntraclashesGroup2","Interaction_Energy","StabilityGroup1","StabilityGroup2"]
    return(temp_ie.Interaction_Energy[0])


def generate_single_individual_list(peptide):
    with open('individual_list.txt','w',newline='\n') as myfile:
        temp_line=[]
        count=2
        for i in peptide:
            if (count==7):
                temp_line.append("A"+"C"+str(count)+i+";\n")
            else:
                    temp_line.append("A"+"C"+str(count)+i+",")
                    count=count+1
        myfile.write("".join(temp_line))
        myfile.close()
    
def generate_FOLDX_analyse_runfiles(originalPDB, nround):
    temp_file="config_analyse_pos"+str(nround)+"_pdb_"+"original"+"_hsp70.cfg"
    with open(temp_file,'w',newline='\n') as configfile:
        text_list_config = ["command=AnalyseComplex\n","pdb-list=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_round"+str(nround)+"/pdb-list.txt\n",
            "rotabaseLocation=/home/krinad/thesisfiles/hsp70/rotabase.txt\n",
            "analyseComplexChains=A,C\n",
            "output-dir=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_round"+str(nround)+"/analyse/\n"]
        configfile.writelines(text_list_config)
        configfile.close()

    temp_file_q="FOLDX_analyse_pos"+str(nround)+"_pdb_"+"original"+"_hsp70.q"
    with open(temp_file_q,'w',newline='\n') as qfile:
        text_list_q = ["#!/bin/bash\n",
                       "#$ -N FOLDX_4po2AC_pos"+str(nround)+"_ana_hc\n",
                       "#$ -cwd\n",
                       "#$ -V\n" ,
                       "#$ -q all.q\n" ,
                       "#$ -l h_vmem=256M\n",
                       "\n",
                       "/switchlab/group/tools/2021_FoldX5_LoopX/foldx -f /home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_round"+str(nround)+"/config_analyse_pos"+str(nround)+"_pdb_"+"original"+"_hsp70.cfg\n"]
        qfile.writelines(text_list_q)
        qfile.close()
    temp_file="pdb-list.txt"
    with open(temp_file,'w',newline='\n') as pdbfile:
        pdbfile.write("/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_round"+str(nround)+"/"+"4po2AC_Repair_1.pdb;\n")
        pdbfile.close()
    
 
def generate_FOLDX_mutate_runfiles(originalPDB,nround):
    temp_file="config_mutate_pos"+str(nround)+"_pdb"+"_original"+"_hsp70.cfg"
    with open(temp_file,'w',newline='\n') as configfile:
                text_list_config = ["command=BuildModel\n","pdb="+originalPDB+"\n",
                    "pdb-dir=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/\n",
                    "rotabaseLocation=/home/krinad/thesisfiles/hsp70/rotabase.txt\n",
                    "mutant-file=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_round"+str(nround)+"/individual_list.txt\n",
                    "output-dir=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_round"+str(nround)+"\n"]
                configfile.writelines(text_list_config)
                configfile.close()
      
    temp_file_q="FOLDX_mutate_pos"+str(nround)+"_pdb"+"_original"+"_hsp70.q"
    with open(temp_file_q,'w',newline='\n') as qfile:
        text_list_q = ["#!/bin/bash\n",
                            "#$ -N FOLDX_4po2AC_pos"+str(nround)+"_hc\n",
                            "#$ -cwd\n",
                            "#$ -V\n" ,
                            "#$ -q all.q\n" ,
                            "#$ -l h_vmem=512M\n",
                            "\n",
                            "/switchlab/group/tools/2021_FoldX5_LoopX/foldx -f /home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_round"+str(nround)+"/"+temp_file+"\n"]
        qfile.writelines(text_list_q)
        qfile.close()
   
###################
#generate 100 random individuals
#randomlist=create_seq(100,6)

crossoverProbability=0.8
nvar=6
maxFitness=-12 ### change to alanine
elitismRate=0.6
populationSize=100
mutationRate=0.015


## 1st INITIALIZATION

randomlist={'RQMSFH': 2,
 'MRYSWF': 4, 'HSAQYG': 2, 'IDTRVM': 3, 'KFPAIQ': 3, 'REKPSA': 1, 'WMQHIY': 4, 'YAQGVE': 3,
 'QHDKAN': 1, 'QYCPGT': 1, 'VCMPWR': 3, 'PNHTWG': 1, 'VFENDP': 2, 'WQNLTC': 2, 'KWNGVE': 2, 
  'NTLGFY': 3, 'REPGAQ': 1, 'DCKIHR': 1, 'CKYHDI': 2, 'FLNGCM': 3, 'PIKLGA': 3, 'QNMVKC': 2, 
  'CRLMDT': 2, 'MKILTE': 3, 'YTCVAP': 3, 'KSNMQH': 1, 'HIERYD': 2, 'DQWCSN': 1, 'LKEIMP': 3, 
  'GDPSQN': 0, 'RCFMGY': 3, 'RIEYKP': 2, 'ARIMES': 3, 'RATPFV': 3, 'NKYSDA': 2, 'AKHDVF': 3,
  'RAPTKY': 2, 'DFEQHC': 1, 'VMCKLG': 3, 'DKCEGM': 1, 'MKVACE': 3, 'QEDTHP': 0, 'RFCHTV': 2, 
  'EVTLCS': 2, 'DNRPIG': 1, 'LYFTRG': 3, 'NFLPAM': 4,
  'ANMKEV': 3,   'HAFSEG': 2, 'RSCEMA': 2, 'EVHKWF': 3, 'ICMQKD': 2, 'EIYTLR': 3, 'SCDWQY': 2, 'NLWPIC': 3,
  'NWPLCT': 2, 'AIMRPT': 3, 'IPTNDV': 2, 'VKQDER': 1, 'WLSIEK': 3, 'FRVGEY': 3, 'IPGKAW': 3, 'WLTDCQ': 2,
  'GWYHQL': 3, 'FNGVDC': 2, 'MKSQCY': 2, 'FNARVT': 3, 'LSDYGC': 2, 'LYMDQE': 3, 'DYGHCP': 1, 'NCTGVD': 1, 'VWRICL': 4,
  'KITDNG': 1, 'MIPVHF': 4, 'VEAMGF': 4, 'IMTRKQ': 2, 'SNPELY': 2, 'FWYTKD': 3, 'RGHFDK': 1, 'CPMSDF': 2, 'VWYQHA': 4, 'SFHMPL': 3,
  'RILMHS': 3, 'VMSWDF': 4, 'FGKCDI': 2, 'HWAGDI': 3, 'WIAYMV': 6,
  'YMTGHD': 2, 'MHICQK': 2, 'NVSIRC': 2, 'MKDLRP': 2, 'TFYLEN': 3, 'IEVQAT': 3,
  'EDSFAG': 2, 'TACNRK': 1, 'RVFPLN': 3, 'NMCWHQ': 2, 'ECDMWP': 2, 'AGVYHP': 3, 'ILMQEN': 3}
#create_indlist(randomlist,"C",2,7) # create the individual list    


# ## EVALUATION
create_indlist(randomlist.keys(),"C",2,7)
generate_FOLDX_mutate_runfiles("4po2AC_Repair_1.pdb",1)
generate_FOLDX_analyse_runfiles("4po2AC_Repair_1.pdb",1)
data=compute_fitness("C:\\Users\\user\\Desktop\\hsp70\\alanine_mut_buildmodel\\GA_round1\\Summary_GA_file.txt")   
# ## elitism?
parents=data.peptideseq
fitness=data.Interaction_Energy
eliteIndividuals = parents[0: math.ceil(populationSize*elitismRate)]
toSelect = populationSize - math.ceil(populationSize*elitismRate)

# # if we want to use the entire data
# ## SELECCTION
select=selection(data.peptideseq, data.Interaction_Energy, toSelect)
# # use the data that are not considered "elite"
# #select=selection(parents[math.ceil(populationSize*elitismRate):],fitness[math.ceil(populationSize*elitismRate):], toSelect)

# ## CROSSOVER
children_new=crossover_half(select)
print(children_new)

# ## MUTATE
mutated_children=mutate(children_new)
print(mutated_children)

# ## CHECK

# ### CONVERGE :y. A gene is said to have converged when 95% of the population share the same value

# count=0
# for i in range(len(data.Interaction_Energy)):
#     if data.Interaction_Energy[i]<maxFitness:
#         print(data.Interaction_Energy[i])
#         count=count+1


    

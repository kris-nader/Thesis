#/usr/local/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 10:45:39 2021

@author: Kristen Michelle Nader
"""
import random
import pandas as pd
import math
import numpy as np
import subprocess
import time



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
    # returns a list of the sequence
    randomlist = []
    for j in range(0,n): 
        temp="".join(random.sample(aa, k))
        randomlist[j]=temp
    return randomlist

## creates individual list for all elements of the new peptides
def create_indlist(randomlist,X,start_position,end_position,nround):  
    name_list="/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"/individual_list.txt"
    with open(name_list,'w',newline='\n') as myfile:
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
               

# reads the analyse complex file and returns a sorted pandas dataframe where the peptideseq and interaction energy are columns
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

# performs a roullette selection and returns n selected peptides
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
            children.append(peptides[i])
            children.append(peptides[i+1])
    return(children)

def run_files(nround):
    temop_q="qsub FOLDX_mutate_round"+str(nround)+"_pdb"+"_original"+"_hsp70.q"
    subprocess.Popen(temop_q, shell=True)
    check_qstat = subprocess.Popen('qstat',stdout=subprocess.PIPE)
    output_qstat = check_qstat.stdout.read()
    while 'FOLDX_' in output_qstat:
        print ("Waiting for all Mutate jobs to finish")
        time.sleep(10)
        check_qstat = subprocess.Popen('qstat',stdout=subprocess.PIPE)
        output_qstat = check_qstat.stdout.read()
    temop_q="FOLDX_analyse_round"+str(nround)+"_pdb_"+"original"+"_hsp70.q"
    subprocess.call(temop_q, shell=True)
    check_qstat = subprocess.Popen('qstat',stdout=subprocess.PIPE)
    output_qstat = check_qstat.stdout.read()
    while 'FOLDX_' in output_qstat:
        print ("Waiting for all Analyse jobs to finish")
        time.sleep(10)
        check_qstat = subprocess.Popen('qstat',stdout=subprocess.PIPE)
        output_qstat = check_qstat.stdout.read()
    subprocess.call("/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/helper_analyze.sh", shell=True)
    
        


# mutates a peptide sequence 
def mutate(peptides):
    mutatedPeptides = []
    local_mutate=random.random()
    for peptide in peptides:
         local_mutate=random.random()
         if local_mutate<=mutationRate:
             position_mutate=random.randint(0,5)
             aa_mutate=random.sample(aa,1)[0]
             temp_peptide = list(peptide)
             temp_peptide[position_mutate]= aa_mutate
             #computeFitness("".join(list(temp_peptide)
             #get_interaction_energy("").join(temp_peptide)
             #random.randint(-15,0)
             if (get_interaction_energy_mutate(("").join(temp_peptide))< maxFitness):
                 peptide = "".join(temp_peptide)
                 mutatedPeptides.append(peptide)
    return (peptides)
   
##############################################################################################################
##############################################################################################################

# private function for the mutate function
def get_interaction_energy_mutate(peptide):
    generate_single_individual_list(peptide)
    generate_single_FOLDX_mutate_runfiles("4po2AC_Repair.pdb",0)
     ### CALL BASH SCRIPT
     #  mkdir analyse
    generate_single_FOLDX_analyse_runfiles("4po2AC_Repair.pdb",0)
    # CALL BASH SCRIPT
    # # cd analyse
    # # mkdir Summary
    # # mv Summary_* ./Summary
    # # cd Summary
    # # tail -n2 Summary_4po2AC_.fxout | head -n1  > Summary_file.txt
    # # tail -n1 -q *.fxout >> Summary_file.txt
    temp_ie=pd.read_csv("/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round0/Summary_4po2AC_Repair_1_AC.fxout",skip=8, sep="\t")
    temp_ie.columns=["Pdb", "Group1", "Group2","IntraclashesGroup1","IntraclashesGroup2","Interaction_Energy","StabilityGroup1","StabilityGroup2"]
    return(temp_ie.Interaction_Energy[0])

# private function for the get_interaction_energy_mutate() function
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

def generate_single_FOLDX_mutate_runfiles(originalPDB,nround):
    temp_file="/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)/+"config_mutate_round"+str(nround)+"_pdb"+"_original"+"_hsp70.cfg"
    with open(temp_file,'w',newline='\n') as configfile:
                text_list_config = ["command=BuildModel\n","pdb="+originalPDB+"\n",
                    "pdb-dir=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/\n",
                    "rotabaseLocation=/home/krinad/thesisfiles/hsp70/rotabase.txt\n",
                    "mutant-file=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"/individual_list.txt\n",
                    "output-dir=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"\n"]
                configfile.writelines(text_list_config)
                configfile.close()
      
    temp_file_q="/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"/FOLDX_mutate_round"+str(nround)+"_pdb"+"_original"+"_hsp70.q"
    with open(temp_file_q,'w',newline='\n') as qfile:
        text_list_q = ["#!/bin/bash\n",
                            "#$ -N FOLDX_4po2AC_round"+str(nround)+"_hc\n",
                            "#$ -cwd\n",
                            "#$ -V\n" ,
                            "#$ -q all.q\n" ,
                            "#$ -l h_vmem=512M\n",
                            "\n",
                            "/switchlab/group/tools/2021_FoldX5_LoopX/foldx -f /home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"/"+temp_file+"\n"]
        qfile.writelines(text_list_q)
        qfile.close()


def generate_single_FOLDX_analyse_runfiles(originalPDB, nround):
    temp_file="/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"/config_analyse_round"+str(nround)+"_pdb_"+"original"+"_hsp70.cfg"
    with open(temp_file,'w',newline='\n') as configfile:
        text_list_config = ["command=AnalyseComplex\n","pdb-list=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"/pdb-list.txt\n",
            "rotabaseLocation=/home/krinad/thesisfiles/hsp70/rotabase.txt\n",
            "analyseComplexChains=A,C\n",
            "output-dir=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"/analyse/\n"]
        configfile.writelines(text_list_config)
        configfile.close()

    temp_file_q="/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"/FOLDX_analyse_round"+str(nround)+"_pdb_"+"original"+"_hsp70.q"
    with open(temp_file_q,'w',newline='\n') as qfile:
        text_list_q = ["#!/bin/bash\n",
                       "#$ -N FOLDX_4po2AC_round"+str(nround)+"_ana_hc\n",
                       "#$ -cwd\n",
                       "#$ -V\n" ,
                       "#$ -q all.q\n" ,
                       "#$ -l h_vmem=512M\n",
                       "\n",
                       "/switchlab/group/tools/2021_FoldX5_LoopX/foldx -f /home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"/config_analyse_round"+str(nround)+"_pdb_"+"original"+"_hsp70.cfg\n"]
        qfile.writelines(text_list_q)
        qfile.close()
        
    temp_file="/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"/pdb-list.txt"
    with open(temp_file,'w',newline='\n') as pdbfile:
        pdbfile.write("/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"/"+"4po2AC_Repair_1pdb;\n")
    pdbfile.close()


##############################################################################################################
##############################################################################################################


## get data for the nexxt round of the GA returns a pandas dataframe with the interaction energy and peptideseq
def get_interaction_energy_full(peptides,count):
    generate_full_individual_list(peptides)
    generate_FOLDX_mutate_runfiles("4po2AC_Repair.pdb",count)
    generate_FOLDX_analyse_runfiles("4po2AC_Repair.pdb",count)
    run_files(count)
    location_file="/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(count)+"/analyse/Summary/Summary_file.txt"
    temp_dataset=compute_fitness(location_file)
    return(temp_dataset)

# generate the individual list ** could reuse the create_indlist(randomlist, X, start_position, end_position)
# but you need to modify it so that it doesnt take a dictionary but a list
def generate_full_individual_list(peptides):
    with open('/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"/individual_list.txt','w',newline='\n') as myfile:
        for peptide in peptides:
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
    
# generate the mutate files 
def generate_FOLDX_mutate_runfiles(originalPDB,nround):
    temp_file="/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"/config_mutate_round"+str(nround)+"_pdb"+"_original"+"_hsp70.cfg"
    with open(temp_file,'w',newline='\n') as configfile:
                text_list_config = ["command=BuildModel\n","pdb="+originalPDB+"\n",
                    "pdb-dir=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/\n",
                    "rotabaseLocation=/home/krinad/thesisfiles/hsp70/rotabase.txt\n",
                    "mutant-file=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"/individual_list.txt\n",
                    "output-dir=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"\n"]
                configfile.writelines(text_list_config)
                configfile.close()
                
    ## create the summary_analuse file for the analyze step
    temp_analyse="/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/Summary_file.txt"
    with open(temp_analyse,'w',newline='\n') as helper_file:
        helper_file.write("Pdb	Group1	Group2	IntraclashesGroup1	IntraclashesGroup2	Interaction Energy	StabilityGroup1	StabilityGroup2\n")
    helper_file.close()
        
    temp_file_q="/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"/FOLDX_mutate_round"+str(nround)+"_pdb"+"_original"+"_hsp70.q"
    with open(temp_file_q,'w',newline='\n') as qfile:
        text_list_q = ["#!/bin/bash\n",
                            "#$ -N FOLDX_4po2AC_round"+str(nround)+"_hc\n",
                            "#$ -cwd\n",
                            "#$ -V\n" ,
                            "#$ -q all.q\n" ,
                            "#$ -l h_vmem=512M\n",
                            "\n",
                            "/switchlab/group/tools/2021_FoldX5_LoopX/foldx -f /home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"/"+temp_file+"\n",
                            "mkdir analyse\n",
                            "mkdir ./analyse/Summary\n",
                            "mv /home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"/Summary_file.txt .analyse/Summary\n"]
        qfile.writelines(text_list_q)
        qfile.close()    
    
    
# generates the analyse runfiles     
def generate_FOLDX_analyse_runfiles(originalPDB, nround):
    temp_file="/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"config_analyse_round"+str(nround)+"_pdb_"+"original"+"_hsp70.cfg"
    with open(temp_file,'w',newline='\n') as configfile:
        text_list_config = ["command=AnalyseComplex\n","pdb-list=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"/pdb-list.txt\n",
            "rotabaseLocation=/home/krinad/thesisfiles/hsp70/rotabase.txt\n",
            "analyseComplexChains=A,C\n",
            "output-dir=/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"/analyse/\n"]
        configfile.writelines(text_list_config)
        configfile.close()

    temp_file_q="/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"FOLDX_analyse_round"+str(nround)+"_pdb_"+"original"+"_hsp70.q"
    with open(temp_file_q,'w',newline='\n') as qfile:
        text_list_q = ["#!/bin/bash\n",
                       "#$ -N FOLDX_4po2AC_round"+str(nround)+"_ana_hc\n",
                       "#$ -cwd\n",
                       "#$ -V\n" ,
                       "#$ -q all.q\n" ,
                       "#$ -l h_vmem=512M\n",
                       "\n",
                       "/switchlab/group/tools/2021_FoldX5_LoopX/foldx -f /home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"/config_analyse_round"+str(nround)+"_pdb_"+"original"+"_hsp70.cfg\n"]
        qfile.writelines(text_list_q)
        qfile.close()
        
    temp_file="/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"/pdb-list.txt"
    with open(temp_file,'w',newline='\n') as pdbfile:
        for i in range(1,101):
            pdbfile.write("/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"/4po2AC_Repair_"+str(i)+".pdb;\n")
    pdbfile.close()
    
# def run_files(nround):
#     temop_q="qsub FOLDX_mutate_round"+str(nround)+"_pdb"+"_original"+"_hsp70.q"
#     subprocess.Popen(temop_q, shell=True)
#     check_qstat = subprocess.Popen('qstat',stdout=subprocess.PIPE)
#     output_qstat = check_qstat.stdout.read()
#     while 'FOLDX_' in output_qstat:
#         print ("Waiting for all Mutate jobs to finish")
#         time.sleep(10)
#         check_qstat = subprocess.Popen('qstat',stdout=subprocess.PIPE)
#         output_qstat = check_qstat.stdout.read()
#     temop_q="FOLDX_analyse_round"+str(nround)+"_pdb_"+"original"+"_hsp70.q"
#     subprocess.call(temop_q, shell=True)
#     check_qstat = subprocess.Popen('qstat',stdout=subprocess.PIPE)
#     output_qstat = check_qstat.stdout.read()
#     while 'FOLDX_' in output_qstat:
#         print ("Waiting for all Analyse jobs to finish")
#         time.sleep(10)
#         check_qstat = subprocess.Popen('qstat',stdout=subprocess.PIPE)
#         output_qstat = check_qstat.stdout.read()
#     subprocess.call("/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/helper_analyze.sh", shell=True)
    
        

    
##############################################################################################################
##############################################################################################################    
 

# creates dummy interaction energy for testing GA
def create_dummy_interactionEnergy(peptides, prev_mean):
    sampl = np.random.uniform(low=-15.5, high=prev_mean, size=(100,))
    temp_list_for_data = list(zip(peptides, sampl))  
    df = pd.DataFrame(temp_list_for_data, 
                  columns = ['peptideseq', 'Interaction_Energy'])  
    return(df)
    
# returns the number of peptides with Interaction energy< max fitness
# can potentially be used for the convergence
def count_over95(data):
    count=0
    for i in range(len(data.Interaction_Energy)):
        if data.Interaction_Energy[i]<maxFitness:
            count=count+1
    return (count)
    



##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
 

#generate 100 random individuals most prob you should start w like 1000
#randomlist=create_seq(100,6)

crossoverProbability=0.8
nvar=6
maxFitness=-12 ### change to alanine
elitismRate=0.6
populationSize=100
mutationRate=0.015


#### initialize the 1st random list. We will remove this later but I need this for testing
## 1st INITIALIZATION

randomlist=['RQMSFH',
 'MRYSWF', 'HSAQYG', 'IDTRVM', 'KFPAIQ', 'REKPSA', 'WMQHIY', 'YAQGVE',
 'QHDKAN', 'QYCPGT', 'VCMPWR', 'PNHTWG', 'VFENDP', 'WQNLTC', 'KWNGVE', 
  'NTLGFY', 'REPGAQ', 'DCKIHR', 'CKYHDI', 'FLNGCM', 'PIKLGA', 'QNMVKC', 
  'CRLMDT', 'MKILTE', 'YTCVAP', 'KSNMQH', 'HIERYD', 'DQWCSN', 'LKEIMP', 
  'GDPSQN', 'RCFMGY', 'RIEYKP', 'ARIMES', 'RATPFV', 'NKYSDA', 'AKHDVF',
  'RAPTKY', 'DFEQHC', 'VMCKLG', 'DKCEGM', 'MKVACE', 'QEDTHP', 'RFCHTV', 
  'EVTLCS', 'DNRPIG', 'LYFTRG', 'NFLPAM',
  'ANMKEV',   'HAFSEG', 'RSCEMA', 'EVHKWF', 'ICMQKD', 'EIYTLR', 'SCDWQY', 'NLWPIC',
  'NWPLCT', 'AIMRPT', 'IPTNDV', 'VKQDER', 'WLSIEK', 'FRVGEY', 'IPGKAW', 'WLTDCQ',
  'GWYHQL', 'FNGVDC', 'MKSQCY', 'FNARVT', 'LSDYGC', 'LYMDQE', 'DYGHCP', 'NCTGVD', 'VWRICL',
  'KITDNG', 'MIPVHF', 'VEAMGF', 'IMTRKQ', 'SNPELY', 'FWYTKD', 'RGHFDK', 'CPMSDF', 'VWYQHA', 'SFHMPL',
  'RILMHS', 'VMSWDF', 'FGKCDI', 'HWAGDI', 'WIAYMV',
  'YMTGHD', 'MHICQK', 'NVSIRC', 'MKDLRP', 'TFYLEN', 'IEVQAT',
  'EDSFAG', 'TACNRK', 'RVFPLN', 'NMCWHQ', 'ECDMWP', 'AGVYHP', 'ILMQEN']
#create_indlist(randomlist,"C",2,7) # create the individual list    


# ## EVALUATION

### need to create directory
create_indlist(randomlist,"C",2,7,1)
generate_FOLDX_mutate_runfiles("4po2AC_Repair.pdb",1)
generate_FOLDX_analyse_runfiles("4po2AC_Repair.pdb",1)
run_files(1)

data1=compute_fitness("/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/GA_round1/Summary_GA_file.txt")   

# ## elitism?
parents=data1.peptideseq
fitness=data1.Interaction_Energy
eliteIndividuals = parents[0: math.ceil(populationSize*elitismRate)]
toSelect = populationSize - math.ceil(populationSize*elitismRate)

# ## SELECCTION
select=selection(data1.peptideseq, data1.Interaction_Energy, toSelect)

# ## CROSSOVER
children_new=crossover(select)

# ## MUTATE
mutated_children=mutate(children_new)
new_gen= list(eliteIndividuals)+ list(mutated_children)
data2=get_interaction_energy_full(new_gen,1)

# ## CHECK
count=2
##data1.Interaction_Energy.mean()>-15
##count_over95(data1)<95
# while (count_over95(data1)<95):
#     parents=data.peptideseq
#     fitness=data.Interaction_Energy
#     eliteIndividuals = parents[0: math.ceil(populationSize*elitismRate)]
#     toSelect = populationSize - math.ceil(populationSize*elitismRate)
#     select=selection(data.peptideseq, data.Interaction_Energy, toSelect)
#     children_new=crossover_half(select)
#     mutated_children=mutate(children_new)
#     new_gen= list(eliteIndividuals)+ list(mutated_children)
#     data1=create_dummy_interactionEnergy(new_gen,data1.Interaction_Energy.mean())
#     count=count+1
#     print(str(count)+" :"+str(count_over95(data1)))

while (count<4):
    ls_bash_command=["mkdir GA_round"+str(count),"cd GA_round"+str(count)]
    subprocess.Popen(ls_bash_command)
    parents=data2.peptideseq
    fitness=data2.Interaction_Energy
    eliteIndividuals = parents[0: math.ceil(populationSize*elitismRate)]
    toSelect = populationSize - math.ceil(populationSize*elitismRate)
    select=selection(data2.peptideseq, data2.Interaction_Energy, toSelect)
    children_new=crossover(select)
    mutated_children=mutate(children_new)
    new_gen= list(eliteIndividuals)+ list(mutated_children)
    data2=get_interaction_energy_full(new_gen,count)
    subprocess.Popen("cd ..")
    count=count+1

    

with open("/home/krinad/thesisfiles/hsp70/buildmodel_AlanineMutagenesis/GA_trial/peptides_bestGA_testing.txt",'w',newline='\n') as best_peptide_file:
    for i in range(0,len(data2.peptideseq)):
        if (data2.Interaction_Energy[i]):
            best_peptide_file.write(" "+data2.peptideseq[i]+" "+str(data2.Interaction_Energy[i])+"\n")
best_peptide_file.close()




# count=0
# for i in range(len(data.Interaction_Energy)):
#     if data.Interaction_Energy[i]<maxFitness:
#         print(data.Interaction_Energy[i])
#         count=count+1

# with open("peptides_best1.txt",'w',newline='\n') as best_peptide_file:
#     for i in range(0,len(data.peptideseq)):
#         if (data.Interaction_Energy[i]<maxFitness):
#             best_peptide_file.write(" "+data.peptideseq[i]+"\n")
# best_peptide_file.close()
#print(find_substate_name_432567("./4po2AC_Repair_19_4_2_5_7_20.pdb","432567",7))

### we have options on where to converge: 
    # A gene is said to have converged when 95% of the population share the same value
    # We converge if the average of the interaction energy is ~~14
    

# while(count<95):
    

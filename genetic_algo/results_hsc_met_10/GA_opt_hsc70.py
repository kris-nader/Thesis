#/software/shared/apps/general/python/3.5.1/bin/python3.5

# """
# Created on Sun Mar  7 10:45:39 2021

# @author: Kristen Michelle Nader
# """


import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import random
import pandas as pd
import math
import os.path
import numpy as np
import subprocess
import time
import sys
import glob
import shutil 


t0 = time.time()


random.seed(10)

aa=["A","R","N","D","C","E","Q","G","H",
    "I","L","K","M","F","P","S","T","W","Y","V"]

hydrophobicaa=["A","V","I","L","M","F","Y","W"]

# Start with n(100 in this case) randomized peptide sequences

def create_seq(n,k):
    # n is the number of sequences to create
    # k is the length of the sequence
    # returns a list of the sequence
    randomlist = []
    for j in range(0,n): 
        temp="".join(random.sample(aa, k))
        randomlist.append(temp)
    return randomlist

## creates individual list for all elements of the new peptides
def create_indlist(randomlist,X,start_position,end_position,nround):
    temp_indlist="individual"+str(nround)+".txt" 
    with open(temp_indlist,"w") as tempfile1:
        chunks = [randomlist[x:x+1] for x in range(0, len(randomlist), 1)]
        for i in range(0,100):
            tempfile1.write(randomlist[i]+"\n")
            path_to_file="./GA_round"+str(nround)+"/R"+str(i+1)
            completeName = os.path.join(path_to_file,"individual_list.txt")
            #print(completeName)
            #name_list="/home/krinad/thesisfiles/dnak/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"/individual_list.txt"
            with open(completeName,"w") as myfile:
                for i in chunks[i]:
                    count=2
                    temp_line=[]
                    for k in list(i):
                        #print("A"+X+str(count)+k+"\n")
                        if (count==7):
                            temp_line.append("A"+X+str(count)+k+";\n")
                            
                        else:
                            temp_line.append("A"+X+str(count)+k+",")
                            
                        count=count+1
                    myfile.write("".join(temp_line))     
                                          
            myfile.close()
    tempfile1.close()
           


# reads the analyse complex file and returns a sorted pandas dataframe where the peptideseq and interaction energy are columns
def compute_fitness(location_file,peptides,nround):
    # based on the ineraction energy
    data=pd.read_csv(location_file, sep="\s+")
    print("gen: \n")
    temp=[]
    for i in range(len(data)) : 
      temp.append(int(data.iloc[i, 0].split("_")[4]))
    data["split"]=temp
    print(data)
    data=data.sort_values(by='split')
    data = data.reset_index(drop=True)
    data=data[['Pdb','Interaction_Energy','split']]
    data['AvgIE'] = data.groupby('split')['Interaction_Energy'].transform('mean')
    temp_avg=[data.AvgIE[i] for i in range(0, len(data)-2, 3)]
    df = pd.DataFrame(list(zip(list(peptides), temp_avg)),
                   columns =['peptideseq', 'Interaction_Energy'])
    df=df.sort_values(by='Interaction_Energy')
    df = df.reset_index(drop=True)


    return(df)


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
            #random.uniform(0,0.1)
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
            child1= peptides[i][0:min(cuttingPoint,cuttingPoint2)] + peptides[i+1][min(cuttingPoint,cuttingPoint2):max(cuttingPoint,cuttingPoint2)] + peptides[i][max(cuttingPoint,cuttingPoint2):]
            child2= peptides[i+1][0:min(cuttingPoint,cuttingPoint2)] + peptides[i][min(cuttingPoint,cuttingPoint2):max(cuttingPoint,cuttingPoint2)] + peptides[i+1][max(cuttingPoint,cuttingPoint2):]
            children.append(child1)
            children.append(child2)                   
        else:
            children.append(peptides[i])
            children.append(peptides[i+1])
    return(children)


# mutates a peptide sequence 
def mutate(peptides):
    local_mutate=random.random()
    new_peptides=[]
    number_mutated=0
    for peptide in peptides:
        local_mutate=random.random()
        if local_mutate<=mutationRate:
            new_peptide=mutate_v1(peptide)
            temp_interaction_e=get_interaction_energy_mutate(new_peptide)
            if (temp_interaction_e< maxFitness):
                new_peptides.append(new_peptide)
                print(new_peptide+"_"+str(temp_interaction_e))
                number_mutated=number_mutated+1
            remove_file()
        else:
            new_peptides.append(peptide)
    print("Number of mutated children: "+str(number_mutated))
    return (peptides)


def mutate_v1(pep):
    pep_list = list(pep)
    mutation_site = random.randint(0, len(pep_list) - 1)
    pep_list[mutation_site] = random.choice(aa)
    return (''.join(pep_list))

   
def remove_file():
    print("removing")
    shutil.rmtree("/switchlab/group/krinad/hsc70/GA_round0/") 
    subprocess.run("mkdir /switchlab/group/krinad/hsc70/GA_round0", shell=True)

##############################################################################################################
##############################################################################################################

# private function for the mutate function
def get_interaction_energy_mutate(peptide):
    generate_single_individual_list(peptide)
    generate_single_FOLDX_mutate_runfiles("hsc70_Repair_1_0.pdb",0)
    generate_single_FOLDX_analyse_runfiles("hsc70_Repair_1_0.pdb",0)


    path_to_file="/switchlab/group/krinad/hsc70/GA_round0"
    temop_mq="qsub -N mutate " + path_to_file+"/FOLDX_mutate_round0_pdb"+"_original"+"_hsc70.q"
    subprocess.run(temop_mq, shell=True)

    while not os.path.exists("/switchlab/group/krinad/hsc70/GA_round0/hsc70_Repair_1_0_1_2.pdb"):
    	time.sleep(1)

    os.chdir(path_to_file)
    temop_aq="qsub -hold_jid mutate -N analyse " +"FOLDX_analyse_round0_pdb_original_hsc70.q"
    subprocess.run(temop_aq, shell=True)   

    while not os.path.exists("/switchlab/group/krinad/hsc70/GA_round0/Summary_hsc70_Repair_1_0_1_2_AC.fxout"):
        time.sleep(1)

    os.chdir("/switchlab/group/krinad/hsc70/GA_round0")
    command="tail -n1 -q Summary_hsc70_Repair_1_0_1_* >> Summary_file.txt"
    subprocess.run(command, shell=True)

    time.sleep(10)


    temp_ie=pd.read_csv("/switchlab/group/krinad/hsc70/GA_round0/Summary_file.txt", delimiter="\s+")
    temp_ie.columns=["Pdb", "Group1", "Group2","IntraclashesGroup1","IntraclashesGroup2","Interaction_Energy","StabilityGroup1","StabilityGroup2"]
    
    os.chdir("/switchlab/group/krinad/hsc70/")
    return(temp_ie.Interaction_Energy.mean())

# private function for the get_interaction_energy_mutate() function
def generate_single_individual_list(peptide):
    with open('/switchlab/group/krinad/hsc70/GA_round0/individual_list.txt','w',newline='\n') as myfile:
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
    temp_file="/switchlab/group/krinad/hsc70/GA_round"+str(nround)+"/config_mutate_round"+str(nround)+"_pdb"+"_original"+"_hsc70.cfg"
    with open(temp_file,'w',newline='\n') as configfile:
                text_list_config = ["command=BuildModel\n","pdb="+originalPDB+"\n",
                    "pdb-dir=/switchlab/group/krinad/hsc70/\n",
                    #"rotabaseLocation/switchlab/group/krinad/dnak/rotabase.txt\n",
                    "mutant-file=/switchlab/group/krinad/hsc70/GA_round"+str(nround)+"/individual_list.txt\n",
                    "output-dir=/switchlab/group/krinad/hsc70/GA_round"+str(nround)+"\n",
                    "numberOfRuns=3\n"]
                configfile.writelines(text_list_config)
                configfile.close()

    temp_analyse="Summary_file.txt"
    with open(temp_analyse,'w') as helper_file:
        helper_file.write("Pdb    Group1    Group2    IntraclashesGroup1    IntraclashesGroup2    Interaction_Energy    StabilityGroup1    StabilityGroup2\n")
    helper_file.close()
      
    temp_file_q="/switchlab/group/krinad/hsc70/GA_round"+str(nround)+"/FOLDX_mutate_round"+str(nround)+"_pdb"+"_original"+"_hsc70.q"
    with open(temp_file_q,'w',newline='\n') as qfile:
        text_list_q = ["#!/bin/bash\n",
                            "#$ -N FOLDX_hsc70_round"+str(nround)+"_hc\n",
                            "#$ -cwd\n",
                            "#$ -V\n" ,
                            "#$ -q all.q\n" ,
                            "#$ -l h_vmem=700M\n",
                            "\n",
                            "/switchlab/group/tools/2021_FoldX5_LoopX/foldx -f /switchlab/group/krinad/hsc70/GA_round"+str(nround)+"/config_mutate_round"+str(nround)+"_pdb"+"_original"+"_hsc70.cfg"+"\n"]
                            #"mv /switchlab/group/krinad/Summary_file.txt /switchlab/group/krinad/GA_round"+str(nround)+"\n"]
        qfile.writelines(text_list_q)
        qfile.close()





def generate_single_FOLDX_analyse_runfiles(originalPDB, nround):
    temp_file="/switchlab/group/krinad/hsc70/GA_round"+str(nround)+"/config_analyse_round"+str(nround)+"_pdb_"+"original"+"_hsc70.cfg"
    with open(temp_file,'w',newline='\n') as configfile:
        text_list_config = ["command=AnalyseComplex\n","pdb-list=/switchlab/group/krinad/hsc70/GA_round"+str(nround)+"/pdb-list.txt\n",
            "analyseComplexChains=A,C\n",
            "output-dir=/switchlab/group/krinad/hsc70/GA_round"+str(nround)+"\n"]
        configfile.writelines(text_list_config)
        configfile.close()

    temp_file_q="/switchlab/group/krinad/hsc70/GA_round"+str(nround)+"/FOLDX_analyse_round"+str(nround)+"_pdb_"+"original"+"_hsc70.q"
    with open(temp_file_q,'w',newline='\n') as qfile:
        text_list_q = ["#!/bin/bash\n",
                        "#$ -N FOLDX_hsc70_round"+str(nround)+"_ana_hc\n",
                        "#$ -cwd\n",
                        "#$ -V\n" ,
                        "#$ -q all.q\n" ,
                        "#$ -l h_vmem=700M\n",
                        "\n",
                        "/switchlab/group/tools/2021_FoldX5_LoopX/foldx -f /switchlab/group/krinad/hsc70/GA_round"+str(nround)+"/config_analyse_round"+str(nround)+"_pdb_"+"original"+"_hsc70.cfg\n"]
        qfile.writelines(text_list_q)
        qfile.close()
        
    temp_file="/switchlab/group/krinad/hsc70/GA_round"+str(nround)+"/pdb-list.txt"
    with open(temp_file,'w',newline='\n') as pdbfile:
        for i in range(0,3):
            pdbfile.write("/switchlab/group/krinad/hsc70/GA_round"+str(nround)+"/"+"hsc70_Repair_1_0_1_"+str(i)+".pdb;\n")
            #pdbfile.write("/switchlab/group/krinad/GA_round"+str(nround)+"/"+"1dkz_Repair_1_0_1.pdb;\n")
    pdbfile.close()


##############################################################################################################
##############################################################################################################


## get data for the nexxt round of the GA returns a pandas dataframe with the interaction energy and peptideseq
def get_interaction_energy_full(peptides,count):
    create_indlist(peptides,"C",2,7,count)
    generate_FOLDX_mutate_runfiles("hsc70_Repair_1_0.pdb",count)
    generate_FOLDX_analyse_runfiles("hsc70_Repair_1_0.pdb",count)
    run_files(count)
    location_file="/switchlab/group/krinad/hsc70/GA_round"+str(count)+"/analyse/Summary_file.txt"
    temp_dataset=compute_fitness(location_file,peptides,count)
    return(temp_dataset)

    
# generate the mutate files 
def generate_FOLDX_mutate_runfiles(originalPDB,nround):
    temp_analyse="/switchlab/group/krinad/hsc70/GA_round"+str(nround)+"/analyse/Summary_file.txt"
    with open(temp_analyse,'w',newline='\n') as helper_file:
        helper_file.write("Pdb    Group1    Group2    IntraclashesGroup1    IntraclashesGroup2    Interaction_Energy    StabilityGroup1    StabilityGroup2\n")
    helper_file.close()
    print("done")
    
    for i in range(1,101):
        temp_str="R"+str(i)
        path_to_file="./GA_round"+str(nround)+"/R"+str(i)
        completeName = os.path.join(path_to_file, "config_mutate_round"+str(nround)+"_pdb"+"_original"+"_hsc70.cfg")
        temp_file="./GA_round"+str(nround)+"/config_mutate_round"+str(nround)+"_R"+str(i)+"_pdb"+"_original"+"_hsc70.cfg"
        with open(completeName,'w') as configfile:
            text_list_config = ["command=BuildModel\n","pdb="+originalPDB+"\n",
                "pdb-dir=/switchlab/group/krinad/hsc70/GA_round"+str(nround)+"/R"+str(i)+"\n",
                #"rotabaseLocation=/switchlab/group/krinad/GA_round"+str(nround)+"/R"+str(i)+"/rotabase.txt\n",
                "mutant-file=/switchlab/group/krinad/hsc70/GA_round"+str(nround)+"/R"+str(i)+"/individual_list.txt\n",
                "output-dir=/switchlab/group/krinad/hsc70/GA_round"+str(nround)+"/R"+str(i)+"\n",
                "numberOfRuns=3\n"]
            configfile.writelines(text_list_config)
        configfile.close()
                    
        ## create the summary_analuse file for the analyze step
                    
        path_to_file="./GA_round"+str(nround)+"/R"+str(i)
        completeName = os.path.join(path_to_file, "FOLDX_mutate_round"+str(nround)+"_R"+str(i)+"_pdb"+"_original"+"_hsc70.q")
        #temp_list_summary=(["mv /switchlab/group/krinad/hsc70/Summary_file.txt /switchlab/group/krinad/hsc70/GA_round"+str(nround)+"/analyse\n"]+["\n"] * 99)
        #temp_file_q="/home/krinad/thesisfiles/dnak/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"/FOLDX_mutate_round"+str(nround)+"_pdb"+"_original"+"_dnak.q"
        with open(completeName,'w',newline='\n') as qfile:
            text_list_q = ["#!/bin/bash\n",
                                "#$ -N FOLDX_hsc70_round"+str(nround)+"_hc\n",
                                "#$ -cwd\n",
                                "#$ -V\n" ,
                                "#$ -q all.q \n",
                                "#$ -l h_vmem=700M\n",
                                "\n",
                                "cp /switchlab/group/krinad/hsc70/"+originalPDB+" /switchlab/group/krinad/hsc70/GA_round"+str(nround)+"/R"+str(i)+"\n",
                                #"cp /home/krinad/thesisfiles/dnak/rotabase.txt /switchlab/group/krinad/GA_round"+str(nround)+"/R"+str(i)+"\n",
                                "/switchlab/group/tools/2021_FoldX5_LoopX/foldx -f /switchlab/group/krinad/hsc70/GA_round"+str(nround)+"/R"+str(i)+"/config_mutate_round"+str(nround)+"_pdb"+"_original"+"_hsc70.cfg"+"\n"]
                                #temp_list_summary[i-1]]
            qfile.writelines(text_list_q)
        qfile.close()    
    
    
#generates the analyse runfiles     
def generate_FOLDX_analyse_runfiles(originalPDB, nround):
        completeName = os.path.join("/switchlab/group/krinad/hsc70/GA_round"+str(nround), "pdb-list.txt")
        #temp_file="/home/krinad/thesisfiles/dnak/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"/pdb-list.txt"
        with open(completeName,'w',newline='\n') as pdbfile:
            for i in range(1,101):
                for j in range(0,3):
                    pdbfile.write("/switchlab/group/krinad/hsc70/GA_round"+str(nround)+"/hsc70_Repair_1_0_"+str(i)+"_"+str(j)+".pdb;\n")
        pdbfile.close()

        path_to_file="./GA_round"+str(nround)
        completeName = os.path.join(path_to_file, "config_analyse_round"+str(nround)+"_pdb_"+"original"+"_hsc70.cfg")
        #temp_file="./GA_round"+str(nround)+"config_analyse_round"+str(nround)+"_pdb_"+"original"+"_dnak.cfg"
        with open(completeName,'w',newline='\n') as configfile:
            text_list_config = ["command=AnalyseComplex\n","pdb-list=/switchlab/group/krinad/hsc70/GA_round"+str(nround)+"/pdb-list.txt\n",
                "analyseComplexChains=A,C\n",
                "output-dir=/switchlab/group/krinad/hsc70/GA_round"+str(nround)+"/analyse\n"]
            configfile.writelines(text_list_config)
        configfile.close()
    
        completeName1 = os.path.join(path_to_file, "FOLDX_analyse_round"+str(nround)+"_pdb_"+"original"+"_hsc70.q")
        #temp_file_q="/home/krinad/thesisfiles/dnak/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"FOLDX_analyse_round"+str(nround)+"_pdb_"+"original"+"_dnak.q"
        with open(completeName1,'w',newline='\n') as qfile:
            text_list_q = ["#!/bin/bash\n",
                            "#$ -N FOLDX_hsc70_round"+str(nround)+"_ana_hc\n",
                            "#$ -cwd\n"
                            "#$ -V\n" ,
                            "#$ -q all.q\n" ,
                            "#$ -l h_vmem=700M\n",
                            "\n",
                            "/switchlab/group/tools/2021_FoldX5_LoopX/foldx -f /switchlab/group/krinad/hsc70/GA_round"+str(nround)+"/config_analyse_round"+str(nround)+"_pdb_"+"original"+"_hsc70.cfg\n"]
            qfile.writelines(text_list_q)
        qfile.close()
            
    
    

def run_files(nround):
    for i in range(1,101):
        path_to_file="/switchlab/group/krinad/hsc70/GA_round"+str(nround)+"/R"+str(i)
        os.chdir(path_to_file)
        temop_mq="qsub -N mutate"+str(i)+" " + path_to_file+"/FOLDX_mutate_round"+str(nround)+"_R"+str(i)+"_pdb"+"_original"+"_hsc70.q"
        subprocess.run(temop_mq, shell=True)

    filelist=["/switchlab/group/krinad/hsc70/GA_round"+str(nround)+"/R" + str(n)+"/hsc70_Repair_1_0_1_"+str(k)+".pdb" for n in range(1,101) for k in range(0,3)]
    k=True
    while k:
       if all([os.path.isfile(f) for f in filelist]):
          k=False
       else:
          time.sleep(10)
          
          
    l=list(range(1,101))
    chunks = [l[x:x+1] for x in range(0, len(l), 1)]
    for i in range(0,100):
        count=1
        for j in chunks[i]:
            for k in range(0,3):
                location="/switchlab/group/krinad/hsc70/GA_round"+str(nround)+"/R"+str(i+1)+"/hsc70_Repair_1_0_1_"+str(k)+".pdb"
                destination="/switchlab/group/krinad/hsc70/GA_round"+str(nround)+"/hsc70_Repair_1_0_"+str(j)+"_"+str(k)+".pdb"
                command='mv '+location+" "+destination
                subprocess.run(command, shell=True)
                count=count+1  

    path_to_file="/switchlab/group/krinad/hsc70/GA_round"+str(nround)        
    os.chdir(path_to_file)
    temop_aq="qsub -N analyse " +"FOLDX_analyse_round"+str(nround)+"_pdb"+"_original"+"_hsc70.q"
    subprocess.run(temop_aq, shell=True)   

    os.chdir("/switchlab/group/krinad/hsc70/GA_round"+str(nround)+"/analyse")

    while not os.path.exists("Summary_hsc70_Repair_1_0_100_2_AC.fxout"):
        time.sleep(1)
    command="tail -n1 -q Summary_hsc70_Repair_1_0_*.fxout >> Summary_file.txt"
    subprocess.run(command, shell=True)

    os.chdir("/switchlab/group/krinad/hsc70/")
     

    #["/switchlab/group/krinad/GA_round"+str(nround)+"/R" + str(n)+"/1dkz_Repair_1_0_1.pdb" for n in range(1,101)]



# # ##############################################################################################################
# # ##############################################################################################################
# # ##############################################################################################################
# # ##############################################################################################################
 

# # #generate 100 random individuals most prob you should start w like 1000
# # #randomlist=create_seq(100,6)

crossoverProbability=0.8
nvar=6
maxFitness=-12 ### change to alanine
elitismRate=0.6
populationSize=100
mutationRate=0.3


# # #### initialize the 1st random list. We will remove this later but I need this for testing
# # ## 1st INITIALIZATION

randomlist=['RQMSFH',
  'MRYSWF', 'HSAQYG', 'IDTRVM', 'KFPAIQ',
  'REKPSA', 'WMQHIY', 'YAQGVE',
  'QHDKAN', 'QYCPGT',
   'VCMPWR', 'PNHTWG', 'VFENDP', 'WQNLTC', 'KWNGVE', 
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
  'EDSFAG', 'TACNRK', 'RVFPLN', 'NMCWHQ', 'ECDMWP', 'AGVYHP', 'ILMQEN']# #create_idlist(randomlist,"B",2,7) # create the individual list    


# ## EVALUATION
command_list=['mkdir ./GA_round1','mkdir ./GA_round0','mkdir ./GA_round1/analyse']+ ["mkdir ./GA_round1/R"+str(n) for n in range(1,101)]
# 'mkdir ./GA_round1/R1','mkdir ./GA_round1/R2','mkdir ./GA_round1/R3','mkdir ./GA_round1/R4','mkdir ./GA_round1/R5',
# 'mkdir ./GA_round1/analyse',
# 'mkdir ./GA_round0']
for command in command_list:
    subprocess.run(command, shell=True)




subprocess.Popen(command_list, shell=True)
time.sleep(2)

#output-file="TAG"

#"mv /home/krinad/thesisfiles/dnak/buildmodel_AlanineMutagenesis/GA_trial/"+"Summary_file.txt "+"/home/krinad/thesisfiles/dnak/buildmodel_AlanineMutagenesis/GA_trial/GA_round"+str(nround)+"/analyse/Summary\n",

### need to create directory
create_indlist(randomlist,"C",2,7,1)
generate_FOLDX_mutate_runfiles("hsc70_Repair_1_0.pdb",1)
generate_FOLDX_analyse_runfiles("hsc70_Repair_1_0.pdb",1)
run_files(1)


data1=compute_fitness("/switchlab/group/krinad/hsc70/GA_round1/analyse/Summary_file.txt",randomlist,1)   

temp_file_1="compute_fitness"+str(1)+".txt"
# with open(temp_file_1, 'w') as f:
#     f.write(data1.to_string(header = True, index = False))
# f.close()
# print(data1)
# # # ## elitism?
parents=data1.peptideseq
fitness=data1.Interaction_Energy
eliteIndividuals = parents[0: math.ceil(populationSize*elitismRate)]
toSelect = populationSize - math.ceil(populationSize*elitismRate)

# # # # ## SELECCTION
select=selection(data1.peptideseq, data1.Interaction_Energy, toSelect)

# # # ## CROSSOVER
children_new=crossover(select)
# # # # ## MUTATE
mutated_children=mutate(children_new)

new_gen=list(eliteIndividuals)+ list(mutated_children)


command_list=['mkdir GA_round2','mkdir ./GA_round2/analyse']+["mkdir ./GA_round2/R"+str(n) for n in range(1,101)]
for command in command_list:
    subprocess.run(command, shell=True)

temp_file_2="compute_fitness"+str(2)+".txt"
data2=get_interaction_energy_full(new_gen,2)
# with open(temp_file_2, 'w') as f:
#     f.write(data2.to_string(header = True, index = False))
# f.close()

# # # ## CHECK
count=3
while (count<11):
    command_list=["mkdir GA_round"+str(count),"mkdir ./GA_round"+str(count)+"/analyse"]+["mkdir ./GA_round"+str(count)+"/R"+str(n) for n in range(1,101)]
    for command in command_list:
        subprocess.run(command, shell=True)
    parents=data2.peptideseq
    fitness=data2.Interaction_Energy
    eliteIndividuals = parents[0: math.ceil(populationSize*elitismRate)]
    toSelect = populationSize - math.ceil(populationSize*elitismRate)
    select=selection(data2.peptideseq, data2.Interaction_Energy, toSelect)
    children_new=crossover(select)
    print("children "+str(count)+": \n")
    print(children_new)
    mutated_children=mutate(children_new)
    print("mutated children "+str(count)+": \n")
    print(mutated_children)
    print("length of mutated children :"+ str(len(mutated_children)))
    new_gen=list(eliteIndividuals)+ list(mutated_children)
    print("new gen "+str(count)+":\n")
    print(new_gen)
    data2=get_interaction_energy_full(new_gen,count)
    temp_file_x="compute_fitness"+str(count)+".txt"
    with open(temp_file_x, 'w') as f:
        f.write(data2.to_string(header = True, index = False))
    f.close()
    count=count+1

t1 = time.time()
total = t1-t0
print("time: "+str(total))
    

# with open("/switchlab/group/krinad/peptides_bestGA_testing.txt",'w',newline='\n') as best_peptide_file:
#     for i in range(0,len(data2.peptideseq)):
#         if (data2.Interaction_Energy[i]):
#             best_peptide_file.write(" "+data2.peptideseq[i]+" "+str(data2.Interaction_Energy[i])+"\n")
# best_peptide_file.close()

# t1 = time.time()
# total = t1-t0

# print("total time :"+ str(total))


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
#print(find_substate_name_432567("./1dkz_Repair_19_4_2_5_7_20.pdb","432567",7))

### we have options on where to converge: 
    # A gene is said to have converged when 95% of the population share the same value
    # We converge if the average of the interaction energy is ~~14
    

# while(count<95):

# peptides=create_seq(120,6)
    
# data=pd.read_csv('C:\\Users\\user\\Desktop\\GA\\new\\ga10rounds\\optimize\\trysomethingnew\\Summary_4po2_file.txt', sep="\t")
# temp=[]
# for i in range(len(data)) : 
#   temp.append(int(data.iloc[i, 0].split("_")[4]))
# data["split"]=temp
# data=data.sort_values(by='split')
# data = data.reset_index(drop=True)
# data=data[['Pdb','Interaction_Energy','split']]
# data['AvgIE'] = data.groupby('split')['Interaction_Energy'].transform('mean')
# temp_avg=[data.AvgIE[i] for i in range(0, len(data)-2, 3)]

# df = pd.DataFrame(list(zip(list(peptides), temp_avg)),
#                columns =['peptideseq', 'Interaction_Energy'])

# df=df.sort_values(by='Interaction_Energy')
# df = df.reset_index(drop=True)


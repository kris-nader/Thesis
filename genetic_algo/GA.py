# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 10:45:39 2021

@author: Kristen Michelle Nader
"""
import random
import pandas as pd
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
     with open('individual_list.txt','w') as myfile:
            for i in randomlist.keys():
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
    return data

def selection(peptides,interaction_energy,n):
    # peptides: the peptides and the pdb names
    # fitness: interaciton energy
    # n is the number of peptides to be selected
    # perform a roulette selection

    nselected=0
    percentages=[f/sum(abs(interaction_energy)) for f in abs(interaction_energy)]
    selectedIndividuals = []
    while nselected!=n:
        for i in range(0, len(interaction_energy)):
            if nselected==n:
                break
            x=random.random()
            if (x<=percentages[i]):
                selectedIndividuals.append(peptides[i])
                nselected=nselected+1
        
    return selectedIndividuals

def crossover(peptides):
    children = []
    #data.peptideseq
    tmp = len(peptides)%2
    for i in range(0,len(peptides)-tmp,2):
        prob = random.random()
        if prob<=crossoverProbability:
            cuttingPoint=random.randint(0, nvar-1)
            cuttingPoint2=random.randint(0, nvar-1)
            child1= peptides[i][0:cuttingPoint] + peptides[i+1][cuttingPoint:cuttingPoint2] + peptides[i][cuttingPoint2:]
            child2= peptides[i+1][0:cuttingPoint] + peptides[i][cuttingPoint:cuttingPoint2] + peptides[i+1][cuttingPoint2:]
            if(len(child1)==6 and len(child2)==6):
                children.append(child1)
                children.append(child2)
                i=i-1
        else:
            children.append(peptides[i])
            children.append(peptides[i+1])
    if tmp ==1:
        children.append(peptides[-1])
    return(children)


def mutate(peptides):
    mutatedPeptides = []
    for peptide in peptides:
        for i in range(0,len(peptide)):
            temp_peptide = list(peptide)
            temp_peptide[i]= random.sample(aa,1)[0]
            if (computeFitness("".join(list(peptide))) > maxFitness):
               peptide = temp_peptide
        mutatedPeptides.append(peptide)
    return (mutatedPeptides)
    



###################
#generate 100 random individuals
#randomlist=create_seq(100,6)

crossoverProbability=0.8
nvar=6
maxFitness=-12
randomlist={'RQMSFH': 2,
 'MRYSWF': 4, 'HSAQYG': 2, 'IDTRVM': 3, 'KFPAIQ': 3, 'REKPSA': 1, 'WMQHIY': 4, 'YAQGVE': 3, 'QHDKAN': 1, 'QYCPGT': 1, 'VCMPWR': 3, 'PNHTWG': 1, 'VFENDP': 2, 'WQNLTC': 2, 'KWNGVE': 2, 'NTLGFY': 3, 'REPGAQ': 1, 'DCKIHR': 1, 'CKYHDI': 2, 'FLNGCM': 3, 'PIKLGA': 3, 'QNMVKC': 2, 'CRLMDT': 2, 'MKILTE': 3, 'YTCVAP': 3, 'KSNMQH': 1, 'HIERYD': 2, 'DQWCSN': 1, 'LKEIMP': 3, 'GDPSQN': 0, 'RCFMGY': 3, 'RIEYKP': 2, 'ARIMES': 3, 'RATPFV': 3, 'NKYSDA': 2, 'AKHDVF': 3, 'RAPTKY': 2, 'DFEQHC': 1, 'VMCKLG': 3, 'DKCEGM': 1, 'MKVACE': 3, 'QEDTHP': 0, 'RFCHTV': 2, 'EVTLCS': 2, 'DNRPIG': 1, 'LYFTRG': 3, 'NFLPAM': 4,
 'ANMKEV': 3,   'HAFSEG': 2, 'RSCEMA': 2, 'EVHKWF': 3, 'ICMQKD': 2, 'EIYTLR': 3, 'SCDWQY': 2, 'NLWPIC': 3,
 'NWPLCT': 2, 'AIMRPT': 3, 'IPTNDV': 2, 'VKQDER': 1, 'WLSIEK': 3, 'FRVGEY': 3, 'IPGKAW': 3, 'WLTDCQ': 2,
 'GWYHQL': 3, 'FNGVDC': 2, 'MKSQCY': 2, 'FNARVT': 3, 'LSDYGC': 2, 'LYMDQE': 3, 'DYGHCP': 1, 'NCTGVD': 1, 'VWRICL': 4,
 'KITDNG': 1, 'MIPVHF': 4, 'VEAMGF': 4, 'IMTRKQ': 2, 'SNPELY': 2, 'FWYTKD': 3, 'RGHFDK': 1, 'CPMSDF': 2, 'VWYQHA': 4, 'SFHMPL': 3,
 'RILMHS': 3, 'VMSWDF': 4, 'FGKCDI': 2, 'HWAGDI': 3, 'WIAYMV': 6,
 'YMTGHD': 2, 'MHICQK': 2, 'NVSIRC': 2, 'MKDLRP': 2, 'TFYLEN': 3, 'IEVQAT': 3,
 'EDSFAG': 2, 'TACNRK': 1, 'RVFPLN': 3, 'NMCWHQ': 2, 'ECDMWP': 2, 'AGVYHP': 3, 'ILMQEN': 3}
#create_indlist(randomlist,"C",2,7) # create the individual list    

## calculate fitness
data=compute_fitness("C:\\Users\\user\\Desktop\\hsp70\\alanine_mut_buildmodel\\GA_round1\\Summary_GA_file.txt")    
select=selection(data.peptideseq, data.Interaction_Energy, 100)
children=crossover(select)
print(children)
#new_pop=mutate(children)
# for i in select:
#     print(data[data.peptideseq==i])
    


    

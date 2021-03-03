# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 12:14:24 2021

@author: Kristen Michelle Nader
create individual list for foldx alanine --> every other amino acid
input: chain X , start position, end position 
output: individual_list.txt
"""


aa=["A","R","N","D","C","E","Q","G",
    "H","I","L","K","M","F","P","S",
    "T","W","Y","V"]

chain_x="C"
start_position=2
end_position=7
count=0

# myfile=open('individual_list_kristen.txt','w')
with open('individual_list.txt','w') as myfile:
    for j in range(start_position,end_position+1):
        for i in range(len(aa)):
                temp="A"+chain_x+str(j)+aa[i]+";"      
                myfile.write(temp)
                myfile.write("\n")

myfile.close()









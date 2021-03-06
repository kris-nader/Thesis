
library(stringr)


diffs_dataset=read.table(file="C:\\Users\\user\\Desktop\\homologymodel\\repair\\build_model\\mutate_each_position\\Dif_homologymodel_4po2AC_hsc70_Repair_1_0.fxout",header=T,skip=8)
average_dataset=read.table(file="C:\\Users\\user\\Desktop\\homologymodel\\repair\\build_model\\mutate_each_position\\average_dataset.fxout",header=T,skip=8)
individual_list=read.table(file="C:\\Users\\user\\Desktop\\homologymodel\\repair\\build_model\\mutate_each_position\\individual_list.txt",header=F)

diffs_dataset$index=1:360
average_dataset$index=1:120
average_z=c(rep(1,20),rep(2,20),rep(3,20),rep(4,20),rep(5,20),rep(6,20))
diffs_z=c(rep(1,60),rep(2,60),rep(3,60),rep(4,60),rep(5,60),rep(6,60))

## Plot Total energy for average and diff
plot(diffs_dataset$index,diffs_dataset$total_energy,main="total energy scatterplot",xlab="PDBs",ylab="total energy",col=c("black","red","blue","green","orange","yellow")[diffs_z])
identify(diffs_dataset$index,diffs_dataset$total_energy,labels=diffs_dataset$Pdb,plot=TRUE)

plot(average_dataset$index,average_dataset$total_energy,main="total energy scatterplot",xlab="PDBs",ylab="total energy",col=c("black","red","blue","green","orange","yellow")[average_z])
identify(average_dataset$index,average_dataset$total_energy,labels=average_dataset$Pdb,plot=TRUE)
identify(average_dataset$index,average_dataset$total_energy,labels=individual_list$V1,plot=TRUE)


## which structure has the lowest average total energy?
which(average_dataset$total_energy==min(average_dataset$total_energy))
average_dataset[73,1]
# "homologymodel_4po2AC_hsc70_Repair_1_0_73"
individual_list$V1[73] #"AC5M;"


## lowest from the black 1:20
# which has the lowest total energy from the black
black_energy=average_dataset[1:20,]
min_energy_black=min(average_dataset[1:20,3])
index_black=which(average_dataset[3]==min_energy_black)

plot(black_energy$total_energy)
points(index_black,black_energy$total_energy[index_black],col=2)
identify(average_dataset$index[1:20],black_energy$total_energy,labels=average_dataset$Pdb,plot=TRUE)

individual_list$V1[index_black] #"AC2M;"

# which has the lowest total energy from the red
red_energy=average_dataset[21:40,]
min_energy_red=min(average_dataset[21:40,3])
index_red=which(average_dataset$total_energy==min_energy_red)

individual_list$V1[index_red] # "AC3D;"


# which has lowest total energy from the blue
min_energy_blue=min(average_dataset[41:60,3])
index_blue=which(average_dataset$total_energy==min_energy_blue)

individual_list$V1[index_blue] #"AC4L;"

# which has the lowest total energy from the green
min_energy_green=min(average_dataset[61:80,3])
index_green=which(average_dataset$total_energy==min_energy_green)

individual_list$V1[index_green]  #"AC5M;"

# which has the lowest total energy from the orange
min_energy_orange=min(average_dataset[81:100,3])
index_orange=which(average_dataset$total_energy==min_energy_orange)

individual_list$V1[index_orange] #"AC6P;"

# which has the lowest total energy from the yellow
min_energy_yellow=min(average_dataset[101:120,3])
index_yellow=which(average_dataset$total_energy==min_energy_yellow)

individual_list$V1[index_yellow] #  "AC7F;"

average_dataset_full_model=read.table("C:\\Users\\user\\Desktop\\homologymodel\\repair\\build_model\\mutate_each_position\\mutate_individually_Rresults\\Average_homologymodel_4po2AC_hsc70_Repair_1_0.fxout",header=T,skip=8)
average_dataset_full_model$index=121

average_dataset[121,]=average_dataset_full_model[1,]

original_repair_model=read.table("C:\\Users\\user\\Desktop\\homologymodel\\repair\\Average_homologymodel_4po2AC_hsc70_Repair.fxout",header=T,skip=8)
original_repair_model$index=122

average_dataset[122,]=original_repair_model[1,]

averagefull_z=c(rep(1,20),rep(2,20),rep(3,20),rep(4,20),rep(5,20),rep(6,20),rep(7,1),rep(8,1))
plot(average_dataset$index,average_dataset$total_energy,main="total energy scatterplot",xlab="PDBs",ylab="total energy",col=c("black","red","blue","green","orange","yellow","purple","pink")[averagefull_z])
identify(average_dataset$index,average_dataset$total_energy,labels=average_dataset$index,plot=TRUE)





### work on summary data-hsc70

summary_dataset=read.table(file="C:\\Users\\user\\Desktop\\homologymodel\\repair\\build_model\\mutate_each_position\\analysecomplex\\analys_complex_all_pdbs\\Summary_file.txt",header=T)
min(summary_dataset$Interaction_Energy)
real_order=c()
for (i in 1:length(summary_dataset$Pdb)){
  temp=str_split(summary_dataset$Pdb[i], "_")[[1]][7]
  real_order=c(real_order,as.integer(temp))
}
summary_dataset$real_order=real_order

summary_dataset=summary_dataset[order(summary_dataset$real_order),]
  


summary_dataset$index=1:360
summary_z=c(rep(1,60),rep(2,60),rep(3,60),rep(4,60),rep(5,60),rep(6,60))
summary_dataset$summary_z=summary_z

plot(summary_dataset$index,summary_dataset$Interaction_Energy,main="Interaction Energy scatterplot",xlab="PDBs",ylab="Interaction Energy",col=c("black","red","blue","green","orange","yellow")[summary_z])
identify(summary_dataset$index,summary_dataset$Interaction_Energy,labels=summary_dataset$Pdb,plot=TRUE)


## black
black_energy1=summary_dataset[summary_dataset$summary_z==1,]
min_energy_black1=min(black_energy1$Interaction_Energy)
index_black1=which(summary_dataset$Interaction_Energy==min_energy_black1)

index_black1

# red
red_energy1=summary_dataset[summary_dataset$summary_z==2,]
min_energy_red1=min(red_energy1$Interaction_Energy)
index_red1=which(summary_dataset$Interaction_Energy==min_energy_red1)

index_red1

# blue
blue_energy1=summary_dataset[summary_dataset$summary_z==3,]
min_energy_blue1=min(blue_energy1$Interaction_Energy)
index_blue1=which(summary_dataset$Interaction_Energy==min_energy_blue1)

index_blue1

# green

green_energy1=summary_dataset[summary_dataset$summary_z==4,]
min_energy_green1=min(green_energy1$Interaction_Energy)
index_green1=which(summary_dataset$Interaction_Energy==min_energy_green1)

index_green1

# orange
orange_energy1=summary_dataset[summary_dataset$summary_z==5,]
min_energy_orange1=min(orange_energy1$Interaction_Energy)
index_orange1=which(summary_dataset$Interaction_Energy==min_energy_orange1)

index_orange1

# yellow
yellow_energy1=summary_dataset[summary_dataset$summary_z==6,]
min_energy_yellow1=min(yellow_energy1$Interaction_Energy)
index_yellow1=which(summary_dataset$Interaction_Energy==min_energy_yellow1)

index_yellow1



### Get the pdb names so we can map them to the aa residue

summary_dataset[index_black1,1]  #"./homologymodel_4po2AC_hsc70_Repair_1_0_13_2.pdb" AC2M
summary_dataset[index_red1,1]    #"./homologymodel_4po2AC_hsc70_Repair_1_0_38_2.pdb" AC3W
summary_dataset[index_blue1,1]   #"./homologymodel_4po2AC_hsc70_Repair_1_0_51_2.pdb" AC4L
summary_dataset[index_green1,1]  #"./homologymodel_4po2AC_hsc70_Repair_1_0_73_0.pdb" AC6M
summary_dataset[index_orange1,1] #"./homologymodel_4po2AC_hsc70_Repair_1_0_95_2.pdb" AC6P
summary_dataset[index_yellow1,1] #"./homologymodel_4po2AC_hsc70_Repair_1_0_114_0.pdb"AC7F

#### Now I want to check the stats of the important aa
summary_dataset[index_black1,]  
summary_dataset[index_red1,]    
summary_dataset[index_blue1,]   
summary_dataset[index_green1,]  
summary_dataset[index_orange1,] 
summary_dataset[index_yellow1,] 

small_dataset_imp=summary_dataset[c(index_black1,index_red1,index_blue1,index_green1,index_orange1,index_yellow1),]
small_dataset_imp$aa=c("AC2M","AC3W","AC4L","AC6M","AC6P","AC7F")

par(mfrow=(c(1,1)))

plot(small_dataset_imp$summary_z,small_dataset_imp$IntraclashesGroup1,main="IntraclashesGroup1 scatterplot",xlab="PDBs",ylab="IntraclashesGroup1",col=c("black","red","blue","green","orange","yellow")[small_dataset_imp$summary_z])
identify(small_dataset_imp$summary_z,small_dataset_imp$IntraclashesGroup1,labels=small_dataset_imp$aa,plot=TRUE)
dev.copy(png,'InterClashG1.png')
dev.off()

plot(small_dataset_imp$summary_z,small_dataset_imp$IntraclashesGroup2,main="IntraclashesGroup2 scatterplot",xlab="PDBs",ylab="IntraclashesGroup2",col=c("black","red","blue","green","orange","yellow")[small_dataset_imp$summary_z])
identify(small_dataset_imp$summary_z,small_dataset_imp$IntraclashesGroup2,labels=small_dataset_imp$aa,plot=TRUE)
dev.copy(png,'InterClashG2.png')
dev.off()

plot(small_dataset_imp$summary_z,small_dataset_imp$Interaction_Energy,main="Interaction Energy scatterplot",xlab="PDBs",ylab="Interaction Energy",col=c("black","red","blue","green","orange","yellow")[small_dataset_imp$summary_z])
identify(small_dataset_imp$summary_z,small_dataset_imp$Interaction_Energy,labels=small_dataset_imp$aa,plot=TRUE)
dev.copy(png,'InteractionEnergy.png')
dev.off()

plot(small_dataset_imp$summary_z,small_dataset_imp$StabilityGroup1,main="StabilityGroup1 scatterplot",xlab="PDBs",ylab="StabilityGroup1",col=c("black","red","blue","green","orange","yellow")[small_dataset_imp$summary_z])
identify(small_dataset_imp$summary_z,small_dataset_imp$StabilityGroup1,labels=small_dataset_imp$aa,plot=TRUE)
dev.copy(png,'StabilityG1.png')
dev.off()

plot(small_dataset_imp$summary_z,small_dataset_imp$StabilityGroup2,main="StabilityGroup2 scatterplot",xlab="PDBs",ylab="StabilityGroup2",col=c("black","red","blue","green","orange","yellow")[small_dataset_imp$summary_z])
identify(small_dataset_imp$summary_z,small_dataset_imp$StabilityGroup2,labels=small_dataset_imp$aa,plot=TRUE)
dev.copy(png,'StabilityG2.png')
dev.off()




## WORKING on HSP70

summary_hsp70_dataset=read.table(file="C:\\Users\\user\\Desktop\\hsp70\\alanine_mut_buildmodel\\mutate_other_aa\\analyse_complex_all_aa\\Summary_4po2_file.txt",header=T)
real_order_hsp70=c()
for (i in 1:length(summary_hsp70_dataset$Pdb)){
  temp=str_split(summary_hsp70_dataset$Pdb[i], "_")[[1]][5]
  real_order_hsp70=c(real_order_hsp70,as.integer(temp))
}
summary_hsp70_dataset$real_order_hsp70=real_order_hsp70

summary_hsp70_dataset=summary_hsp70_dataset[order(summary_hsp70_dataset$real_order_hsp70),]


summary_hsp70_dataset$index=1:360
summary_hsp70_z=c(rep(1,60),rep(2,60),rep(3,60),rep(4,60),rep(5,60),rep(6,60))
summary_hsp70_dataset$summary_hsp70_z=summary_hsp70_z

plot(summary_hsp70_dataset$index,summary_hsp70_dataset$Interaction_Energy,main="Interaction Energy scatterplot",xlab="PDBs",ylab="Interaction Energy",col=c("black","red","blue","green","orange","yellow")[summary_hsp70_z])
# identify(summary_dataset$index,summary_dataset$Interaction_Energy,labels=summary_dataset$Pdb,plot=TRUE)


## black
black_energy2=summary_hsp70_dataset[summary_hsp70_dataset$summary_hsp70_z==1,]
min_energy_black2=min(black_energy2$Interaction_Energy)
index_black2=which(summary_hsp70_dataset$Interaction_Energy==min_energy_black2)

index_black2

# red
red_energy2=summary_hsp70_dataset[summary_hsp70_dataset$summary_hsp70_z==2,]
min_energy_red2=min(red_energy2$Interaction_Energy)
index_red2=which(summary_hsp70_dataset$Interaction_Energy==min_energy_red2)

index_red2

# blue
blue_energy2=summary_hsp70_dataset[summary_hsp70_dataset$summary_hsp70_z==3,]
min_energy_blue2=min(blue_energy2$Interaction_Energy)
index_blue2=which(summary_hsp70_dataset$Interaction_Energy==min_energy_blue2)

index_blue2

# green

green_energy2=summary_hsp70_dataset[summary_hsp70_dataset$summary_hsp70_z==4,]
min_energy_green2=min(green_energy2$Interaction_Energy)
index_green2=which(summary_hsp70_dataset$Interaction_Energy==min_energy_green2)

index_green1

# orange
orange_energy2=summary_hsp70_dataset[summary_hsp70_dataset$summary_hsp70_z==5,]
min_energy_orange2=min(orange_energy2$Interaction_Energy)
index_orange2=which(summary_hsp70_dataset$Interaction_Energy==min_energy_orange2)

index_orange2

# yellow
yellow_energy2=summary_hsp70_dataset[summary_hsp70_dataset$summary_hsp70_z==6,]
min_energy_yellow2=min(yellow_energy2$Interaction_Energy)
index_yellow2=which(summary_hsp70_dataset$Interaction_Energy==min_energy_yellow2)

index_yellow2


### Get the pdb names so we can map them to the aa residue

summary_hsp70_dataset[index_black2,1]  #"./4po2AC_Repair_1_0_13_0.pdb"
summary_hsp70_dataset[index_red2,1]    #"./4po2AC_Repair_1_0_38_0.pdb"
summary_hsp70_dataset[index_blue2,1]   #"./4po2AC_Repair_1_0_51_1.pdb"
summary_hsp70_dataset[index_green2,1]  #"./4po2AC_Repair_1_0_73_1.pdb"
summary_hsp70_dataset[index_orange2,1] #"./4po2AC_Repair_1_0_95_1.pdb"
summary_hsp70_dataset[index_yellow2,1] #"./4po2AC_Repair_1_0_114_0.pdb"

#### Now I want to check the stats of the important aa
summary_hsp70_dataset[index_black2,]  
summary_hsp70_dataset[index_red2,]    
summary_hsp70_dataset[index_blue2,]   
summary_hsp70_dataset[index_green2,]  
summary_hsp70_dataset[index_orange2,] 
summary_hsp70_dataset[index_yellow2,] 


small_dataset_imp_hsp70=summary_hsp70_dataset[c(index_black2,index_red2,index_blue2,index_green2,index_orange2,index_yellow2),]
small_dataset_imp_hsp70$aa=c("AC2M","AC3W","AC4L","AC6M","AC6P","AC7F")



setwd("C:\\Users\\user\\Desktop\\hsp70\\alanine_mut_buildmodel\\mutate_other_aa\\analyse_complex_all_aa")

par(mfrow=(c(1,1)))

plot(small_dataset_imp_hsp70$summary_hsp70_z,small_dataset_imp_hsp70$IntraclashesGroup1,main="IntraclashesGroup1 scatterplot",xlab="PDBs",ylab="IntraclashesGroup1",col=c("black","red","blue","green","orange","yellow")[small_dataset_imp_hsp70$summary_hsp70_z])
identify(small_dataset_imp_hsp70$summary_hsp70_z,small_dataset_imp_hsp70$IntraclashesGroup1,labels=small_dataset_imp_hsp70$aa,plot=TRUE)
dev.copy(png,'InterClashG1_hsp70.png')
dev.off()

plot(small_dataset_imp_hsp70$summary_hsp70_z,small_dataset_imp_hsp70$IntraclashesGroup2,main="IntraclashesGroup2 scatterplot",xlab="PDBs",ylab="IntraclashesGroup2",col=c("black","red","blue","green","orange","yellow")[small_dataset_imp_hsp70$summary_hsp70_z])
identify(small_dataset_imp_hsp70$summary_hsp70_z,small_dataset_imp_hsp70$IntraclashesGroup2,labels=small_dataset_imp_hsp70$aa,plot=TRUE)
dev.copy(png,'InterClashG2_hsp70.png')
dev.off()

plot(small_dataset_imp_hsp70$summary_hsp70_z,small_dataset_imp_hsp70$Interaction_Energy,main="Interaction Energy scatterplot",xlab="PDBs",ylab="Interaction Energy",col=c("black","red","blue","green","orange","yellow")[small_dataset_imp_hsp70$summary_hsp70_z])
identify(small_dataset_imp_hsp70$summary_hsp70_z,small_dataset_imp_hsp70$Interaction_Energy,labels=small_dataset_imp_hsp70$aa,plot=TRUE)
dev.copy(png,'InteractionEnergy_hsp70.png')
dev.off()

plot(small_dataset_imp_hsp70$summary_hsp70_z,small_dataset_imp_hsp70$StabilityGroup1,main="StabilityGroup1 scatterplot",xlab="PDBs",ylab="StabilityGroup1",col=c("black","red","blue","green","orange","yellow")[small_dataset_imp_hsp70$summary_hsp70_z])
identify(small_dataset_imp_hsp70$summary_hsp70_z,small_dataset_imp_hsp70$StabilityGroup1,labels=small_dataset_imp_hsp70$aa,plot=TRUE)
dev.copy(png,'StabilityG1_hsp70.png')
dev.off()

plot(small_dataset_imp_hsp70$summary_hsp70_z,small_dataset_imp_hsp70$StabilityGroup2,main="StabilityGroup2 scatterplot",xlab="PDBs",ylab="StabilityGroup2",col=c("black","red","blue","green","orange","yellow")[small_dataset_imp_hsp70$summary_hsp70_z])
identify(small_dataset_imp_hsp70$summary_hsp70_z,small_dataset_imp_hsp70$StabilityGroup2,labels=small_dataset_imp_hsp70$aa,plot=TRUE)
dev.copy(png,'StabilityG2_hsp70.png')
dev.off()




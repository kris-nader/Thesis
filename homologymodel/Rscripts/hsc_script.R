
diffs_dataset=read.table(file="C:\\Users\\user\\Desktop\\homologymodel\\repair\\build_model\\mutate_each_position\\Dif_homologymodel_4po2AC_hsc70_Repair_1_0.fxout",header=T,skip=8)
average_dataset=read.table("C:\\Users\\user\\Desktop\\homologymodel\\repair\\build_model\\mutate_each_position\\Average_homologymodel_4po2AC_hsc70_Repair_1_0.fxout",header=T,skip=8)
individual_list=read.table("C:\\Users\\user\\Desktop\\homologymodel\\repair\\build_model\\mutate_each_position\\individual_list.txt",header=F)

diffs_dataset$index=1:360
average_dataset$index=1:120
average_z=c(rep(1,20),rep(2,20),rep(3,20),rep(4,20),rep(5,20),rep(6,20))
diffs_z=c(rep(1,60),rep(2,60),rep(3,60),rep(4,60),rep(5,60),rep(6,60))

## Plot Total energy for average and diff
plot(diffs_dataset$index,diffs_dataset$total_energy,main="total energy scatterplot",xlab="PDBs",ylab="total energy",col=c("black","red","blue","green","orange","yellow")[diffs_z])
identify(diffs_dataset$index,diffs_dataset$total_energy,labels=diffs_dataset$Pdb,plot=TRUE)

plot(average_dataset$index,average_dataset$total_energy,main="total energy scatterplot",xlab="PDBs",ylab="total energy",col=c("black","red","blue","green","orange","yellow")[average_z])
identify(average_dataset$index,average_dataset$total_energy,labels=average_dataset$index,plot=TRUE)


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

# before you start this script, the files need to be edited so that 
# you remove the 1st 8 lines of the 1st fxout and then the 
# 1st 9 lines of the rest because you will be appending to the 1st
# sed -i '1,8d' Interaction_4po2_Repair_1_1_0_AC.fxout
# sed -i '1,9d' Interaction_4po2_Repair_1_* rest of them 
# C:\Users\user\Desktop\fall2020\thesis\thesisfiles\analyseresults_compare\Interaction_4po2_Repair_1_ready_for_analysis.fxout
set_pdbs=read.table(file=file.choose(),header=T,sep="\t")
#C:\Users\user\Desktop\fall2020\thesis\thesisfiles\analyse_complex_for_repaired
repair_ref=read.table(file=file.choose(),header=T,sep="\t")

# this is a list of all pdbs and then a list of all the mutations 
# C:\Users\user\Desktop\fall2020\thesis\thesisfiles\new_muatte\analyzeComplex\pdb_list.txt
# C:\Users\user\Desktop\fall2020\thesis\thesisfiles\new_muatte\individual_list.txt
raw_data=read.table(file=file.choose(),header=F,sep="/")
raw_mutations=read.table(file=file.choose(),header=F)
                    
set_pdbs[365,]=repair_ref[1,]
set_pdbs$index=1:365
y=c(rep(99999,365))
for(i in 1:360){
  temp_parse=as.integer(strsplit(set_pdbs[i,1],"_")[[1]][4])
  if(temp_parse>=1 && temp_parse<=20){
    y[i]=1
  }
  if(temp_parse>20 && temp_parse<=40){
    y[i]=2
  }
  if(temp_parse>40 && temp_parse<=60){
    y[i]=3
  }
  if(temp_parse>60 && temp_parse<=80){
    y[i]=4
  }
  if(temp_parse>80 && temp_parse<=100){
    y[i]=5
  }
  if(temp_parse>100 && temp_parse<=120){
    y[i]=6
  }
  
}
y[361]=7
y[362]=8
y[363]=8
y[364]=8
y[365]=9
set_pdbs$y=y

z=c(rep(1,360),rep(2,1),rep(3,3))
plot(set_pdbs$index,set_pdbs$Interaction.Energy,main="Interaction energy scatterplot",xlab="PDBs",ylab="Interaction energy",col=c("black","green","red")[z])

plot(set_pdbs$index,set_pdbs$Interaction.Energy,main="Interaction energy scatterplot",xlab="PDBs",ylab="Interaction energy",col=c("black","red","blue","green","orange","yellow","cyan","purple","pink")[y],legend=TRUE)

legend("bottomleft",c("1-20", "21-40","41-60","61-80","81-100","101-120","reference_alanine","new pep","repair"),cex=.8,col=c("black","red","blue","green","orange","yellow","cyan","purple","pink"),pch=c(1))
       


reference_interaction=set_pbds[361,]$Interaction.Energy
set_pdbs$centered_energy=set_pdbs$Interaction.Energy+abs(reference_interaction)
plot(set_pdbs$index,set_pdbs$centered_energy,main="Centered Interaction energy scatterplot",xlab="PDBs",ylab="Centered Interaction energy",col=c("black","red","blue","green","orange","yellow","cyan","purple","pink")[y])
legend("bottomleft",c("1-20", "21-40","41-60","61-80","81-100","101-120","reference_alanine","new pep","repair"),cex=.8,col=c("black","red","blue","green","orange","yellow","cyan","purple","pink"),pch=c(1))


# before what we got: AC2M,AC3D,AC4L,AC5M,AC6P,AC7F;
# from which files: 13,24,51,73,95,114
# try and find min of each subset
subset_1=subset(set_pdbs,set_pdbs$y==1)
dim(subset_1)
min_energy_1=min(subset_1$centered_energy)
index_1=which(subset_1$centered_energy==min_energy_1)

a1=substring(subset_1[index_1,1],3)
substring_a1=as.integer(strsplit(a1,"_")[[1]][4])
raw_mutations[substring_a1,]
plot(subset_1$centered_energy,main="Subset 1 centered energy")
points(index_1,subset_1$centered_energy[index_1],col=2)
# add points from the build model
# check the actual index from the dataset
index_t1=which(subset_1$index==76)
index_t2=which(subset_1$index==77)
index_t3=which(subset_1$index==78)
points(index_t1,subset_1$centered_energy[index_t1],col=3)
points(index_t2,subset_1$centered_energy[index_t2],col=3)
points(index_t3,subset_1$centered_energy[index_t3],col=3)

# "AC2I;" hydrophobic

subset_2=subset(set_pdbs,set_pdbs$y==2)
dim(subset_2)
min_energy_2=min(subset_2$centered_energy)
index_2=which(subset_2$centered_energy==min_energy_2)
subset_2[index_2,]

a2=substring(subset_2[index_2,1],3)
substring_a2=as.integer(strsplit(a2,"_")[[1]][4])
raw_mutations[substring_a2,]

plot(subset_2$centered_energy,main="Subset 2 centered energy")
points(index_2,subset_2$centered_energy[index_2],col=2)
# add points from the build model
# check the actual index from the dataset
index_t4=which(subset_2$index==112)
index_t5=which(subset_2$index==113)
index_t6=which(subset_2$index==114)
points(index_t4,subset_2$centered_energy[index_t4],col=3)
points(index_t5,subset_2$centered_energy[index_t5],col=3)
points(index_t6,subset_2$centered_energy[index_t6],col=3)

# "AC3W;" hydrophobic

subset_3=subset(set_pdbs,set_pdbs$y==3)
dim(subset_3)
min_energy_3=min(subset_3$centered_energy)
index_3=which(subset_3$centered_energy==min_energy_3)
subset_3[index_3,]

a3=substring(subset_3[index_3,1],3)
substring_a3=as.integer(strsplit(a3,"_")[[1]][4])
raw_mutations[substring_a3,]

plot(subset_3$centered_energy,main="Subset 3 centered energy")
points(index_3,subset_3$centered_energy[index_3],col=2)
points(31,subset_3$centered_energy[31],col=3)
points(32,subset_3$centered_energy[32],col=4)
subset_3[33,1]
subset_3[31,1]
subset_3[32,1]

# add points from the build model
# check the actual index from the dataset
index_t7=which(subset_3$index==200)
index_t8=which(subset_3$index==201)
index_t9=which(subset_3$index==202)
points(index_t7,subset_3$centered_energy[index_t7],col=3)
points(index_t8,subset_3$centered_energy[index_t8],col=3)
points(index_t9,subset_3$centered_energy[index_t9],col=3)
# "AC4L;" hydrophobic # check more here there are 3 or 4 other points RESOLVED ITS THE 3 NUMBER OF TRIES

subset_4=subset(set_pdbs,set_pdbs$y==4)
dim(subset_4)
min_energy_4=min(subset_4$centered_energy)
index_4=which(subset_4$centered_energy==min_energy_4)
subset_4[index_4,]

a4=substring(subset_4[index_4,1],3)
substring_a4=as.integer(strsplit(a4,"_")[[1]][4])
raw_mutations[substring_a4,]

plot(subset_4$centered_energy,main="Subset 4 centered energy")
points(index_4,subset_4$centered_energy[index_4],col=2)

# "AC5M;" hydrophobic 

subset_5=subset(set_pdbs,set_pdbs$y==5)
dim(subset_5)
min_energy_5=min(subset_5$centered_energy)
index_5=which(subset_5$centered_energy==min_energy_5)
subset_5[index_5,]

a5=substring(subset_5[index_5,1],3)
substring_a5=as.integer(strsplit(a5,"_")[[1]][4])
raw_mutations[substring_a5,]

plot(subset_5$centered_energy,main="Subset 5 centered energy")
points(index_5,subset_5$centered_energy[index_5],col=2)
#"AC6P;" special cases

subset_6=subset(set_pdbs,set_pdbs$y==6)
dim(subset_6)
min_energy_6=min(subset_6$centered_energy)
index_6=which(subset_6$centered_energy==min_energy_6)
subset_6[index_6,]

a6=substring(subset_6[index_6,1],3)
substring_a6=as.integer(strsplit(a6,"_")[[1]][4])
raw_mutations[substring_a6,]

plot(subset_6$centered_energy,main="Subset 6 centered energy")
points(index_6,subset_6$centered_energy[index_6],col=2)
#"AC7F;" hydrophobic

# now we are going to look at the new peptide
#C:\Users\user\Desktop\fall2020\thesis\thesisfiles\try__new_pep_ac\ac\Interaction_4po2_Repair_1_merged_newpep_AC.txt
new_pep=read.table(file=file.choose(),header=T,sep="\t")
new_pep$index=c(366,367,368)
new_pep$y=rep(0,3)
new_pep$centered_energy=new_pep$Interaction.Energy+abs(reference_interaction)
set_pdbs_new=set_pdbs
set_pdbs_new=rbind(set_pdbs_new,new_pep)

# now we do the same plot
plot(set_pdbs$index,set_pdbs$centered_energy,main="Centered Interaction energy scatterplot",xlab="PDBs",ylab="Centered Interaction energy",col=c("black","red","blue","green","orange","yellow","cyan","purple","pink")[y])

plot(set_pdbs_new$index,set_pdbs_new$centered_energy,main=" centered Interaction energy scatterplot",xlab="PDBs",ylab="Centered Interaction energy",col=c("black","red","blue","green","orange","yellow","cyan","purple","pink","teal")[y])
points(366,set_pdbs_new$centered_energy[366],pch=19,col="darkgreen")
points(367,set_pdbs_new$centered_energy[367],pch=19,col="darkgreen")
points(368,set_pdbs_new$centered_energy[368],pch=19,col="darkgreen")
legend("bottomleft",c("1-20", "21-40","41-60","61-80","81-100","101-120","reference_alanine","new_pep_bm","repair","new_pep_ac"),cex=.8,col=c("black","red","blue","green","orange","yellow","cyan","purple","pink","brown"),pch=c(1))

# look at individual effects
for (i in c(seq(4,32),34)){
  temp_colnames=colnames(set_pdbs_new)[i]
  temp_save_name=paste(temp_colnames,".jpeg",sep="")
  jpeg(temp_save_name, width = 736, height = 481)
  plot(set_pdbs_new$index,set_pdbs_new[,i],xlab="PDBs",ylab=temp_colnames,col=c("black","red","blue","green","orange","yellow","cyan","purple","pink","teal")[y])
  points(366,set_pdbs_new[366,i],pch=19,col="darkgreen")
  points(367,set_pdbs_new[367,i],pch=19,col="darkgreen")
  points(368,set_pdbs_new[368,i],pch=19,col="darkgreen")
  dev.off()
  
}

par(mfrow=c(2,3))
for (i in c(seq(4,32),34)){
  temp_colnames=colnames(set_pdbs_new)[i]
  plot(set_pdbs_new$index,set_pdbs_new[,i],xlab="PDBs",ylab=temp_colnames,col=c("black","red","blue","green","orange","yellow","cyan","purple","pink","teal")[y])
  points(366,set_pdbs_new[366,i],pch=19,col="darkgreen")
  points(367,set_pdbs_new[367,i],pch=19,col="darkgreen")
  points(368,set_pdbs_new[368,i],pch=19,col="darkgreen")
}

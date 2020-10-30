# C:\Users\user\Desktop\fall2020\thesis\thesisfiles\new_muatte  average edit
data <- read.table(file=file.choose(), sep = "\t" , header = T)
raw_data=read.table(file=file.choose(),header=F)

data$index=1:120
z=c(rep(1,20),rep(2,20),rep(3,20),rep(4,20),rep(5,20),rep(6,20))
# plot total energy scatterplot
plot(data$index,data$total.energy,main="total energy scatterplot",xlab="PDBs",ylab="total energy",col=c("black","red","blue","green","orange","yellow")[z])
# whats the maximum total energy? which pdb
which(data[3] > 5)
data[78,1]
raw_data[78,1]

# whats the min total energy? which pdb
which(data[3] < -4)
data[70,1]
data[71,1]
data[73,1]
raw_data[70,1]
raw_data[71,1]
raw_data[73,1]
# which has the lowest total energy from the black
black_energy=data[1:20,]
min_energy_black=min(data[1:20,3])
index_black=which(data[3]==min_energy_black)
raw_data[index_black,1]

plot(black_energy$total.energy)
points(index_black,black_energy$total.energy[index_black],col=2)

# which has the lowest total energy from the red
min_energy_red=min(data[21:40,3])
index_red=which(data[3]==min_energy_red)
raw_data[index_red,1]
# which has the lowest total energy from the blue
min_energy_blue=min(data[41:60,3])
index_blue=which(data[3]==min_energy_blue)
raw_data[index_blue,1]
# which has the lowest total energy from the green
min_energy_green=min(data[61:80,3])
index_green=which(data[3]==min_energy_green)
raw_data[index_green,1]
# which has the lowest total energy from the orange
min_energy_orange=min(data[81:100,3])
index_orange=which(data[3]==min_energy_orange)
raw_data[index_orange,1]
# which has the lowest total energy from the yellow
min_energy_yellow=min(data[101:120,3])
index_yellow=which(data[3]==min_energy_yellow)
raw_data[index_yellow,1]

# sequence: AC2M,AC3D,AC4L,AC5M,AC6P,AC7F;
new_data=read.table(file=file.choose(), sep = "\t" , header = T)
new_data[1]="4po2_afterR_analysis_mutateoneindex"
new_data$index=121
data[121,]=new_data[1,]
z[121]=7
plot(data$index,data$total.energy,main="total energy scatterplot",xlab="PDBs",ylab="total energy",col=c("black","red","blue","green","orange","yellow","purple")[z])

# whats the minimum total energ and which pdb
min_energy=min(data[3])
which(data[3]==min_energy)
raw_data[73,1]
# plot backbone H bonds
plot(data$index,data$Backbone.Hbond,main="total energy scatterplot",xlab="PDBs",ylab="total energy",col=c("black","red","blue","green","orange","yellow","purple")[z])
# plot sidechainH bond
plot(data$index,data$Sidechain.Hbond,main="total energy scatterplot",xlab="PDBs",ylab="total energy",col=c("black","red","blue","green","orange","yellow","purple")[z])
# plot vanderwaals
plot(data$index,data$Van.der.Waals,main="total energy scatterplot",xlab="PDBs",ylab="total energy",col=c("black","red","blue","green","orange","yellow","purple")[z])


data[index_black,1]
data[index_red,1]
data[index_blue,1]
data[index_green,1]
data[index_orange,1]
data[index_yellow,1]




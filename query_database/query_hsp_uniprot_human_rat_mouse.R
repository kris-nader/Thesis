# Query Uniprot
# Install packages
install.packages("biomaRt")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("UniProt.ws")
install.packages("Rpdb")
# Load library
library("biomaRt")
library("UniProt.ws")
library("Rpdb")

templateHSP <- read.table("C:\\Users\\user\\Desktop\\thesis_files\\gene_list.txt",header=T,sep="\t")

# check the available datasets
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)

ensembl_human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_rat <- useMart("ensembl", "rnorvegicus_gene_ensembl")
ensembl_mouse <- useMart("ensembl", "mmusculus_gene_ensembl")


## HUMAN
structureData_human <- getBM(mart = ensembl_human,
                             attributes = c("external_gene_name", "uniprot_gn_id", "pdb"),
                             filters = "external_gene_name",
                             values = unique(templateHSP$Gene[which(templateHSP$Organism=="Human")]) )
structureData_human$external_gene_name <- toupper (structureData_human$external_gene_name)

## RAT
structureData_rat <- getBM(mart = ensembl_rat,
                             attributes = c("external_gene_name", "uniprot_gn_id", "pdb"),
                             filters = "external_gene_name",
                             values = unique(templateHSP$Gene[which(templateHSP$Organism=="Rat")]) )
structureData_rat$external_gene_name <- toupper (structureData_rat$external_gene_name)

## MOUSE
structureData_mouse <- getBM(mart = ensembl_mouse,
                           attributes = c("external_gene_name", "uniprot_gn_id", "pdb"),
                           filters = "external_gene_name",
                           values = unique(templateHSP$Gene[which(templateHSP$Organism=="Mouse")]) )
structureData_mouse$external_gene_name <- toupper (structureData_mouse$external_gene_name)


# Create Directory
dir.create("C:\\Users\\user\\Desktop\\thesis_files\\Structures_human_PDB")
dir.create("C:\\Users\\user\\Desktop\\thesis_files\\Structures_rat_PDB")
dir.create("C:\\Users\\user\\Desktop\\thesis_files\\Structures_mouse_PDB")

# QUERY HUMAN
for(human_protein in unique(structureData_human$pdb)){
    file <- paste("C:\\Users\\user\\Desktop\\thesis_files\\Structures_human_PDB\\", human_protein, ".pdb", sep = "")
    if (!file.exists(file)){
      download.file(url = paste("https://files.rcsb.org/download/",human_protein,".pdb", sep=""), destfile = file)
    }
}

# QUERY RAT
for(rat_protein in unique(structureData_rat$pdb)){
  file <- paste("C:\\Users\\user\\Desktop\\thesis_files\\Structures_rat_PDB\\", rat_protein, ".pdb", sep = "")
  if (!file.exists(file)){
    download.file(url = paste("https://files.rcsb.org/download/",rat_protein,".pdb", sep=""), destfile = file)
  }
}

# QUERY MOUSE
for(mouse_protein in unique(structureData_mouse$pdb)){
  file <- paste("C:\\Users\\user\\Desktop\\thesis_files\\Structures_mouse_PDB\\", mouse_protein, ".pdb", sep = "")
  if (!file.exists(file)){
    download.file(url = paste("https://files.rcsb.org/download/",mouse_protein,".pdb", sep=""), destfile = file)
  }
}

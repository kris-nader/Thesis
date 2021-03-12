## sequence logo
require(ggplot2)
require(ggseqlogo)

seq_aa_6=read.table("C:\\Users\\user\\Desktop\\hill_climbing\\peptides_best7.txt",header=F)

ggseqlogo( seq_aa_6$V1, seq_type='aa' )

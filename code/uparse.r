# Import data
otutab <- read.table("data/processed/16s_both_uparse/otutab.txt", comment.char="", header=T, row.names=1)
print_names <- c("CP_Deg", "HI_Deg", "MkB_Deg", "MkP_Deg", "OP_Deg", "CP_Fresh", "MkB_Fresh", "MkP_Fresh", "OP_Fresh", "CP_Sed", "Mk_Sed", "OP_Sed")
names(print_names) <- c("CP_Meta_old", "HI_Meta_Old", "MkB_Meta_Old", "MkP_Meta_Old", "OP_Meta_Old", "C.win1", "MkB.win1", "MkP.win1", "O.win1", "C.sed1", "Mk.sed1", "O.sed1")
print_cols <- c("#a46cb7", "#7aa457", "#cb6a49")[c(1,1,1,1,1,2,2,2,2,3,3,3)]
names(print_cols) <- print_names

metadata = data.frame(Sets=print_names, Type=c("Deg", "Deg", "Deg", "Deg", "Deg", "Fresh", "Fresh", "Fresh", "Fresh", "Sed", "Sed", "Sed"), Location=c("CP", "HI", "MkB", "MkP", "OP", "CP", "MkB", "MkP", "OP", "CP", "Mk", "OP"), Location2=c("CP", "HI", "Mk", "Mk", "OP", "CP", "Mk", "Mk", "OP", "CP", "Mk", "OP"))

# Add the sample pairs together
pair_names <- colnames(otutab)[grepl("1$", colnames(otutab))]
pair_sums <- list()
for(pair_name in pair_names){
    pair_sums[[pair_name]] <- apply(otutab[,grepl(pair_name, colnames(otutab))], 1, sum)
}
pair_sums <- do.call(cbind, pair_sums)

# Combine with 'old' samples into final otutab
f_otutab <- cbind(otutab[,grepl("Meta", colnames(otutab))], pair_sums)
colnames(f_otutab) <- print_names[colnames(f_otutab)]
f_otutab <- f_otutab[,print_names]

# Rarefy table
library(vegan)
rr_otutab <- rrarefy(t(f_otutab), 2500)

# Nice order
nice_order <- c("CP_Fresh", "CP_Deg", "CP_Sed", "MkB_Fresh", "MkB_Deg", "Mk_Sed", "MkP_Fresh", "MkP_Deg", "OP_Fresh", "OP_Deg", "OP_Sed", "HI_Deg")

# Alpha analyses
rich <- specnumber(rr_otutab)
rich <- rich[nice_order]
shannon <- diversity(rr_otutab)
shannon <- shannon[nice_order]
pielou <- shannon/log(rich)
pielou <- pielou[nice_order]

svg("data/figures/alpha_measures.svg", width=18, height=6)
par(mfrow=c(1,3), mar=c(7.1,4.1,4.1,2.1))

barplot(rich, xaxt="n", yaxt="n", ylim=c(0,2000), col=print_cols[names(rich)], space=c(1,0,0,1,0,0,1,0,1,0,0,1), names.arg=NA)
abline(h=1:3*500, lty=2)
par(new=T)
barplot(rich, ylim=c(0,2000), col=print_cols[names(rich)], space=c(1,0,0,1,0,0,1,0,1,0,0,1), las=2, main="Richness")

barplot(shannon, xaxt="n", yaxt="n", ylim=c(0,7), col=print_cols[names(shannon)], space=c(1,0,0,1,0,0,1,0,1,0,0,1), las=2, names.arg=NA)
abline(h=1:6, lty=2)
par(new=T)
barplot(shannon, ylim=c(0,7), col=print_cols[names(shannon)], space=c(1,0,0,1,0,0,1,0,1,0,0,1), las=2, main="Shannon Diversity")

barplot(pielou, xaxt="n", yaxt="n", ylim=c(0,1), col=print_cols[names(pielou)], space=c(1,0,0,1,0,0,1,0,1,0,0,1), las=2, names.arg=NA)
abline(h=1:4*0.2, lty=2)
par(new=T)
barplot(pielou, ylim=c(0,1), col=print_cols[names(pielou)], space=c(1,0,0,1,0,0,1,0,1,0,0,1), las=2, main="Pielou Evenness")

dev.off()

# Beta analyses
shared_species <- apply(rr_otutab, 1, function(x) apply(rr_otutab, 1, function(y) sum((x>0) & (y>0))))
bc_dm <- vegdist(rr_otutab)
nmds <- monoMDS(bc_dm, pc=T, k=3)

# Significance (post-review)
sig = adonis(bc_dm ~ Type + Location, metadata, permutations=1e4)

# NMDS plot
svg("data/figures/nmds.svg", width=6, height=6)
plot(nmds$points[,1:2], pch=21, cex=2, col=1, bg=print_cols[rownames(nmds$points)], xlim=c(-1.1,1.1), ylim=c(-1.1,1.1), sub=paste("Stress:",round(nmds$stress,4)))
text(nmds$points[,1:2], rownames(nmds$points), pos=1)
dev.off()

# Do everything below twice, for 90% and 80% identity to taxonomic db
for(id in c(80,90)){

# Import OTU taxonomy
otutax <- read.table(paste("data/processed/16s_both_uparse/otus_", id, ".lca", sep=""), sep="\t", row.names=1)
otutax <- otutax[rownames(f_otutab),]
colnames(otutax) <- c("Taxonomy", "Identity")
split_otutax <- t(sapply(strsplit(otutax[,1], ";"), function(x) x[1:6]))
rownames(split_otutax) <- rownames(otutax)
colnames(split_otutax) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
split_otutax <- as.data.frame(split_otutax)
split_otutax[is.na(split_otutax)] <- "Unknown"

# Function to summarise table by taxonomic level
summarise_taxon <- function(tab, level){
    count <- aggregate(tab, by=list(split_otutax[rownames(tab),level]), FUN=sum)
    rownames(count) <- count[,1]
    count <- count[,-1]
    pc <- apply(count, 2, function(x) x/sum(x))
    sum <- sort(apply(count, 1, sum), decreasing=T)
    pc <- pc[names(sum), nice_order]

    write.table(pc, paste("data/processed/", level, "_percent_table_", id, ".tsv", sep=""))
    return(pc)
}

# Summarise each sample by phyla
#phyla_count <- aggregate(f_otutab, by=list(split_otutax[rownames(f_otutab),]$Phylum), FUN=sum)
#rownames(phyla_count) <- phyla_count[,1]
#phyla_count <- phyla_count[,-1]
#phyla_pc <- apply(phyla_count, 2, function(x) x/sum(x))
#
#phyla_sum <- sort(apply(phyla_count, 1, sum), decreasing=T)
#phyla_pc <- phyla_pc[names(phyla_sum), nice_order]
#
#write.table(phyla_pc, paste("data/processed/phyla_percent_table_", id, ".tsv", sep=""))

phyla_pc = summarise_taxon(f_otutab, "Phylum")

# Phyla plot
svg(paste("data/figures/phyla_", id, ".svg", sep=""), width=6, height=6)
par(mar=c(7,4,4,10)+0.1, xpd=T)
hili_cols <- c(Firmicutes="#5fba5a",
               Proteobacteria="#ab56c3",
               Bacteroidota="#b1b041",
               Campylobacterota="#6279c9",
               Planctomycetota="#c28541",
               Desulfobacterota="#c080c0",
               Spirochaetota="#5f803f",
               Chloroflexi="#c6517f",
               Fusobacteriota="#4eb8b0",
               Acidobacteriota="#cc5343",
               Unknown="white")
bar_cols <- c(hili_cols[rownames(phyla_pc)[1:10]],gray.colors(nrow(phyla_pc)-10))
barplot(phyla_pc, yaxt="n", beside=F, col=bar_cols, las=2, space=c(1,0.2,0.2,1,0.2,0.2,1,0.2,1,0.2,0.2,1))
axis(2, at=c(0,0.25,0.5,0.75,1), labels=c("0%", "25%", "50%", "75%", "100%"), las=2)
legend("topright", inset=c(-0.65,0), xjust=0, legend=rownames(phyla_pc)[10:1], fill=hili_cols[rownames(phyla_pc)[10:1]])
dev.off()

# Proteobacteria breakdown
proteo_otutab <- f_otutab[split_otutax[rownames(f_otutab),]$Phylum=="Proteobacteria",]
class_pc = summarise_taxon(proteo_otutab, "Class")
order_pc = summarise_taxon(proteo_otutab, "Order")
family_pc = summarise_taxon(proteo_otutab, "Family")
genus_pc = summarise_taxon(proteo_otutab, "Genus")

summary_plot <- function(pc, filename){
    svg(filename, width=6, height=6)
    par(mar=c(7,4,4,10)+0.1, xpd=T)
    bar_cols <- rainbow(nrow(pc))
    barplot(pc, yaxt="n", beside=F, col=bar_cols, las=2, space=c(1,0.2,0.2,1,0.2,0.2,1,0.2,1,0.2,0.2,1))
    axis(2, at=c(0,0.25,0.5,0.75,1), labels=c("0%", "25%", "50%", "75%", "100%"), las=2)
    legend("topright", inset=c(-0.65,0), xjust=0, legend=rownames(pc), fill=bar_cols)
    dev.off()
}

summary_plot(class_pc, paste("data/figures/proteo_class_", id, ".svg", sep=""))
summary_plot(order_pc, paste("data/figures/proteo_order_", id, ".svg", sep=""))
summary_plot(family_pc, paste("data/figures/proteo_family_", id, ".svg", sep=""))
summary_plot(genus_pc, paste("data/figures/proteo_genus_", id, ".svg", sep=""))

}

# Upset plot attempt
library(UpSetR)
b_otutab <- f_otutab
b_otutab[b_otutab>1] <- 1
b_otutab <- b_otutab[,rev(c("CP_Fresh", "CP_Deg", "CP_Sed", "MkB_Fresh", "MkB_Deg", "Mk_Sed", "MkP_Fresh", "MkP_Deg", "OP_Fresh", "OP_Deg", "OP_Sed", "HI_Deg"))]
svg("data/figures/upset.svg", width=18, height=6)
upset(b_otutab, sets=colnames(b_otutab), nintersects=185, order.by=c("degree"), keep.order=T, set.metadata=list(data=metadata, plots=list(list(type="matrix_rows", column="Type", colors=c(Deg="#a46cb7", Fresh="#7aa457", Sed="#cb6a49"), alpha=0.5))))
dev.off()

# Venn diagrams
library(limma)
vennCP <- vennCounts(b_otutab[,grepl("CP", colnames(b_otutab))])
vennCP[1,ncol(vennCP)] <- 0
svg("data/figures/vennCP.svg", width=6, height=6)
vennDiagram(vennCP, circle.col=print_cols[colnames(vennCP)])
dev.off()

vennOP <- vennCounts(b_otutab[,grepl("OP", colnames(b_otutab))])
vennOP[1,ncol(vennOP)] <- 0
svg("data/figures/vennOP.svg", width=6, height=6)
vennDiagram(vennOP, circle.col=print_cols[colnames(vennOP)])
dev.off()

vennMkB <- vennCounts(b_otutab[,c("Mk_Sed", "MkB_Deg", "MkB_Fresh")])
vennMkB[1,ncol(vennMkB)] <- 0
svg("data/figures/vennMkB.svg", width=6, height=6)
vennDiagram(vennMkB, circle.col=print_cols[colnames(vennMkB)])
dev.off()

vennMkP <- vennCounts(b_otutab[,c("Mk_Sed", "MkP_Deg", "MkP_Fresh")])
vennMkP[1,ncol(vennMkP)] <- 0
svg("data/figures/vennMkP.svg", width=6, height=6)
vennDiagram(vennMkP, circle.col=print_cols[colnames(vennMkP)])
dev.off()

vennMk <- vennCounts(b_otutab[,grepl("Mk", colnames(b_otutab))])
vennMk[1,ncol(vennMk)] <- 0
svg("data/figures/vennMk.svg", width=6, height=6)
vennDiagram(vennMk, circle.col=as.vector(print_cols[colnames(vennMk)][c(2:5,1)]))
dev.off()

# Shared taxonomies
sub_otutab <- b_otutab[,grepl("CP", colnames(b_otutab))]
sub_otutab <- cbind(sub_otutab, otutax)
sub_otutab <- sub_otutab[apply(sub_otutab, 1, function(x) sum(as.numeric(x[1:3]))>0),]
sub_otutab <- sub_otutab[order(sub_otutab[,1],sub_otutab[,2],sub_otutab[,3]),]
write.table(sub_otutab, "data/results/CP.txt")

sub_otutab <- b_otutab[,grepl("OP", colnames(b_otutab))]
sub_otutab <- cbind(sub_otutab, otutax)
sub_otutab <- sub_otutab[apply(sub_otutab, 1, function(x) sum(as.numeric(x[1:3]))>0),]
sub_otutab <- sub_otutab[order(sub_otutab[,1],sub_otutab[,2],sub_otutab[,3]),]
write.table(sub_otutab, "data/results/OP.txt")

sub_otutab <- b_otutab[,c("Mk_Sed", "MkB_Deg", "MkB_Fresh")]
sub_otutab <- cbind(sub_otutab, otutax)
sub_otutab <- sub_otutab[apply(sub_otutab, 1, function(x) sum(as.numeric(x[1:3]))>0),]
sub_otutab <- sub_otutab[order(sub_otutab[,1],sub_otutab[,2],sub_otutab[,3]),]
write.table(sub_otutab, "data/results/MkB.txt")

sub_otutab <- b_otutab[,c("Mk_Sed", "MkP_Deg", "MkP_Fresh")]
sub_otutab <- cbind(sub_otutab, otutax)
sub_otutab <- sub_otutab[apply(sub_otutab, 1, function(x) sum(as.numeric(x[1:3]))>0),]
sub_otutab <- sub_otutab[order(sub_otutab[,1],sub_otutab[,2],sub_otutab[,3]),]
write.table(sub_otutab, "data/results/MkP.txt")

sub_otutab <- b_otutab[,grepl("Mk", colnames(b_otutab))]
sub_otutab <- cbind(sub_otutab, otutax)
sub_otutab <- sub_otutab[apply(sub_otutab, 1, function(x) sum(as.numeric(x[1:5]))>0),]
sub_otutab <- sub_otutab[order(sub_otutab[,1],sub_otutab[,2],sub_otutab[,3],sub_otutab[,4],sub_otutab[,5]),]
write.table(sub_otutab, "data/results/Mk.txt")

sub_otutab <- b_otutab[,grepl("HI", colnames(b_otutab)),F]
sub_otutab <- cbind(sub_otutab, otutax)
sub_otutab <- sub_otutab[apply(sub_otutab, 1, function(x) sum(as.numeric(x[1]))>0),]
sub_otutab <- sub_otutab[order(sub_otutab[,1],sub_otutab[,2]),]
write.table(sub_otutab, "data/results/HI.txt")

# Cyanobacteria specific abundances
cyano_tab <- read.table("data/processed/cyanobacteria.tax", sep="\t")
cyano_otus <- cyano_tab[,1]
cyano_otutab <- f_otutab[cyano_otus,]
write.table(cyano_otutab, "data/processed/cyanobacteria_otutab.txt")

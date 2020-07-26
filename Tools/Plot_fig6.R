library(data.table)
library(ggplot2)


##### Fig. 6a

load("../Data/Isogenic_Screenings/screen.results.RData") # dat, dat.crispri.kras

dat.kras <- as.data.table(dat$kras$genelevel)
setnames(dat.kras, c("gene","rank","group","type","method"))
dat.kras[, method:=ifelse(method=="CRISPR-Cas9","CRISPR-KO","CRISPRi")]
dat.kras[, group:=factor(ifelse(group=="KRAS CDE Pos","CDE+","CDE-"), levels=c("CDE+","CDE-"))]
dat.kras[, type:=factor(ifelse(type=="WildType","KRAS WT","KRAS mutant"), levels=c("KRAS WT","KRAS mutant"))]
dat.kras1 <- rbind(dat.kras, dat.crispri.kras)

dat.kras.s1 <- rbind(
  dat.kras1[group=="CDE+", .(group="CDE+",type="KRAS WT",p=sprintf("P=%.3f",wilcox.test(rank[type=="KRAS WT"],rank[type=="KRAS mutant"], paired=TRUE, alternative="less")$p.value)), by=method],
  dat.kras1[group=="CDE-", .(group="CDE-",type="KRAS WT",p=sprintf("P=%.3f",wilcox.test(rank[type=="KRAS WT"],rank[type=="KRAS mutant"], paired=TRUE, alternative="greater")$p.value)), by=method]
) # the `type` here is just a place-holder for the labeling in ggplot() to work properly
dat.kras.s1[, x:=ifelse(method=="CRISPR-KO",1,2)]
dat.kras.s1[, group:=factor(group, levels=c("CDE+","CDE-"))]
dat.kras.s1[, type:=factor(type, levels=c("KRAS WT","KRAS mutant"))]
dat.kras.s1[, y:=6200]

ggplot(dat.kras1, aes(x=method, y=rank, fill=type)) +
  ylab("Rank of Fold-Change") + ylim(0, 6400) +
  facet_wrap(~group, scale="free_y", ncol=2) +
  geom_boxplot() +
  geom_text(data=dat.kras.s1, aes(x=x, y=y, label=p), size=4) +
  theme_classic() +
  scale_fill_manual(values=c("grey","lightcoral")) +
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_text(size=12, angle=30, hjust=1),
    axis.ticks.x=element_blank(),
    axis.title.y=element_text(size=14, vjust=2),
    axis.text.y=element_text(size=10),
    strip.text.x=element_text(size=12),
    legend.title=element_blank(),
    legend.text=element_text(size=12),
    legend.position="bottom"
  )


##### Fig. 6c

dat <- fread("../Data/competition.assay.kras.tsv")
dat <- dat[, .(gene=Gene, rep=Rep, y=`TdTomato%`, t=Timepoint, sgRNA=Sequence)]

tmp <- dat[, .(n=uniqueN(sgRNA)), by=gene][n==1, gene]
tmp1 <- dat[gene %in% tmp, as.list(coef(summary(lm(y ~ t)))["t",c("Estimate","Std. Error","Pr(>|t|)")]), by=gene]
tmp2 <- dat[!gene %in% tmp, as.list(coef(summary(lm(y ~ t + sgRNA)))["t",c("Estimate","Std. Error","Pr(>|t|)")]), by=gene]
res <- rbind(tmp1, tmp2)
setnames(res, c("gene","coef","se","pval"))

tmp <- fread("../Data/competition.assay.kras.top.cde.genes.csv")
tmp <- tmp[Gene.group=="KRAS CDE Pos", Gene]
res1 <- res[gene %in% c(tmp,"NTC")]

res1 <- res1[order(-coef)]
res1[, gene:=factor(gene, levels=gene)]
cs <- c(colorRampPalette(c("red3","grey"))(8), colorRampPalette(c("grey","blue3"))(2))

ggplot(data=res1[gene!="NTC"], aes(x=gene, y=15*coef, fill=gene)) +
  xlab("CDE+ genes") +
  ylab("% KRAS-mutated cells\n Day 15 - Day 0") +
  geom_bar(stat="identity") +
  scale_fill_manual(values=cs) +
  geom_errorbar(aes(ymin=15*(coef-se), ymax=15*(coef+se)), width=.2, position=position_dodge(.9), color="grey20", size=0.4) +
  geom_hline(yintercept=res1[gene=="NTC",15*coef], linetype="dashed", size=0.7, color="black") +
  annotate("text", x=1, y=res1[gene=="NTC",15*coef]-1.1, hjust=0.3, label="Level of NTC", size=4.5, color="black") +
  theme_classic() +
  theme(axis.title.x=element_text(size=14),
        axis.text.x=element_text(size=12, angle=45, hjust=1),
        axis.title.y=element_text(size=14, vjust=2),
        axis.text.y=element_text(size=12),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        legend.position="none")


##### Fig. 6f

dat0 <- fread("..Data/Enache.mutations.tsv")
dat <- dat0[, .(freq.parent=Allele_Fraction[match("PARENTAL",`Parental/Cas9`)], freq.cas9=Allele_Fraction[match("CAS9",`Parental/Cas9`)]), by=.(gene=Hugo_Symbol, cell=Cell_Line_Origin, mut=Canonical_cDNA_Change, mut.type=Canonical_Variant_Classification)]
dat[is.na(freq.parent), freq.parent:=0]
dat[is.na(freq.cas9), freq.cas9:=0]

cells <- fread("../Data/Enache.cell.line.cas9.level.tsv")
setnames(cells, c("cell","cas9","p53.mut"))

res1 <- dat[mut.type!="synonymous" & freq.parent!=1 & !cell %in% cells[cas9<20,cell] & gene=="KRAS"]

res1[, mutation:=paste(cell, mut, sep=":")]
res1 <- res1[order(freq.cas9-freq.parent)]
res1[, mutation:=factor(mutation, levels=mutation)]
res1[, cl:=(freq.cas9-freq.parent)>0]

ggplot() +
  geom_segment(data=res1, aes(x=freq.parent, y=mutation, xend=freq.cas9, yend=mutation, color=cl), size=1, arrow=arrow(length=unit(0.1, "inches")))+
  xlab("Mutant allele frequency") + ylab("Cell line : Mutation") +
  scale_color_manual(values=c("blue2","red2")) +
  theme_classic() +
  theme(axis.title.x=element_text(size=14),
        axis.text.x=element_text(size=12),
        axis.title.y=element_text(size=14, vjust=2),
        axis.text.y=element_text(size=12),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        legend.position="none")


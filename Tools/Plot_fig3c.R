library(data.table)
library(ggplot2)


##### Fig. 3c

dat <- fread("../Data/competition.assay.p53.tsv")
dat <- dat[, .(gene=Gene, rep=Rep, y=`TdTomato`, t=Timepoint, sgRNA=Sequence)]

tmp <- dat[, .(y=y[t==15]-y[t==0]), by=.(gene, sgRNA, rep)]
dat3 <- tmp[, .(y=mean(y), se=sd(y)/sqrt(length(y))), by=.(gene, sgRNA)]
dat3 <- dat3[, .(id=paste0(gene,".",1:(.N)), y, se), by=gene][order(-y)]
dat3[, id:=factor(id, levels=id)]
cs <- c(colorRampPalette(c("red3","grey"))(15), colorRampPalette(c("grey","blue3"))(13))

ggplot(data=dat3[id!="NTC.1"], aes(x=id, y=y, fill=id)) +
  xlab("CDE+ genes (sgRNAs)") +
  ylab("% p53-mutated cells\n Day 15 - Day 0") +
  geom_bar(stat="identity") +
  scale_fill_manual(values=cs) +
  geom_errorbar(aes(ymin=y-se, ymax=y+se), width=.2, position=position_dodge(.9), color="grey20", size=0.4) +
  geom_hline(yintercept=dat3[id=="NTC.1",y], linetype="dashed", size=0.7, color="black") +
  annotate("text", x=22, y=dat3[id=="NTC.1",y]+3, hjust=-0.95, label="Level of NTC", size=4.5, color="black") +
  theme_classic() +
  theme(axis.title.x=element_text(size=16),
        axis.text.x=element_text(size=14, angle=45, hjust=1),
        axis.title.y=element_text(size=16, vjust=2),
        axis.text.y=element_text(size=14),
        legend.title=element_blank(),
        legend.text=element_text(size=14),
        legend.position="none")


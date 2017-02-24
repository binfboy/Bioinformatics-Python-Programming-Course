###################################
#Joshua Rhoades
#77-jrhoades
#13DEC16
###################################
#MAST697 
#Python Programming
###################################
#FINAL PROJECT R-script
#
#Plot three plots
#Combined plots 1+2 and 3+4
#Used Multiplot script to combine three plots
#
###################################
library(ggplot2)

# data - - - - - - - - - - - - - 
d <- read.table("77-Rhoades.Final.OUTPUT.txt",header=TRUE,sep="\t")
xx <- read.table("77-jrhoades-Final-ExpObsCxxC2mer-Output.txt",header=TRUE,sep="\t")

# 1. CpG motif counts vs. %G+C content 
# 2. ApT motif counts vs. %A+T content
#Use ggplot, but use the different aesthetics on two geom_points. 

atcg <- ggplot(d) + 
  geom_point(aes(ATfrac, ATobs, colour = "red"), shape = 3) + 
  geom_point(aes(GCfrac, CGobs, colour = "blue"), shape = 1)+
  labs(x="Complementary Nucleotide Frequency Within Genes", y="Observed Counts of\n ApT and CpG Per Gene")+ 
  theme(axis.line = element_line(colour="black")) +
  ggtitle("Bacillus_thuringiensis_YBT_1518_uid229419 Gene Compostion") +
  labs(colour="Nucleotides")+  scale_colour_discrete(labels=c("%C+G & #CpG", "%A+T & #ApT"))
atcg

# save plot . . . . . . . . . . . . . . . . 
ggsave(atcg, file="RPLOT-77-jrhoades.final.FREQ-vs-MOTIF.pdf",dpi=300)


#3. Expected CpG motif counts vs. Observed CpG motif counts
#4. Expected ApT motif counts vs. Observed ApT motif counts
#same technique as above, but having trouble with regression line(s), so not included....
count <-ggplot(d) +
  geom_point(aes(ATobs, ATexp, colour = "red"), shape = 3) +
  geom_point(aes(CGobs, CGexp, colour = "blue"), shape = 1) + 
  labs(x="Observed Count", y="Expected Count") +
  theme(axis.line = element_line(colour="black")) +
  ggtitle("Bacillus_thuringiensis_YBT_1518_uid229419 Genome Observed vs Expected Motif Counts") +
  labs(colour="Motifs")+  scale_colour_discrete(labels=c("CpG","ApT"))
count

# save plot . . . . . . . . . . . . . . . . 
ggsave(count, file="RPLOT-77-jrhoades.final.MOTIF-PLOT.pdf",dpi=300)

#put a color gradient to make the high observation points pop out a bit more....
CxxCplot <- ggplot(xx, aes(EXPfreq, CxxC)) + geom_point(aes(colour = CxxC)) + 
  scale_colour_gradient(limits=c(0, 9), low="blue", high="red") +
  stat_smooth(method="lm", se=TRUE) + 
  labs(x="Expected Frequency", y="Observed Counts") +
  theme(axis.line = element_line(colour="black"), legend.position="none") + 
  ggtitle("Expected Frequency vs Observed Counts of all possible  2mer Amino Acid Combinations\n within CxxC Motifs in the Bacillus_thuringiensis_YBT_1518_uid229419 Genome") +
  labs(colour="Observed Count")
CxxCplot

# save plot . . . . . . . . . . . . . . . . 
ggsave(CxxCplot, file="RPLOT-77-jrhoades.final.CxxC-PLOT.pdf",dpi=300)


#plot graphs together . . . . . . . . . . 
source("multiplot.R")
multiplot(atcg, count, CxxCplot, cols=1)
png("77-jrhoades-FINAL-MULTIPLOT.png", width = 2500, height=3500, res=300)
multiplot(atcg, count, CxxCplot, cols=1)
dev.off()




























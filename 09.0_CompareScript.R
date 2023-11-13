# code to link up systematic review betas with direct paper betas

# load libraries
library('tidyverse')
library('readxl')
library('ggrepel')
library('Hmisc')

# capture output
date <- Sys.Date()
sink(file=paste0(date,"_CompareScript.out"),split=TRUE)

# read in systematic review
sysrev <- readxl::read_excel( "risk_betas.xlsx" , sheet = 3 )
head(sysrev)

# separate out beta column
sysrev_new <- sysrev %>% separate(`SRR (95% CI)`, c("SRR","95%CI"), sep = "(?=[(])")

# read in my results
mymsresults <- readxl::read_excel("../../../Submission3/20230725_Corbin-DiRECTmetabolomics_ESM2_Tables.xlsx", sheet=3, skip=5, na = c("NA"))
dim(mymsresults)

mynmrresults <- readxl::read_excel("../../../Submission3/20230725_Corbin-DiRECTmetabolomics_ESM2_Tables.xlsx", sheet=2, skip=5, na = c("NA"))
dim(mynmrresults)

# process results and combine
mymsresults <- mymsresults[,c("Biochemical.name", "Source.platform", "HMDB.id", "Allocation.beta","Allocation.p","Allocation.assoc.flag" )]
mynmrresults <- mynmrresults[,c("Biochemical.name", "Source.platform", "Allocation.beta","Allocation.p", "Allocation.assoc.flag")]
mynmrresults$HMDB.id <- NA
myresults <- bind_rows(mymsresults,mynmrresults)
dim(myresults)

# restrict to significant in my analysis
myresults_sig <- myresults[myresults$Allocation.p < 0.05,]
dim(myresults_sig) 

# make nmr lower case letter to start
myresults_sig$Biochemical.name <- str_to_lower(myresults_sig$Biochemical.name)

# order by p then drop duplicates
myresults_sig <- myresults_sig[order(myresults_sig$Allocation.p),]
myresults_sig <- myresults_sig[!duplicated(myresults_sig$Biochemical.name),]

# merge
# match on HMDB
myresults_sig$SRR1 <- log(as.numeric(sysrev_new$SRR[match(myresults_sig$HMDB.id,sysrev_new$`HMDB/ChEBI\r\nID`)]))
myresults_sig$SRR1p <- (as.numeric(sysrev_new$P[match(myresults_sig$HMDB.id,sysrev_new$`HMDB/ChEBI\r\nID`)]))

# match on biochemical name
myresults_sig$SRR2 <- log(as.numeric(sysrev_new$SRR[match(myresults_sig$Biochemical.name,sysrev_new$Compound)]))
myresults_sig$SRR2p <- (as.numeric(sysrev_new$P[match(myresults_sig$Biochemical.name,sysrev_new$Compound)]))

# match on alternative name with capital
myresults_sig$Biochemical.name.upper <- str_to_title(myresults_sig$Biochemical.name)
myresults_sig$SRR3 <- log(as.numeric(sysrev_new$SRR[match(myresults_sig$Biochemical.name.upper,sysrev_new$Compound)]))
myresults_sig$SRR3p <- (as.numeric(sysrev_new$P[match(myresults_sig$Biochemical.name.upper,sysrev_new$Compound)]))

# combine matches
myresults_sig$SRR <- myresults_sig$SRR1 
myresults_sig$p <- myresults_sig$SRR1p 

myresults_sig <- myresults_sig %>% mutate(`SRR` = replace(`SRR`, which(is.na(`SRR1`)),myresults_sig$SRR2[which(is.na(`SRR1`))]))
myresults_sig <- myresults_sig %>% mutate(`p` = replace(`p`, which(is.na(`SRR1p`)),myresults_sig$SRR2p[which(is.na(`SRR1p`))]))

myresults_sig <- myresults_sig %>% mutate(`SRR` = replace(`SRR`, which(is.na(`SRR1`) & is.na(`SRR2`)),myresults_sig$SRR3[which(is.na(`SRR1p`) & is.na(`SRR2p`))]))
myresults_sig <- myresults_sig %>% mutate(`p` = replace(`p`, which(is.na(`SRR1p`) & is.na(`SRR2p`)),myresults_sig$SRR3p[which(is.na(`SRR1p`) & is.na(`SRR2p`))]))


# check matches
print("Total nominally associated in DiRECT after removal of duplicates:")
nrow(myresults_sig)

print("Total matched is:")
length(which(!is.na(myresults_sig$p)))

print("Number associated (after correction) in DiRECT and nominally associated in SR:")
nrow(myresults_sig[myresults_sig$Allocation.assoc.flag == 1 & myresults_sig$p <= 0.05 & !is.na(myresults_sig$p),])

myresults_sig_toplot <- myresults_sig[!is.na(myresults_sig$SRR),] 

# make first letter caps
myresults_sig_toplot$Biochemical.name.for.figure <- str_to_title(myresults_sig_toplot$Biochemical.name)
# manually replace biochem names where upper case added in second part of name
myresults_sig_toplot <- myresults_sig_toplot %>% mutate(`Biochemical.name.for.figure` = replace(`Biochemical.name.for.figure`, `Biochemical.name.for.figure` == "3-Methyl-2-Oxovalerate" , "3-Methyl-2-oxovalerate")) 
myresults_sig_toplot <- myresults_sig_toplot %>% mutate(`Biochemical.name.for.figure` = replace(`Biochemical.name.for.figure`, `Biochemical.name.for.figure` == "3-Methyl-2-Oxobutyrate" , "3-Methyl-2-oxobutyrate")) 
myresults_sig_toplot <- myresults_sig_toplot %>% mutate(`Biochemical.name.for.figure` = replace(`Biochemical.name.for.figure`, `Biochemical.name.for.figure` == "N-Acetylglycine" , "N-acetylglycine")) 
myresults_sig_toplot <- myresults_sig_toplot %>% mutate(`Biochemical.name.for.figure` = replace(`Biochemical.name.for.figure`, `Biochemical.name.for.figure` == "Alpha-Hydroxyisovalerate" , "\u03b1-hydroxyisovalerate")) 

# make plot
pdf(file=paste0(date,"_ComparePlot.pdf"),width = 9, height = 6)
#postscript(file=paste0(date,"_ComparePlot.eps"),width = 9, height = 6)
myplot <- ggplot(myresults_sig_toplot, aes(x=Allocation.beta,y=SRR,label=Biochemical.name.for.figure)) +
  geom_smooth(method="lm",linetype="dashed", color = "black") +
  geom_point(aes(colour=-log10(p))) +
  labs(color = bquote("-log"[10]*"("*italic(p)*")")) +
  geom_text_repel(data = subset(myresults_sig_toplot, p < 0.05 & Allocation.assoc.flag == 1), size = 3) +
  scale_y_continuous(breaks=seq(-0.25,1.0,0.25),limits=c(-0.25,1.0)) +
  scale_x_continuous(breaks=seq(-1,1,0.25),limits=c(-0.75,0.75)) +
  xlab("Difference in metabolite level (normalised SD units) after DiRECT intervention") +
  ylab(bquote(atop("Log"[e] ~ "(relative risk) of type 2 diabetes" ,"per 1-SD increase in metabolite"))) +
  geom_vline(xintercept = 0, lty="dotted") +
  geom_hline(yintercept = 0, lty="dotted") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(myplot)
dev.off()

# evaluate relationship
cor.test(myresults_sig$Allocation.beta, myresults_sig$SRR)
summary(lm(myresults_sig$SRR~myresults_sig$Allocation.beta))
# write out table of matches
write.table(myresults_sig[!is.na(myresults_sig$SRR),c(1:6, 14:15)], "TableS8_srr_comparison.txt", sep="\t", row.names=F, quote=F)

##################
## QUIT SESSION ##

# capture session info##################

print("Session information:")
sessionInfo()

# save output
sink()

# remove data
rm(list = ls())



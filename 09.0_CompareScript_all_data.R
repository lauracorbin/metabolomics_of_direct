# code to link up systematic review betas with direct paper betas

# load libraries
library('tidyverse')
library('readxl')
library('ggrepel')
library('Hmisc')

# capture output
date <- Sys.Date()
sink(file=paste0(date,"_CompareScript_all_data.out"),split=TRUE)

# read in systematic review
sysrev <- readxl::read_excel( "risk_betas.xlsx" , sheet = 3 )
head(sysrev)

# separate out beta column
sysrev_new <- sysrev %>% separate(`SRR (95% CI)`, c("SRR","95%CI"), sep = "(?=[(])")

# read in my results
mymsresults <- readxl::read_excel("../../../Submission3_AAM/20230725_Corbin-DiRECTmetabolomics_ESM2_Tables.xlsx", sheet=3, skip=5, na = c("NA"))
dim(mymsresults)

mynmrresults <- readxl::read_excel("../../../Submission3_AAM/20230725_Corbin-DiRECTmetabolomics_ESM2_Tables.xlsx", sheet=2, skip=5, na = c("NA"))
dim(mynmrresults)

# process results and combine
mymsresults <- mymsresults[,c("Biochemical.name", "Source.platform", "HMDB.id", "Allocation.beta","Allocation.p","Allocation.assoc.flag" )]
mynmrresults <- mynmrresults[,c("Biochemical.name", "Source.platform", "Allocation.beta","Allocation.p", "Allocation.assoc.flag")]
mynmrresults$HMDB.id <- NA
myresults <- bind_rows(mymsresults,mynmrresults)
dim(myresults)

# make nmr lower case letter to start
myresults$Biochemical.name <- str_to_lower(myresults$Biochemical.name)

# order by p then drop duplicates
myresults <- myresults[order(myresults$Allocation.p),]
myresults <- myresults[!duplicated(myresults$Biochemical.name),]

# merge
# match on HMDB
myresults$SRR1 <- log(as.numeric(sysrev_new$SRR[match(myresults$HMDB.id,sysrev_new$`HMDB/ChEBI\r\nID`)]))
myresults$SRR1p <- (as.numeric(sysrev_new$P[match(myresults$HMDB.id,sysrev_new$`HMDB/ChEBI\r\nID`)]))

# match on biochemical name
myresults$SRR2 <- log(as.numeric(sysrev_new$SRR[match(myresults$Biochemical.name,sysrev_new$Compound)]))
myresults$SRR2p <- (as.numeric(sysrev_new$P[match(myresults$Biochemical.name,sysrev_new$Compound)]))

# match on alternative name with capital
myresults$Biochemical.name.upper <- str_to_title(myresults$Biochemical.name)
myresults$SRR3 <- log(as.numeric(sysrev_new$SRR[match(myresults$Biochemical.name.upper,sysrev_new$Compound)]))
myresults$SRR3p <- (as.numeric(sysrev_new$P[match(myresults$Biochemical.name.upper,sysrev_new$Compound)]))

# combine matches
myresults$SRR <- myresults$SRR1 
myresults$p <- myresults$SRR1p 

myresults <- myresults %>% mutate(`SRR` = replace(`SRR`, which(is.na(`SRR1`)),myresults$SRR2[which(is.na(`SRR1`))]))
myresults <- myresults %>% mutate(`p` = replace(`p`, which(is.na(`SRR1p`)),myresults$SRR2p[which(is.na(`SRR1p`))]))

myresults <- myresults %>% mutate(`SRR` = replace(`SRR`, which(is.na(`SRR1`) & is.na(`SRR2`)),myresults$SRR3[which(is.na(`SRR1p`) & is.na(`SRR2p`))]))
myresults <- myresults %>% mutate(`p` = replace(`p`, which(is.na(`SRR1p`) & is.na(`SRR2p`)),myresults$SRR3p[which(is.na(`SRR1p`) & is.na(`SRR2p`))]))

print("Total matched is:")
length(which(!is.na(myresults$p)))

print("Number associated (after correction) in DiRECT and nominally associated in SR:")
nrow(myresults[myresults$Allocation.assoc.flag == 1 & myresults$p <= 0.05 & !is.na(myresults$p),])

myresults_toplot <- myresults[!is.na(myresults$SRR),] 

# make first letter caps
myresults_toplot$Biochemical.name.for.figure <- str_to_title(myresults_toplot$Biochemical.name)
# manually replace biochem names where upper case added in second part of name
myresults_toplot <- myresults_toplot %>% mutate(`Biochemical.name.for.figure` = replace(`Biochemical.name.for.figure`, `Biochemical.name.for.figure` == "3-Methyl-2-Oxovalerate" , "3-Methyl-2-oxovalerate")) 
myresults_toplot <- myresults_toplot %>% mutate(`Biochemical.name.for.figure` = replace(`Biochemical.name.for.figure`, `Biochemical.name.for.figure` == "3-Methyl-2-Oxobutyrate" , "3-Methyl-2-oxobutyrate")) 
myresults_toplot <- myresults_toplot %>% mutate(`Biochemical.name.for.figure` = replace(`Biochemical.name.for.figure`, `Biochemical.name.for.figure` == "N-Acetylglycine" , "N-acetylglycine")) 
myresults_toplot <- myresults_toplot %>% mutate(`Biochemical.name.for.figure` = replace(`Biochemical.name.for.figure`, `Biochemical.name.for.figure` == "Alpha-Hydroxyisovalerate" , "\u03b1-hydroxyisovalerate")) 

# make plot
pdf(file=paste0(date,"_ComparePlot_all_data.pdf"),width = 9, height = 6)
#postscript(file=paste0(date,"_ComparePlot.eps"),width = 9, height = 6)
myplot <- ggplot(myresults_toplot, aes(x=Allocation.beta,y=SRR,label=Biochemical.name.for.figure)) +
  geom_smooth(method="lm",linetype="dashed", color = "black") +
  geom_point(aes(colour=-log10(p))) +
  labs(color = bquote("-log"[10]*"("*italic(p)*")")) +
  geom_text_repel(data = subset(myresults_toplot, p < 0.05 & Allocation.assoc.flag == 1), size = 3) +
  scale_y_continuous(breaks=seq(-0.5,1.0,0.25),limits=c(-0.50,1.0)) +
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
cor.test(myresults$Allocation.beta, myresults$SRR)
summary(lm(myresults$SRR~myresults$Allocation.beta))

# rename SRR column to reflect fact it is log(SRR)
names(myresults)[14] <- "log(SRR)"

# write out table of matches
write.table(myresults[!is.na(myresults$`log(SRR)`),c(1:6, 14:15)], "TableS8_srr_comparison.txt", sep="\t", row.names=F, quote=F)

##################
## QUIT SESSION ##

# capture session info##################

print("Session information:")
sessionInfo()

# save output
sink()

# remove data
rm(list = ls())



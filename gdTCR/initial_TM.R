
# initial look at gdT compiled data set
# testing data format and upload

require(reshape2)
require(plyr)
require(ggplot2)
require(RColorBrewer)
require(dplyr)
require(grid)
require(vegan)
require(knitr)
require(scales)
require(lattice)
require(latticeExtra)
require(ggdendro)
require(ggfortify)
require(erer)

setwd('~/Documents/Dinner/gd TCR/')


# Themes

paper <- theme(text=element_text(size=20), axis.text=element_text(size = 20), axis.title=element_text(size = 20), plot.title=element_text(size = 20), 
               strip.text=element_text(size=20), legend.key.size=unit(1,"cm"))

ppt <- theme(text=element_text(size=30), axis.text=element_text(size = 30), axis.title=element_text(size = 30), plot.title=element_text(size = 30), 
             strip.text=element_text(size=30), legend.key.size=unit(1.5,"cm"))

## folder tree organizes patients and corresponding TCR sequence data
# patient group -> patient -> sequencing tables
# for each patient, we have (IEL,PBL) x (TRD,TRG) tables
# all are from Vd1 sorted gdT

# set plot folder
plot_folder = 'plots_nochallenge/'

# with uncategorized
#patient_groups = c('Control','Active','GFD','Uncategorized')
#patient_ids = list(c(7,13,40,53,110),c(35,46,47,51,81),c(4,9,28,33,41,43),c(1,49,72,78,111,22,3))

# with key 3 subsets
patient_groups = c('Control','Active','GFD')
patient_ids = list(c(7,13,40,53,106,110,144,111),c(35,46,47,51,81,112,143,22),c(4,9,28,33,41,43,113,3))


# given patient group and id, return tcr sequencing data table
read_patient = function(group,id) {
  # tables to retrieve per individual
  tissues = c('IEL','PBL')
  chains = c('TRG','TRD')
  
  # want to add group,id,tissue,chain columns to data table
  # explicitly written to keep types consistent
  patient.df = data.frame(Group=character(0),ID=numeric(0),Tissue=character(0),Chain=character(0),
                          TRV=character(0),CDR3=character(0),TRJ=character(0),Freq=numeric(0),
                          Count=numeric(0),NT=numeric(0),stringsAsFactors=FALSE)
  
  for (i in 1:length(tissues)) {
    for (j in 1:length(chains)) {
      table_file = paste('data/parsed/',group,'/Chicago #',id,
                         '/Chicago #',id,'_',tissues[i],'_Vd1_',chains[j],'.txt',sep='')
      if (file.exists(table_file)) {
        seqs.df = read.table(table_file,header=T,sep='\t',stringsAsFactors=FALSE)
        colnames(seqs.df) = c('TRV','CDR3','TRJ','Freq','Count','NT')
        
        # add columns and merge with current patient data table
        patient.df = rbind(patient.df,cbind(data.frame(Group=group,ID=id,Tissue=tissues[i],
                                                       Chain=chains[j],stringsAsFactors=FALSE),seqs.df))
      }
    }
  }
  return(patient.df)
}

# mostly works, there are some file warnings -- should check
## iterate through all groups and ids and read patient tables
all_patients.df = data.frame(matrix(nrow=0,ncol=10))
for (i in 1:length(patient_groups)) {
  for (j in 1:length(patient_ids[[i]])) {
    all_patients.df = rbind(all_patients.df,read_patient(patient_groups[i],patient_ids[[i]][j]))
  }
}


## Editing all_patients.df but with patient info on sex, age, cell number etc done in excel
## So to generate the new all_patients.df I need to run above code then edit it by hand then reimort it here
## to keep going

# write.csv(all_patients.df, "all_patients.csv")


all_patients.df <- read.csv("data/data_summary.csv")

## slightly modify TRJ labeling => make '2 or 1' equivalent to '1 or 2'
all_patients.df[all_patients.df$TRJ=='2 or 1',]$TRJ = '1 or 2'

# set order of patient group
all_patients.df$Group = factor(all_patients.df$Group,levels=c('Control','Active','GFD','Gluten Challenge Active', 'Uncategorized'))

## finding duplicate chains ##
all_patients.df[duplicated(all_patients.df[,c('Chain','TRV','CDR3','TRJ')]),]

# there are only 4 exact duplicates (matching chain,V,CDR3,J)!
all_patients.df[with(all_patients.df,Chain=='TRD' & TRV=='1' & CDR3=='CALGDQRVPIPWTGGYRHTDKLIF' & TRJ=='1'),]
all_patients.df[with(all_patients.df,Chain=='TRD' & TRV=='1' & CDR3=='CALGEYGRGSWGISHTDKLIF' & TRJ=='1'),]
all_patients.df[with(all_patients.df,Chain=='TRG' & TRV=='2' & CDR3=='CATWDGPNYYKKLF' & TRJ=='2'),]
all_patients.df[with(all_patients.df,Chain=='TRG' & TRV=='9' & CDR3=='CALWEVLYKKLF' & TRJ=='1 or 2'),]

## getting all duplicated CDR3 sequences ##
dup_cdr3_indices = which(duplicated(all_patients.df[,c('Chain','CDR3')]))
dup_cdr3 = list()
for (di in dup_cdr3_indices) {
  dup_cdr3 = c(dup_cdr3,list(all_patients.df[with(all_patients.df,Chain==Chain[di] & CDR3==CDR3[di]),]))
}

## plotting chain distributions by patient ##
trv_summary = ddply(all_patients.df,.(Group,ID,Tissue,Chain,TRV),summarize,count=sum(Count),freq=sum(Freq))

trv_summary = ddply(all_patients.df,.(Group,ID,Tissue,Chain,TRV,Cell_count,Age,Sex,TCRgd_freq,Sort_day,
                                      Cell_count_category,Age_category),summarize,count=sum(Count),freq=sum(Freq))

# create barplot for each combination of tissue and chain
tissues = c('IEL','PBL')
chains = c('TRG','TRD')

# two separate data frames by chain
trv_summary_trd = trv_summary[trv_summary$Chain=='TRD',]
trv_summary_trd$TRV = factor(trv_summary_trd$TRV)
trv_summary_trg = trv_summary[trv_summary$Chain=='TRG',]
trv_summary_trg$TRV = factor(trv_summary_trg$TRV)
trv_summary_split = list(trv_summary_trg,trv_summary_trd)
trv_summary_trg_control <- filter(trv_summary_trg, Group == "Control", Tissue == "PBL" | Tissue == "IEL")

write.csv(trv_summary_trg, "trv_summary_trg.csv")


# looking at top ranked clones in each sample 

topclones <- read.csv("data/Top_ranked_clones.csv")

topclones$Group <- factor(topclones$Group, levels = c('Control', "Active", 'GFD'))

g = ggplot(topclones, aes(factor(Group), Freq))
g = g + geom_boxplot(outlier.size=0)
g = g + geom_point(aes(fill=Group), position = position_jitterdodge(), size = 1.5)
g = g + geom_text(aes(label = Freq, fill=Group), size = 4, hjust = .6, vjust = 1, position = "jitter") 
g = g + facet_grid(Tissue~Chain)
g
ggsave(paste(plot_folder,'Top_clones.png',sep=''),width=15,height=10)

# Looking at if similar numbers of Gammas and Delta clones were sequenced (total transcripts/sample)

g = ggplot(trv_summary, aes(factor(ID), count, fill=Tissue))
g = g + geom_bar(stat="Identity")
#g = g + scale_fill_manual(values=c("TRD"="grey75", "TRG"="grey100"))
g = g + facet_grid(Chain~Group, scales='free_x',space='free_x')
g
ggsave(paste(plot_folder,'Transcripts_Total.png',sep=''),width=15,height=10)


# I wanted to count how many clones per sample so made new column clones
sample_summary_clonecount = ddply(all_patients.df,.(Group,ID,Tissue,Chain),summarize,Clones=n_distinct(CDR3))
sample_summary_clonalsize = ddply(all_patients.df,.(Group,ID,Tissue,Chain),summarize,ClonalSize=n_distinct(CDR3)/sum(Count))

#### gamma to delta clonal size pairwise for controls ###  Figure 2C using clonal size

sample_summary_clonalsize$Tissue <- factor(sample_summary_clonalsize$Tissue, levels = c("PBL", "IEL"))
sample_summary_clonalsize$Chain <- factor(sample_summary_clonalsize$Chain, levels = c("TRG", "TRD"))

g = ggplot(sample_summary_clonalsize[sample_summary_clonalsize$Group=="Control",], aes(factor(Chain), ClonalSize))
g = g + geom_line(aes(group=ID))
g = g + geom_point()
#g = g + geom_boxplot(aes(color = Chain), alpha = 0.1, coef = 0)
g = g + scale_color_manual(values=c("TRG"="grey10", "TRD"="grey50"))
g = g + labs(x ="", y="Clonal Size")
g = g + theme_classic()
g = g + theme(text = element_text(size=20), legend.key.size = unit (1.5, "cm")) + theme(legend.position="none")
g = g + facet_wrap(~Tissue)
g
ggsave(paste(plot_folder,'clonalsize_pairwise.G_D_controls.png',sep=''),width=10,height=5)

#### gamma to delta clonal size pairwise for controls ###  Figure 2C using shannon

divshannon_summary$Tissue <- factor(divshannon_summary$Tissue, levels = c("PBL", "IEL"))
divshannon_summary$Chain <- factor(divshannon_summary$Chain, levels = c("TRG", "TRD"))

g = ggplot(divshannon_summary[divshannon_summary$Group=="Control",], aes(factor(Chain), median))
g = g + geom_line(aes(group=ID))
g = g + geom_point()
#g = g + geom_boxplot(aes(color = Chain), alpha = 0.1, coef = 0)
g = g + scale_color_manual(values=c("TRG"="grey10", "TRD"="grey50"))
g = g + labs(x ="", y="Median Shannnon Diversity")
g = g + theme_classic()
g = g + ppt
g = g + facet_wrap(~Tissue)
g

ggsave(paste(plot_folder,'shannon_pairwise.G_D_controls.pdf',sep=''),width=12,height=8)


# p-value = 0.03876 when you compare Chains irrespective of tissue


#### gamma to delta clonal size pairwise for Active and GFD ###  Figure 4C

sample_summary_clonalsize_A_G <- filter(sample_summary_clonalsize, Group == "Active" | Group == "GFD")
sample_summary_clonalsize_A_G$Tissue <- factor(sample_summary_clonalsize_A_G$Tissue, levels = c("PBL", "IEL"))
sample_summary_clonalsize_A_G$Chain <- factor(sample_summary_clonalsize_A_G$Chain, levels = c("TRG", "TRD"))

g = ggplot(sample_summary_clonalsize_A_G, aes(factor(Chain), ClonalSize))
g = g + geom_line(aes(group=ID))
g = g + geom_point()
#g = g + geom_boxplot(aes(color = Chain), alpha = 0.1, coef = 0)
g = g + scale_color_manual(values=c("TRG"="grey10", "TRD"="grey50"))
g = g + labs(x ="", y="Clonal Size", title="")
g = g + theme_classic()
g = g + theme(panel.background = element_rect(fill="NA", color ="black"))
g = g + theme(text = element_text(size=30), legend.key.size = unit (1.5, "cm")) + theme(legend.position="none")
g = g + facet_grid(Group~Tissue)
g
ggsave(paste(plot_folder,'clonalsize_pairwise_G_D_A_G.png',sep=''),width=13,height=10)

#### gamma to delta shannon pairwise for Active and GFD ###  Figure 4C

divshannon_summary_A_G <- filter(divshannon_summary, Group == "Active" | Group == "GFD")
divshannon_summary_A_G$Tissue <- factor(divshannon_summary_A_G$Tissue, levels = c("PBL", "IEL"))
divshannon_summary_A_G$Chain <- factor(divshannon_summary_A_G$Chain, levels = c("TRG", "TRD"))

g = ggplot(divshannon_summary_A_G, aes(factor(Chain), median))
g = g + geom_line(aes(group=ID))
g = g + geom_point()
#g = g + geom_boxplot(aes(color = Chain), alpha = 0.1, coef = 0)
g = g + scale_color_manual(values=c("TRG"="grey10", "TRD"="grey50"))
g = g + labs(x ="", y="Median Shannon Diversity", title="")
g = g + theme_classic()
g = g + theme(panel.background = element_rect(fill="NA", color ="black"))
g = g + ppt
g = g + facet_grid(Group~Tissue)
g
ggsave(paste(plot_folder,'shannon_pairwise_G_D_A_G.pdf',sep=''),width=12,height=8)

## sampling summary ##
sample_summary = ddply(all_patients.df,.(Group,ID,Tissue,Chain),summarize,samples=sum(Count))
g = ggplot(sample_summary,aes(factor(ID),samples))
g = g + geom_point(aes(color=Tissue,shape=Chain),size=3)
g = g + facet_grid(~Group,scales='free_x',space='free_x')
g = g + xlab('Patient ID')
g = g + ylab('Clones Sequenced')
g = g + scale_color_manual(values=c('IEL'='orange','PBL'='mediumorchid2'))
g = g + theme(text = element_text(size=25))
g
ggsave(paste(plot_folder,'sample_sizes.png',sep=''),width=10,height=5)


# IEL vs. PBL gamma chain usage # supp

Freq_gamma_usage <- read.csv("data/trv_summary_trg.csv")
Freq_gamma_4_usage <- filter(Freq_gamma_usage, TRV == "4")

g = ggplot(Freq_gamma_usage, aes(factor(TRV), freq))
g = g + geom_boxplot()
g = g + geom_point() + geom_jitter(position = position_jitter(width = 0))
g = g + facet_wrap(~Tissue)
g
ggsave(paste(plot_folder,'IEL vs. PBL gamma chain usage.png'), width=12,height=8)

# Control vs. CeD VG4 vs VG3 usage # Figure 

Freq_gamma_usage <- read.csv("data/trv_summary_trg.csv")
Freq_gamma_3_4_usage <- filter(Freq_gamma_usage, TRV == "4" | TRV == "3")

Freq_gamma_3_4_usage$Tissue <- factor(Freq_gamma_3_4_usage$Tissue, levels =c("PBL", "IEL"))
Freq_gamma_3_4_usage$Group_2 <- factor(Freq_gamma_3_4_usage$Group_2, levels =c("Control", "CeD"))


g = ggplot(Freq_gamma_3_4_usage, aes(factor(Tissue), freq, fill = Group_2))
g = g + geom_boxplot(outlier.size = 0)
g = g + geom_point(aes(fill=Group_2), position = position_jitterdodge(), size = 2)
g = g + theme_bw()
#g = g + geom_jitter()
g = g + labs(title="Gamma Usage", x="", y="Frequency of total gammas present")
g = g + scale_fill_manual(values=c("Control"="grey75", "CeD"="grey100"))
g = g + theme(text=element_text(size = 30), legend.key.size = unit (1.5, "cm"))
g = g + facet_wrap(~TRV)
g
ggsave(paste(plot_folder,'Control vs. CeD VG4 vs VG3 usage.png'), width=12,height=8)

# Control vs. CeD Uncharged vs. Positive Charge better # Figure 3C
# Here I made a new file where the freq of pos/neut for each patient is represented so even if its 0%

Freq_gamma_usage_charge_full <- read.csv("data/VG_PosvsNeutral.csv")

Freq_gamma_usage_charge_full_positive <- filter(Freq_gamma_usage_charge_full, Charge=="Positive")

Freq_gamma_usage_charge_full_positive$Tissue <- factor(Freq_gamma_usage_charge_full_positive$Tissue, levels =c("PBL", "IEL"))
Freq_gamma_usage_charge_full_positive$Group <- factor(Freq_gamma_usage_charge_full_positive$Group, levels =c("Control", "Active", "GFD"))

g = ggplot(Freq_gamma_usage_charge_full_positive, aes(factor(Group), Freq, fill = Group))
g = g + geom_boxplot(outlier.size=0)
g = g + geom_point(aes(fill=Group), position = position_jitterdodge(), size = 1.5)
g = g + theme_classic()
g = g + labs(title="Sum of gamma chain CDRs", x="", y="Freq with (+) charges")
g = g + scale_fill_manual(values=c("Control"="grey55", "Active"="grey75", "GFD"="grey95"))
g = g + ppt
g = g + ylim(0, 105)
g = g + theme(legend.title=element_blank())
g = g + facet_wrap(~Tissue)
g

ggsave(paste(plot_folder,'Control vs. CeD Positive Charge.pdf'), width=12,height=8)

stats: 
  balls <- filter(Freq_gamma_usage_charge_full_positive, Tissue=="IEL", Group=="Control" | Group=="GFD")
  wilcox.test(Freq~Group, balls)
  
# Control Charge figure # Figure 2B
Freq_gamma_usage_charge_full <- read.csv("data/VG_PosvsNeutral.csv")

Freq_gamma_usage_charge_full_positive <- filter(Freq_gamma_usage_charge_full, Charge=="Positive")
Freq_gamma_usage_charge_full_positive_Control <- filter(Freq_gamma_usage_charge_full_positive, Group=="Control")

Freq_gamma_usage_charge_full_positive_Control$Tissue <- factor(Freq_gamma_usage_charge_full_positive_Control$Tissue, levels =c("PBL", "IEL"))

g = ggplot(Freq_gamma_usage_charge_full_positive_Control, aes(factor(Tissue), Freq, fill = Tissue))
g = g + geom_boxplot(outlier.size=0)
g = g + geom_point(aes(fill=Tissue), position = position_jitterdodge(), size = 2)
g = g + theme_classic()
g = g + labs(title="Sum of gamma chain CDRs", x="", y="Freq with (+) charges")
g = g + scale_fill_manual(values=c("PBL"="grey70", "IEL"="grey95"))
g = g + ppt
g = g +  theme(legend.title=element_blank())
g = g + ylim(0, 105)
g = g + facet_wrap(~Group)
g

ggsave(paste(plot_folder,'Control Positive Charge.pdf'), width=8,height=8)

Stats:
  wilcox.test(Freq~Tissue, Freq_gamma_usage_charge_full_positive_Control)
  
# Vg usage in PBL vs. IEL # VD1 Paper Figure 3B

vcolors_control <-colorRampPalette(brewer.pal(8,'Set3'))(length(levels(trv_summary_trg_IEL$TRV)))

# IEL
trv_summary_trg_IEL <- filter(trv_summary_trg,Tissue == "IEL", Group=="Control" | Group=="Active" | Group=="GFD")
trv_summary_trg_IEL$Group <- factor(trv_summary_trg_IEL$Group, levels =c("Control", "Active", "GFD"))
trv_summary_trg_IEL$TRV <- factor(trv_summary_trg_IEL$TRV, levels =c("1", "2", "3", "4", "5", "8", "9", "10"))

#PBL
trv_summary_trg_PBL <- filter(trv_summary_trg, Tissue == "PBL", Group == "Control" | Group == "Active" | Group == "GFD")
trv_summary_trg_PBL$Group <- factor(trv_summary_trg_PBL$Group, levels =c("Control", "Active", "GFD"))
trv_summary_trg_PBL$TRV <- factor(trv_summary_trg_PBL$TRV, levels =c("1", "2", "3", "4", "5", "8", "9", "10"))


# IEL_PBL 
trv_summary_trg_IEL_PBL <- filter(trv_summary_trg, Group=="Control" | Group=="Active" | Group=="GFD")
trv_summary_trg_IEL_PBL$Tissue <- factor(trv_summary_trg_IEL_PBL$Tissue, levels =c("PBL", "IEL"))
trv_summary_trg_IEL_PBL$Group <- factor(trv_summary_trg_IEL_PBL$Group, levels =c("Control", "Active", "GFD"))
trv_summary_trg_IEL_PBL$TRV <- factor(trv_summary_trg_IEL_PBL$TRV, levels =c("1", "2", "3", "4", "5", "8", "9", "10"))
trv_summary_trg_IEL_PBL$ID <- factor(trv_summary_trg_IEL_PBL$ID, levels =c("7", "13", "40", "53", "110", "111", "144", "106", "22", "35", "46", "47", "51"
                                                                           , "81", "112", "143", "3", "4", "28", "33", "41", "113", "9", "43"))

g = ggplot(trv_summary_trg_IEL_PBL,aes(factor(ID),freq,fill=TRV))
g = g + geom_bar(stat='identity')
g = g + scale_fill_manual(values=cbbPalette,limits=levels(trv_summary_trg_IEL_PBL$TRV))
g = g + facet_grid(Tissue~Group, scales='free_x',space='free_x')
g = g + xlab('Patient ID')
g = g + ylab('Frequency')
g = g + labs(title="Vgamma Usage")
g = g + theme_classic()
g = g + ppt
#g = g + theme(strip.background = element_blank(),strip.text.x = element_blank())
#g = g + theme(panel.background = element_rect(fill="NA", color ="black"))
g

#IEL
ggsave(paste(plot_folder,'Vg usage in IEL.pdf'),width=14,height=7)
#PBL
ggsave(paste(plot_folder,'Vg usage in PBL.pdf'),width=14,height=7)
#IEL and PBL together
ggsave(paste(plot_folder,'Vg usage in IEL_PBL.pdf'),width=18,height=12)

# Control Vg usage in PBL vs. IEL # VD1 Paper Figure 2A

# To set the order of PBL first before IEL
trv_summary_trg_control$Tissue <- factor(trv_summary_trg_control$Tissue, levels =c("PBL", "IEL"))
trv_summary_trg_control$TRV <- factor(trv_summary_trg_control$TRV, levels =c("1", "2", "3", "4", "5", "8", "9", "10"))

vcolors_control <-colorRampPalette(brewer.pal(8,'Set3'))(length(levels(trv_summary_trg_control$TRV)))

cbbPalette <- c("#FFFFE0", "#B0E0E6", "#F08080", "#BEBEBE", "#EED2EE", "#FFF68F", "#7F7F7F", "#90EE90")

g = ggplot(trv_summary_trg_control,aes(factor(ID),freq,fill=TRV))
g = g + geom_bar(stat='identity')
g = g + scale_fill_manual(values=cbbPalette,limits=levels(trv_summary_trg_control$TRV))
g = g + facet_grid(~Tissue, scales='free_x',space='free_x')
g = g + xlab('Patient ID')
g = g + ylab('Frequency')
g = g + labs(title="Vgamma Usage")
g = g + theme_classic()
g = g + theme(panel.background = element_rect(fill="NA", color ="black"))
g = g + ppt
g

ggsave(paste(plot_folder,'Control Vg usage in PBL vs. IEL.pdf'),width=13,height=8.66)


# Control Vg 4 vs. 3_8_9 in IEL vs. PBL # Figure 2

trv_summary_trg_control_4_3_8_9 <- read.csv("data/trv_summary_trg_control_4_3_8_9.csv")
trv_summary_trg_control_4_3_8_9$Tissue <- factor(trv_summary_trg_control_4_3_8_9$Tissue, levels =c("PBL", "IEL"))


g = ggplot(trv_summary_trg_control_4_3_8_9, aes(factor(Tissue), Frequency))
#g = g + geom_line(aes(group=ID), position = "jitter")
g = g + geom_boxplot(aes(fill = Tissue), outlier.size = 0)
g = g + geom_point(aes(fill=Tissue), position = position_jitterdodge(), size = 2)
g = g + scale_fill_manual(values=c("PBL"="grey75", "IEL"="grey100"))
g = g + labs(x ="", y="Frequency of total gammas present")
g = g + theme_bw()
g = g + theme(text = element_text(size=20), legend.key.size = unit (1.5, "cm")) + theme(legend.position="none")
g = g + facet_wrap(~Chain_Group)
g
ggsave(paste(plot_folder,'Control_VG-3_4_8_9.png',sep=''),width=10,height=5)


# For Figure 3 

for (i in 1:length(tissues)) {
  for (j in 1:length(chains)) {
    trv = trv_summary_split[[j]][with(trv_summary_split[[j]],Tissue==tissues[i] & Chain==chains[j]),]
    
    if (j == 1) {
      vcolors = colorRampPalette(brewer.pal(8,'Accent'))(length(levels(trv$TRV)))
    }
    else {
      vcolors = colorRampPalette(brewer.pal(8,'Set1'))(length(levels(trv$TRV)))
    }
    
    # grouped barplot showing proportions
    g = ggplot(trv,aes(factor(ID),freq,fill=TRV))
    g = g + geom_bar(stat='identity')
    g = g + scale_fill_manual(values=vcolors,limits=levels(trv$TRV))
    g = g + facet_grid(~Group, scales='free_x',space='free_x')
    g = g + xlab('Patient ID')
    g = g + ylab('Frequency')
    g = g + theme(text = element_text(size=20), legend.key.size = unit (1, "cm"))
    g = g + theme(strip.background = element_rect(fill="lavender"))
    g

    ggsave(paste(plot_folder,'TRV_freq_',tissues[i],'_',chains[j],'.png',sep=''),width=10,height=5)
  }
}

## barplot all clones ##
clone_summary = all_patients.df[sample(dim(all_patients.df)[1]),]
clone_summary$Clone = rownames(clone_summary)

# create barplot for each combination of tissue and chain
tissues = c('IEL','PBL')
chains = c('TRG','TRD')

# Control PBL vs. IEL Gamma vs. Delta CDR3 freq bar plots # Figure 2C

clone_summary_control <- filter(clone_summary, Group == "Control")

clone_summary_control$Tissue <- factor(clone_summary_control$Tissue, levels =c("PBL", "IEL"))
clone_summary_control$Chain <- factor(clone_summary_control$Chain, levels =c("TRG", "TRD"))


clonecolors_control = sample(colorRampPalette(brewer.pal(8,'Greys'))(dim(clone_summary_control)[1]))


g = ggplot(clone_summary_control,aes(factor(ID),Freq,fill=Clone))
g = g + geom_bar(stat='identity')
g = g + scale_fill_manual(values=clonecolors_control)
g = g + facet_grid(Chain ~ Tissue,scales='free_x',space='free_x')
g = g + theme_classic()
g = g + theme(panel.background = element_rect(fill="NA", color ="black"))
g = g + theme(legend.position="none")
g = g + xlab('Patient ID')
g = g + ylab('Frequency')
g = g + labs(title="Control CDR3s")
g = g + ppt
g

ggsave(paste(plot_folder,'Control PBL vs. IEL Gamma vs. Delta CDR3 freq bar plots.pdf'),width=12,height=10)

# Active vs. GFD for IEL and PBL separatley compairing Gamma vs. Delta CDR3 freq bar plots # Figure 4A

clone_summary_Active_GFD_IEL <- filter(clone_summary, Tissue == "IEL", Group == "Active" | Group == "GFD" | Group == "Control")
clone_summary_Active_GFD_IEL$Group <- factor(clone_summary_Active_GFD_IEL$Group, levels =c("Control", "Active", "GFD"))

clone_summary_Active_GFD_PBL <- filter(clone_summary, Tissue == "PBL", Group == "Active" | Group == "GFD" | Group == "Control")
clone_summary_Active_GFD_PBL$Group <- factor(clone_summary_Active_GFD_PBL$Group, levels =c("Control", "Active", "GFD"))

clone_summary_Active_GFD_IEL$Chain <- factor(clone_summary_Active_GFD_IEL$Chain, levels =c("TRG", "TRD"))
clone_summary_Active_GFD_PBL$Chain <- factor(clone_summary_Active_GFD_PBL$Chain, levels =c("TRG", "TRD"))


clonecolors_Active_GFD_PBL = sample(colorRampPalette(brewer.pal(8,'Greys'))(dim(clone_summary_Active_GFD_PBL)[1]))
clonecolors_Active_GFD_IEL = sample(colorRampPalette(brewer.pal(8,'Greys'))(dim(clone_summary_Active_GFD_IEL)[1]))


g = ggplot(clonecolors_Active_GFD_IEL,aes(factor(ID),Freq,fill=Clone))
g = g + geom_bar(stat='identity')
g = g + scale_fill_manual(values=clonecolors_Active_GFD_IEL)
g = g + facet_grid(Chain ~ Group,scales='free_x',space='free_x')
g = g + theme_classic()
g = g + theme(panel.background = element_rect(fill="NA", color ="black"))
g = g + theme(legend.position="none")
g = g + xlab('Patient ID')
g = g + ylab('Frequency')
g = g + labs(title="IEL CDR3s")
g = g + ppt
g

ggsave(paste(plot_folder,'Active vs. GFD IEL Gamma vs. Delta CDR3 freq bar plots.pdf'),width=15,height=10)

ggsave(paste(plot_folder,'Active vs. GFD PBL Gamma vs. Delta CDR3 freq bar plots.pdf'),width=15,height=10)

for (i in 1:length(tissues)) {
  for (j in 1:length(chains)) {
    clone = clone_summary[with(clone_summary,Tissue==tissues[i] & Chain==chains[j]),]
    
    clonecolors = sample(colorRampPalette(brewer.pal(8,'Set2'))(dim(clone)[1]))
    # grouped barplot showing proportions
    g = ggplot(clone,aes(factor(ID),Freq,fill=Clone))
    g = g + geom_bar(stat='identity')
    g = g + scale_fill_manual(values=clonecolors)
    g = g + facet_grid(~Group,scales='free_x',space='free_x')
    g = g + theme(legend.position="none")
    g = g + xlab('Patient ID')
    g = g + ylab('Frequency')
    g = g + theme(text = element_text(size=25))
    g
    
    ggsave(paste(plot_folder,'clone_freq_',tissues[i],'_',chains[j],'.png',sep=''),width=10,height=5)
  }
}

## calculate shannon diversity from cdr3 frequencies ##
require(vegan)

diversity_summary = ddply(all_patients.df,.(Group,ID,Tissue,TRV),summarize,
                    shannon=diversity(Count,index='shannon'),simpson=diversity(Count,index='simpson'))
diversity_summary_control_IEL <- filter(diversity_summary, Group == 'Control', Tissue == 'IEL')
# pairwise comparison by chain shannon

g = ggplot(diversity_summary_control_IEL,aes(Chain,shannon,color=Chain))
g = g + geom_point() 
g = g + geom_boxplot()
g = g + geom_line(aes(group=ID))
g = g + scale_color_manual(values=c('TRD'='deepskyblue','TRG'='tomato2'))
g
ggsave(paste(plot_folder,'pairwise_shannon_bychain.png'),width=10,height=5)

# pairwise comparison by chain simpson

g = ggplot(diversity_summary,aes(Group,shannon,color=Tissue))
#g = g + geom_point() 
g = g + geom_boxplot() # + geom_jitter(position=position_jitter(width=0.5, height=1))
#g = g + geom_line(aes(group=ID))
#g = g + scale_color_manual(values=c('Control'='deepskyblue','Active'='tomato2','GFD'='lawngreen'))
g = g + scale_color_manual(values=c('IEL' = 'orange', 'PBL' = 'red'))
g
ggsave(paste(plot_folder,'pairwise_simpson_bychain.png'),width=10,height=5)

# scatter comparison by patient and group
## using shannon
g = ggplot(diversity_summary,aes(factor(ID),shannon))
g = g + geom_point(aes(color=Tissue,shape=Chain),size=3)
g = g + facet_grid(~Group,scales='free_x',space='free_x')
g = g + xlab('Patient ID')
g = g + ylab('Shannon Diversity Index')
g = g + scale_color_manual(values=c('IEL'='orange','PBL'='mediumorchid2'))
g = g + theme(text = element_text(size=25))
g
ggsave(paste(plot_folder,'shannon_summary.png'),width=10,height=5)

shannon_tissue = dcast(diversity_summary,Group + ID + Chain ~ Tissue,value.var='shannon')
g = ggplot(shannon_tissue,aes(IEL,PBL))
g = g + geom_point(aes(color=Group,shape=Chain),size=3)
g = g + geom_abline(a=1)
g = g + scale_color_manual(values=c('Control'='deepskyblue','Active'='tomato2','GFD'='lawngreen'))
g = g + theme(text = element_text(size=25))
g
ggsave(paste('plots_nochallenge/shannon_tissue.png',sep=''),width=10,height=5)

shannon_chain = dcast(diversity_summary,Group + ID + Tissue ~ Chain,value.var='shannon')
g = ggplot(shannon_chain,aes(TRG,TRD))
g = g + geom_point(aes(color=Group,shape=Tissue),size=3)
g = g + geom_abline(a=1)
g = g + scale_color_manual(values=c('Control'='deepskyblue','Active'='tomato2','GFD'='lawngreen'))
g = g + theme(text = element_text(size=25))
g
ggsave(paste(plot_folder,'shannon_chain.png',sep=''),width=10,height=5)

## using simpson
g = ggplot(diversity_summary,aes(factor(ID),simpson))
g = g + geom_point(aes(color=Tissue,shape=Chain),size=3)
g = g + facet_grid(~Group,scales='free_x',space='free_x')
g = g + xlab('Patient ID')
g = g + ylab('Simpson Diversity Index')
g = g + scale_color_manual(values=c('IEL'='orange','PBL'='mediumorchid2'))
g = g + theme(text = element_text(size=25))
g
ggsave(paste(plot_folder,'simpson_summary.png',sep=''),width=10,height=5)

simpson_tissue = dcast(diversity_summary,Group + ID + Chain ~ Tissue,value.var='simpson')
g = ggplot(simpson_tissue,aes(IEL,PBL))
g = g + geom_point(aes(color=Group,shape=Chain),size=3)
g = g + geom_abline(a=1)
g = g + scale_color_manual(values=c('Control'='deepskyblue','Active'='tomato2','GFD'='lawngreen'))
g = g + theme(text = element_text(size=25))
g
ggsave(paste(plot_folder,'simpson_tissue.png',sep=''),width=10,height=5)

simpson_chain = dcast(diversity_summary,Group + ID + Tissue ~ Chain,value.var='simpson')
g = ggplot(simpson_chain,aes(TRG,TRD))
g = g + geom_point(aes(color=Group,shape=Tissue),size=3)
g = g + geom_abline(a=1)
g = g + scale_color_manual(values=c('Control'='deepskyblue','Active'='tomato2','GFD'='lawngreen'))
g = g + theme(text = element_text(size=25))
g
ggsave(paste(plot_folder,'simpson_chain.png',sep=''),width=10,height=5)

## DIVERSITY SECTION

### rarefaction of diversities ###
divsample = function(counts.df,subsamples,nsamples,divindex) {
  counts_expand.df = counts.df[rep(1:nrow(counts.df),counts.df[['Count']]),]
  if (nrow(counts_expand.df) <= subsamples){
    # return regular diversity if sample size is smaller than subsample size
    divresults = diversity(counts.df$Count,index=divindex)
  } else {
    divresults = rep(NA,nsamples)
    for (i in 1:nsamples) {
      subcounts_expand.df = counts_expand.df[sample.int(nrow(counts_expand.df),subsamples),]
      subcounts.df = ddply(subcounts_expand.df,.(Group,ID,Tissue,Chain,TRV,CDR3,TRJ,Cell_count,Age,Sex,TCRgd_freq,Sort_day,
                                               Cell_count_category,Age_category),
                           summarize,Count=length(Count))
      divresults[i] = diversity(subcounts.df$Count,index=divindex)
    }
  }
  return(data.frame(median=median(divresults),lq=quantile(divresults,0.25),
                    uq=quantile(divresults,0.75)))
}

## summarize diversity resampling for Vgamma combined ##
divshannon_summary = ddply(all_patients.df,.(Group,ID,Tissue,Chain,Cell_count,Age,Sex,TCRgd_freq,Sort_day,
                                    Cell_count_category,Age_category),function(x){return(divsample(x,50,100,'shannon'))})

divsimpson_summary = ddply(all_patients.df,.(Group,ID,Tissue,Chain,Cell_count,Age,Sex,TCRgd_freq,Sort_day,
                                             Cell_count_category,Age_category),
                           function(x){return(divsample(x,50,100,'simpson'))})

divshannon_summary_1 = ddply(all_patients.df,.(Group,ID,Tissue,Chain,Cell_count,Age,Sex,TCRgd_freq,Sort_day,
                                               Cell_count_category,Age_category),
                           function(x){return(divsample(x,50,1000,'shannon'))})

write.csv(divshannon_summary, "divshannon_summary.csv")

## diversity for individual gammas 
divshannon_summary_trv = ddply(all_patients.df,.(Group,ID,Tissue,TRV,Cell_count,Age,Sex,TCRgd_freq,Sort_day,
                                             Cell_count_category,Age_category),function(x){return(divsample(x,50,100,'shannon'))})

## cumulative curve shannon diversity (Figure 4B)

divshannon_summary_IEL <- filter(divshannon_summary_minus_uncat, Tissue=="IEL")
divshannon_summary_PBL <- filter(divshannon_summary_minus_uncat, Tissue=="PBL")

g = ggplot(divshannon_summary_minus_uncat, aes(median, colour = Group))
g = g + stat_ecdf()
g = g + theme_classic()
g = g + labs(title = "", x = "Median Shannon Diversity", y = "Cumulative Probability")
g = g + scale_colour_manual(values=c("Control"="grey30", "Active"="firebrick2", "GFD"="blue"))
g = g + ppt
g = g + facet_grid(~Chain)
g

ggsave(paste(plot_folder,'cumulative_dist_IELvsPBL.pdf'),width=12,height=8)

ggsave(paste(plot_folder,'cumulative_dist_TRDvsTRG.pdf'),width=12,height=8)

# sort day effect on diversity 

sortdayeffect <- read.csv("data/divshannon_summary_sortdays.csv")

g = ggplot(sortdayeffect,aes(Sort_day, median))
g = g + geom_boxplot()
g = g + geom_jitter()
#g = g + facet_grid(Chain~Group)
g
ggsave(paste(plot_folder,'sort day diversity effect.png'),width=10,height=5)

# age effect on diversity 

g = ggplot(divshannon_summary,aes(Age_category, median))
g = g + geom_boxplot()
g = g + geom_jitter()
#g = g + facet_grid(Chain~Group)
g
ggsave(paste(plot_folder,'Age diversity effect.png'),width=10,height=5)

# cell count effect on diversity 

g = ggplot(divshannon_summary,aes(Cell_count_category, median))
g = g + geom_boxplot(outlier.size=0)
g = g + geom_jitter()
g
ggsave(paste(plot_folder,'Cell count diversity effect.png'),width=10,height=5)

## Testing for correlations between factors in data set (focussed on shannon diversity here)

data="divshannon_summary"

balls <- divshannon_summary

IEL <- filter(divshannon_summary, Tissue == "IEL")
PBL <- filter(divshannon_summary, Tissue == "PBL")
Active <- filter(divshannon_summary, Group=="Active", Tissue=="IEL")

a=Active$median
b=Active$Cell_count
c=Active$Age
d=Active$TCRgd_freq
e=Active$Sex
f=Active$Sort_day
g=Active$Tissue
h=Active$Group
i=Active$Chain
j=Active$TRV

plot(a,f)

qqnorm()
qqline()
qqplot(a,d)

fit <-lm(a ~ e)
fit <-lm(a ~ h+i+b+c+d+e+f)
summary(fit)
plot(fit)

## Testing for correlations between factors in data set (focussed on shannon diversity of individual gammas here)

data="div_shannon_summary_trv_gammas"

balls <- div_shannon_trv_gammas

IEL <- filter(div_shannon_trv_gammas, Tissue == "IEL")
PBL <- filter(div_shannon_trv_gammas, Tissue == "PBL")
Active <- filter(div_shannon_trv_gammas, Group=="Active", Tissue=="IEL")

a=IEL$median
b=IEL$Cell_count
c=IEL$Age
d=IEL$TCRgd_freq
e=IEL$Sex
f=IEL$Sort_day
g=IEL$Tissue
h=IEL$Group
i=IEL$Chain
j=IEL$TRV

plot(h,a)

qqnorm()
qqline()
qqplot(a,d)

fit <-lm(a ~ e)
fit <-lm(a ~ h+i+b+c+d+e+f)
summary(fit)
plot(fit)

## Testing for correlations between factors in data set (focussed on Vgamma usage here)

data="trv_summary_trg"

balls <- filter(trv_summary_trg, Group=="GFD")

a=balls$freq
b=balls$Cell_count
c=balls$Age
d=balls$TCRgd_freq
e=balls$Sex
f=balls$Sort_day
g=balls$Tissue
h=balls$Group
i=balls$Chain
j=balls$TRV

p = plot(j,g)

qqnorm()
qqline()
qqplot()

fit <-lm(a ~ h)
summary(fit)
plot(fit)

# PBL vs. IEL for Active vs. GFD Median Subsampled Shannon or Simpson # Figure 4b panel I

divshannon_summary_minus_uncat <- filter(divshannon_summary, Group == "Control" | Group == "Active" | Group == "GFD")
divshannon_summary_minus_uncat$Tissue <- factor(divshannon_summary_minus_uncat$Tissue, levels =c("PBL", "IEL"))
divshannon_summary_minus_uncat$Group <- factor(divshannon_summary_minus_uncat$Group, levels =c("Control", "Active", "GFD"))

divshannon_summary_A_G <- filter(divshannon_summary, Group == "Active" | Group == "GFD")
divshannon_summary_A_G$Tissue <- factor(divshannon_summary_A_G$Tissue, levels =c("PBL", "IEL"))

divsimpson_summary_A_G <- filter(divsimpson_summary, Group == "Active" | Group == "GFD")
divsimpson_summary_A_G$Tissue <- factor(divsimpson_summary_A_G$Tissue, levels =c("PBL", "IEL"))

g = ggplot(divshannon_summary_minus_uncat,aes(Tissue,median,fill=Group))
g = g + geom_violin(adjust=.5)
g = g + geom_point(aes(fill=Group), position = position_jitterdodge(dodge.width = 1), size = 1)
g = g + theme_classic()
g = g + scale_fill_manual(values=c('Control'='grey50','Active'='grey70','GFD'='grey95'))
g = g + labs(y = "Median Shannon Diversity", x="", title="Gamma and Delta CDR3 diversity")
g = g + ppt
g = g + theme(title = element_text(vjust=2))
#g = g + facet_wrap(~Group)
g

#A vs G
ggsave(paste(plot_folder,'PBL vs. IEL for Active vs. GFD Median Subsampled Shannon.pdf'),width=12,height=8)

#Everyone by disease
ggsave(paste(plot_folder,'PBL vs. IEL for All Median Subsampled Shannon.png'),width=12,height=8)

stats: 
  IEL <- filter(divshannon_summary_A_G, Tissue=="IEL")
  wilcox.test(median~Group, IEL, alternative = c("two.sided"))
  compute.p.values(median~Group, IEL)

# gamma vs delta for Active vs. GFD Median Subsampled Shannon # Figure 4b panel 2

divshannon_summary_minus_uncat <- filter(divshannon_summary, Group == "Control" | Group == "Active" | Group == "GFD")
divshannon_summary_minus_uncat$Tissue <- factor(divshannon_summary_minus_uncat$Tissue, levels =c("PBL", "IEL"))
divshannon_summary_minus_uncat$Chain <- factor(divshannon_summary_minus_uncat$Chain, levels =c("TRG", "TRD"))

divshannon_summary_A_G <- filter(divshannon_summary, Group == "Active" | Group == "GFD")
divshannon_summary_A_G$Tissue <- factor(divshannon_summary_A_G$Tissue, levels =c("PBL", "IEL"))
divshannon_summary_A_G$Chain <- factor(divshannon_summary_A_G$Chain, levels =c("TRG", "TRD"))

g = ggplot(divshannon_summary_minus_uncat,aes(Tissue,median,fill=Group))
g = g + geom_violin(adjust=.5)
g = g + geom_point(aes(fill=Group), position = position_jitterdodge(dodge.width = 1), size = 1)
g = g + theme_classic()
g = g + scale_fill_manual(values=c('Control'='grey50','Active'='grey70','GFD'='grey95'))
g = g + labs(y = "Median Shannon Diversity", x="", title="CDR3 diversity")
g = g + ppt
g = g + theme(title = element_text(vjust=2))
g = g + facet_wrap(~Chain)
g

#A vs. G
ggsave(paste(plot_folder,'gamma vs delta for Active vs. GFD Median Subsampled Shannon.pdf'),width=12,height=8)

#All
ggsave(paste(plot_folder,'gamma vs delta for All Median Subsampled Shannon.png'),width=12,height=8)

# Looking at TRV diversity at level of individual delta or gamma chains 

div_shannon_trv_gammas <- filter(divshannon_summary_trv, TRV=="2" | TRV=="3" | TRV=="4" | TRV=="5" | TRV=="8" | TRV=="9" | TRV=="10")
div_shannon_trv_gammas$TRV <- factor(div_shannon_trv_gammas$TRV, levels =c("1", "2", "3", "4", "5", "8", "9", "10"))


g = ggplot(div_shannon_trv_gammas,aes(TRV,median))
g = g + geom_boxplot(outlier.size=0)
g = g + geom_point(aes(fill=TRV), position = position_jitterdodge(), size = 1.5)
g = g + theme_classic()
#g = g + scale_fill_manual(values=c('Control'='grey50','Active'='grey70','GFD'='grey95'))
g = g + labs(y = "Median Shannon Diversity", x="", title="CDR3 diversity")
g = g + paper
g = g + theme(title = element_text(vjust=2))
g = g + facet_grid(Group~Tissue)
g

#A vs. G
ggsave(paste(plot_folder,'gamma vs delta for Active vs. GFD Median Subsampled Shannon.pdf'),width=12,height=8)

#All
ggsave(paste(plot_folder,'gamma vs delta for All Median Subsampled Shannon.png'),width=12,height=8)

# subsampled summary plot
g = ggplot(divshannon_summary,aes(factor(ID),median,color=Tissue,shape=Chain))
g = g + geom_point(size=3)
g = g + geom_errorbar(aes(ymin=lq,ymax=uq))
g = g + facet_grid(~Group,scales='free_x',space='free_x')
g = g + xlab('Patient ID')
g = g + ylab('Shannon Diversity Index')
g = g + scale_color_manual(values=c('IEL'='orange','PBL'='mediumorchid2'))
g = g + theme(text = element_text(size=25))
g
ggsave(paste(plot_folder,'shannon_subsample_summary.png',sep=''),width=10,height=5)

## tissue and chain subsampled
shannon_tissue = dcast(divshannon_summary,Group + ID + Chain ~ Tissue,value.var='median')
g = ggplot(shannon_tissue,aes(IEL,PBL))
g = g + geom_point(aes(color=Group,shape=Chain),size=4)
g = g + geom_abline(a=1)
g = g + scale_color_manual(values=c('Control'='deepskyblue','Active'='tomato2','GFD'='lawngreen'))
g = g + theme(text = element_text(size=25))
g
ggsave(paste(plot_folder,'shannon_subsample_tissue.png',sep=''),width=9,height=5)

shannon_chain = dcast(divshannon_summary,Group + ID + Tissue ~ Chain,value.var='median')
g = ggplot(shannon_chain,aes(TRG,TRD))
g = g + geom_point(aes(color=Group,shape=Tissue),size=4)
g = g + geom_abline(a=1)
g = g + scale_color_manual(values=c('Control'='deepskyblue','Active'='tomato2','GFD'='lawngreen'))
g = g + theme(text = element_text(size=25))
g
ggsave(paste(plot_folder,'shannon_subsample_chain.png',sep=''),width=9,height=5)

## simpson subsampled
g = ggplot(divsimpson_summary,aes(factor(ID),median,color=Tissue,shape=Chain))
g = g + geom_point(size=3)
g = g + geom_errorbar(aes(ymin=lq,ymax=uq))
g = g + facet_grid(~Group,scales='free_x',space='free_x')
g = g + xlab('Patient ID')
g = g + ylab('Simpson Diversity Index')
g = g + scale_color_manual(values=c('IEL'='orange','PBL'='mediumorchid2'))
g = g + theme(text = element_text(size=25))
g
ggsave(paste(plot_folder,'simpson_subsample_summary.png',sep=''),width=10,height=5)

## tissue and chain subsampled
simpson_tissue = dcast(divsimpson_summary,Group + ID + Chain ~ Tissue,value.var='median')
g = ggplot(simpson_tissue,aes(IEL,PBL))
g = g + geom_point(aes(color=Group,shape=Chain),size=3)
g = g + geom_abline(a=1)
g = g + scale_color_manual(values=c('Control'='deepskyblue','Active'='tomato2','GFD'='lawngreen'))
g = g + theme(text = element_text(size=25))
g
ggsave(paste(plot_folder,'simpson_subsample_tissue.png',sep=''),width=10,height=5)

simpson_chain = dcast(divsimpson_summary,Group + ID + Tissue ~ Chain,value.var='median')
g = ggplot(simpson_chain,aes(TRG,TRD))
g = g + geom_point(aes(color=Group,shape=Tissue),size=3)
g = g + geom_abline(a=1)
g = g + scale_color_manual(values=c('Control'='deepskyblue','Active'='tomato2','GFD'='lawngreen'))
g = g + theme(text = element_text(size=25))
g
ggsave(paste(plot_folder,'simpson_subsample_chain.png',sep=''),width=10,height=5)


### some CDR3 sequence statistics ###
# define amino acid groupings
# -- simply basic,acidic,polar,nonpolar
hydrophobic_AA = c('A','V','L','I','P','F','M','W','G','C')
hydrophilic_AA = c('N','Q','S','T','Y')
positive_AA = c('K','R','H')
negative_AA = c('D','E')

# function for testing number of group representatives
test_AA = function(seq,AA){
  seq_vector = unlist(strsplit(seq,""))
  return(sum(seq_vector %in% AA)/nchar(seq))
}

test_hscales = function(seq,hscales,scale){
  seq_vector = unlist(strsplit(seq,""))
  return(sum(hscales[[scale]][match(seq_vector,hscales$Residue)])/nchar(seq))
}

## more refined hydrophobicity scales
hscales = read.table('data/hydrophobicity.txt',sep='\t',header=T)

# add cdr3 properties to complete data set
cdr3_properties = all_patients.df
cdr3_properties$CDR3.length = nchar(all_patients.df$CDR3)
cdr3_properties$hydrophobic = unlist(lapply(all_patients.df$CDR3,function(x){return(test_AA(x,hydrophobic_AA))}))
cdr3_properties$hydrophilic = unlist(lapply(all_patients.df$CDR3,function(x){return(test_AA(x,hydrophilic_AA))}))
cdr3_properties$positive = unlist(lapply(all_patients.df$CDR3,function(x){return(test_AA(x,positive_AA))}))
cdr3_properties$negative = unlist(lapply(all_patients.df$CDR3,function(x){return(test_AA(x,negative_AA))}))
## better hydrophobic scales
cdr3_properties$KD = unlist(lapply(all_patients.df$CDR3,function(x){return(test_hscales(x,hscales,'KD'))}))
cdr3_properties$WW = unlist(lapply(all_patients.df$CDR3,function(x){return(test_hscales(x,hscales,'WW'))}))
cdr3_properties$HH = unlist(lapply(all_patients.df$CDR3,function(x){return(test_hscales(x,hscales,'HH'))}))

# expand by count data so boxplot statistics are correct
cdr3_prop_expand = cdr3_properties[rep(1:nrow(cdr3_properties),cdr3_properties[['Count']]),]

## plot properties by group
for (i in 1:length(tissues)){
  # length
  g = ggplot(cdr3_prop_expand[cdr3_prop_expand$Tissue==tissues[i],],aes(factor(ID),CDR3.length))
  g = g + geom_boxplot(aes(color=Chain))
  g = g + facet_grid(~Group,scales='free_x',space='free_x')
  g = g + xlab('Patient ID')
  g = g + ylab('CDR3 Length')
  g = g + scale_color_manual(values=c('TRD'='deeppink','TRG'='slateblue'))
  g = g + theme(text = element_text(size=25))
  g
  ggsave(paste(plot_folder,'cdr3length_',tissues[i],'.png',sep=''),width=10,height=5)
  
   # hydrophobic
  g = ggplot(cdr3_prop_expand[cdr3_prop_expand$Tissue==tissues[i],],aes(factor(ID),hydrophobic))
  g = g + geom_boxplot(aes(color=Chain))
  g = g + facet_grid(~Group,scales='free_x',space='free_x')
  g = g + xlab('Patient ID')
  g = g + scale_color_manual(values=c('TRD'='deeppink','TRG'='slateblue'))
  g = g + theme(text = element_text(size=25))
  g
  ggsave(paste(plot_folder,'hydrophobic_',tissues[i],'.png',sep=''),width=10,height=5)

  # hydrophilic
  g = ggplot(cdr3_prop_expand[cdr3_prop_expand$Tissue==tissues[i],],aes(factor(ID),hydrophilic))
  g = g + geom_boxplot(aes(color=Chain))
  g = g + facet_grid(~Group,scales='free_x',space='free_x')
  g = g + xlab('Patient ID')
  g = g + scale_color_manual(values=c('TRD'='deeppink','TRG'='slateblue'))
  g = g + theme(text = element_text(size=25))
  g
  ggsave(paste(plot_folder,'hydrophilic_',tissues[i],'.png',sep=''),width=10,height=5)
  
  # positive
  g = ggplot(cdr3_prop_expand[cdr3_prop_expand$Tissue==tissues[i],],aes(factor(ID),positive))
  g = g + geom_boxplot(aes(color=Chain))
  g = g + facet_grid(~Group,scales='free_x',space='free_x')
  g = g + xlab('Patient ID')
  g = g + scale_color_manual(values=c('TRD'='deeppink','TRG'='slateblue'))
  g = g + theme(text = element_text(size=25))
  g
  ggsave(paste(plot_folder,'positive_',tissues[i],'.png',sep=''),width=10,height=5)
  
  # negative
  g = ggplot(cdr3_prop_expand[cdr3_prop_expand$Tissue==tissues[i],],aes(factor(ID),negative))
  g = g + geom_boxplot(aes(color=Chain))
  g = g + facet_grid(~Group,scales='free_x',space='free_x')
  g = g + xlab('Patient ID')
  g = g + scale_color_manual(values=c('TRD'='deeppink','TRG'='slateblue'))
  g = g + theme(text = element_text(size=25))
  g
  ggsave(paste(plot_folder,'negative_',tissues[i],'.png',sep=''),width=10,height=5)
  
  # more refined hydrophobic scales
  # Kyte Doolittle (KD)
  g = ggplot(cdr3_prop_expand[cdr3_prop_expand$Tissue==tissues[i],],aes(factor(ID),KD))
  g = g + geom_boxplot(aes(color=Chain))
  g = g + facet_grid(~Group,scales='free_x',space='free_x')
  g = g + xlab('Patient ID')
  g = g + scale_color_manual(values=c('TRD'='deeppink','TRG'='slateblue'))
  g = g + theme(text = element_text(size=25))
  g
  ggsave(paste(plot_folder,'hydrophobic_KD_',tissues[i],'.png',sep=''),width=10,height=5)
  
  # Whimley White (WW)
  g = ggplot(cdr3_prop_expand[cdr3_prop_expand$Tissue==tissues[i],],aes(factor(ID),WW))
  g = g + geom_boxplot(aes(color=Chain))
  g = g + facet_grid(~Group,scales='free_x',space='free_x')
  g = g + xlab('Patient ID')
  g = g + scale_color_manual(values=c('TRD'='deeppink','TRG'='slateblue'))
  g = g + theme(text = element_text(size=25))
  g
  ggsave(paste(plot_folder,'hydrophobic_WW_',tissues[i],'.png',sep=''),width=10,height=5)
  
  # Hessa von Heigne (HH)
  g = ggplot(cdr3_prop_expand[cdr3_prop_expand$Tissue==tissues[i],],aes(factor(ID),HH))
  g = g + geom_boxplot(aes(color=Chain))
  g = g + facet_grid(~Group,scales='free_x',space='free_x')
  g = g + xlab('Patient ID')
  g = g + scale_color_manual(values=c('TRD'='deeppink','TRG'='slateblue'))
  g = g + theme(text = element_text(size=25))
  g
  ggsave(paste(plot_folder,'hydrophobic_HH_',tissues[i],'.png',sep=''),width=10,height=5)
  
  }


#################################
## multiple sequence alignment ##
#################################
# split into gamma and delta chains
all_patients_TRG.df = all_patients.df[all_patients.df$Chain=='TRG',]
all_patients_TRD.df = all_patients.df[all_patients.df$Chain=='TRD',]
all_patients_TRG_IEL.df = all_patients_TRG.df[all_patients_TRG.df$Tissue=='IEL',]
all_patients_TRD_IEL.df = all_patients_TRD.df[all_patients_TRD.df$Tissue=='IEL',]
all_patients_TRG_PBL.df = all_patients_TRG.df[all_patients_TRG.df$Tissue=='PBL',]
all_patients_TRD_PBL.df = all_patients_TRD.df[all_patients_TRD.df$Tissue=='PBL',]

## first try pairwise alignment between CDR3 sequences of all unique clones ##
require(Biostrings)

calc_pairwise = function(patients.df,subMat){
  ncdr3 = nrow(patients.df)
  all_scores = matrix(NA,nrow=ncdr3,ncol=ncdr3)
  for (i in 1:ncdr3){
    print(i)
    for (j in 1:ncdr3){
      alm_score = pairwiseAlignment(patients.df[i,]$CDR3,patients.df[j,]$CDR3,
                              substitutionMatrix=subMat,scoreOnly=T)
      all_scores[i,j] = alm_score
    }
  }
  return(all_scores)
}

## this takes a long time ##
# TRD_scores = calc_pairwise(all_patients_TRD.df,"BLOSUM62")
# TRG_scores = calc_pairwise(all_patients_TRG.df,"BLOSUM62")
TRG_scores_IEL = calc_pairwise(all_patients_TRG_IEL.df,"BLOSUM62")
TRD_scores_IEL = calc_pairwise(all_patients_TRD_IEL.df,"BLOSUM62")
TRG_scores_PBL = calc_pairwise(all_patients_TRG_PBL.df,"BLOSUM62")
TRD_scores_PBL = calc_pairwise(all_patients_TRD_PBL.df,"BLOSUM62")



# save these to file
# save(TRD_scores,TRG_scores,all_patients_TRD.df,all_patients_TRG.df,
#      file='data/pairwise_scores.RData')

save(TRD_scores_IEL,TRG_scores_IEL,
           file='data/pairwise_scores.RData')

# plot distances as heatmap
require(NMF)
quakescale = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), space="Lab")

plot_pairwise = function(scores,patients.df,chain,cwidth){
  ann_colors = list(Disease_State=c('red','blue','lawngreen'),
                    Patient_ID=colorRampPalette(brewer.pal(8,'Set2'))(length(unique(patients.df$ID))),
                    Tissue=c('orange'),
                    V_chain=colorRampPalette(brewer.pal(8,'Set1'))(length(unique(patients.df$TRV))),
                    CDR3_freq=c('darkgreen'))
  
  annotation = data.frame(Disease_State=patients.df$Group,Patient_ID=factor(patients.df$ID),Tissue=patients.df$Tissue,
                          V_Chain=factor(patients.df$TRV),CDR3_freq=patients.df$Freq)
  
  return(aheatmap(scores,color=quakescale(100),annCol=annotation,annColors=ann_colors,annRow=annotation,
           filename=paste('plots/pairwise_',chain,'.png',sep=''),cellwidth=cwidth,cellheight=cwidth,
           labCol=NA,labRow=NA,Colv='Rowv'))
}

TRD_ahm = plot_pairwise(TRD_scores,all_patients_TRD.df,'TRD',1.2)
TRG_ahm = plot_pairwise(TRG_scores,all_patients_TRG.df,'TRG',1.5)
TRD_IEL_ahm = plot_pairwise(TRD_scores_IEL,all_patients_TRD_IEL.df,'TRD',1.2)
TRG_IEL_ahm = plot_pairwise(TRG_scores_IEL,all_patients_TRG_IEL.df,'TRG',1.5)

# select particular heatmap clusters
## generic function to traverse dendrogram and retrieve leaf labels
get_branch_labels = function(dend,bvector){
  # tree base case, bvector empty
  if (length(bvector) == 0) {
    # if we reach a leaf
    if (attributes(dend)$members == 1) {
      return(attributes(dend)$label)
    }
    return(c(get_branch_labels(dend[[1]],bvector),get_branch_labels(dend[[2]],bvector)))
  }
  return(get_branch_labels(dend[[bvector[1]]],bvector[-1]))
}

## first TRD
## returns indices of all_patients_TRD.df
TRD_A_indices = get_branch_labels(TRD_ahm$Colv,c(2,1,1)) 
TRD_B_indices = get_branch_labels(TRD_ahm$Colv,c(2,1,2))
TRD_C_indices = get_branch_labels(TRD_ahm$Colv,c(2,2,1))

# add label to data frame for clusters
all_patients_TRD_clust.df = all_patients_TRD.df
all_patients_TRD_clust.df$Cluster = 'rest'
all_patients_TRD_clust.df$Cluster[TRD_A_indices] = 'A'
all_patients_TRD_clust.df$Cluster[TRD_B_indices] = 'B'
all_patients_TRD_clust.df$Cluster[TRD_C_indices] = 'C'

TRD_clustgroup_summary = ddply(all_patients_TRD_clust.df,.(Group,Cluster),summarize,
                               counts=sum(Count))
g = ggplot(TRD_clustgroup_summary,aes(factor(Cluster),counts,fill=Group))
g = g + geom_bar(stat='identity')
# g = g + scale_fill_manual(values=vcolors,limits=levels(trv$TRV))
# g = g + facet_grid(~Group,scales='free_x',space='free_x')
g

TRD_clustgroup_table = dcast(TRD_clustgroup_summary,Group ~ Cluster, value.var='counts')

## now TRG
## returns indices of all_patients_TRG.df
TRG_A_indices = get_branch_labels(TRG_ahm$Colv,c(2,1)) 
TRG_B_indices = get_branch_labels(TRG_ahm$Colv,c(2,2,1)) 
TRG_C_indices = get_branch_labels(TRG_ahm$Colv,c(2,2,2,1)) 
TRG_D_indices = get_branch_labels(TRG_ahm$Colv,c(2,2,2,2,2)) 

# add label to data frame for clusters
all_patients_TRG_clust.df = all_patients_TRG.df
all_patients_TRG_clust.df$Cluster = 'rest'
all_patients_TRG_clust.df$Cluster[TRG_A_indices] = 'A'
all_patients_TRG_clust.df$Cluster[TRG_B_indices] = 'B'
all_patients_TRG_clust.df$Cluster[TRG_C_indices] = 'C'
all_patients_TRG_clust.df$Cluster[TRG_D_indices] = 'D'

TRG_clustgroup_summary = ddply(all_patients_TRG_clust.df,.(Group,Cluster),summarize,
                               counts=sum(Count))
g = ggplot(TRG_clustgroup_summary,aes(factor(Cluster),counts,fill=Group))
g = g + geom_bar(stat='identity')
# g = g + scale_fill_manual(values=vcolors,limits=levels(trv$TRV))
# g = g + facet_grid(~Group,scales='free_x',space='free_x')
g

TRG_clustgroup_table = dcast(TRG_clustgroup_summary,Group ~ Cluster, value.var='counts')


## write sequences to fasta for playing around
write_group_fasta = function(patients.df,chain){
  for (i in 1:length(patient_groups)){
    subpatients.df = patients.df[patients.df$Group==patient_groups[i],]
    cdr3_names = paste(subpatients.df$Group,subpatients.df$ID,subpatients.df$Tissue,
                       subpatients.df$Chain,subpatients.df$TRV,sep='_')
    cdr3.ss = AAStringSet(subpatients.df$CDR3)
    names(cdr3.ss) = cdr3_names
    writeXStringSet(cdr3.ss,paste('data/fasta/',patient_groups[i],'_',chain,'.fa',sep=''))
  }
}

# groups by chain
write_group_fasta(all_patients_TRG.df,'TRG')
write_group_fasta(all_patients_TRD.df,'TRD')

# groups by chain and tissue
write_group_fasta(all_patients_TRG.df[all_patients_TRG.df$Tissue=='IEL',],'TRG_IEL')
write_group_fasta(all_patients_TRG.df[all_patients_TRG.df$Tissue=='PBL',],'TRG_PBL')
write_group_fasta(all_patients_TRD.df[all_patients_TRD.df$Tissue=='IEL',],'TRD_IEL')
write_group_fasta(all_patients_TRD.df[all_patients_TRD.df$Tissue=='PBL',],'TRD_PBL')

write_fasta = function(patients.df,label){
  cdr3_names = paste(patients.df$Group,patients.df$ID,patients.df$Tissue,
                     patients.df$Chain,patients.df$TRV,sep='_')
  cdr3.ss = AAStringSet(patients.df$CDR3)
  names(cdr3.ss) = cdr3_names
  writeXStringSet(cdr3.ss,paste('data/fasta/',label,'.fa',sep=''))
}

# just by chain and tissue
write_fasta(all_patients_TRG.df[all_patients_TRG.df$Tissue=='IEL',],'TRG_IEL')
write_fasta(all_patients_TRG.df[all_patients_TRG.df$Tissue=='PBL',],'TRG_PBL')
write_fasta(all_patients_TRD.df[all_patients_TRD.df$Tissue=='IEL',],'TRD_IEL')
write_fasta(all_patients_TRD.df[all_patients_TRD.df$Tissue=='PBL',],'TRD_PBL')

## write fasta files based on pairwise clusters
# 3 TRD clusters
write_fasta(all_patients_TRD_clust.df[all_patients_TRD_clust.df$Cluster=='A',],'TRD_clust_A')
write_fasta(all_patients_TRD_clust.df[all_patients_TRD_clust.df$Cluster=='B',],'TRD_clust_B')
write_fasta(all_patients_TRD_clust.df[all_patients_TRD_clust.df$Cluster=='C',],'TRD_clust_C')

# 4 TRG clusters
write_fasta(all_patients_TRG_clust.df[all_patients_TRG_clust.df$Cluster=='A',],'TRG_clust_A')
write_fasta(all_patients_TRG_clust.df[all_patients_TRG_clust.df$Cluster=='B',],'TRG_clust_B')
write_fasta(all_patients_TRG_clust.df[all_patients_TRG_clust.df$Cluster=='C',],'TRG_clust_C')
write_fasta(all_patients_TRG_clust.df[all_patients_TRG_clust.df$Cluster=='D',],'TRG_clust_D')

## perform full multiple sequence alignment ##
require(msa) # going to try this new multiple sequence alignment package






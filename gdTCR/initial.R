
# initial look at gdT compiled data set
# testing data format and upload

setwd('~/Documents/Dinner/gd TCR/')

## folder tree organizes patients and corresponding TCR sequence data
# patient group -> patient -> sequencing tables
# for each patient, we have (IEL,PBL) x (TRD,TRG) tables
# all are from Vd1 sorted gdT

# set plot folder
plot_folder = 'plots_nochallenge/'

# define groups and patient ids
# patient_groups = c('Active','Challenge','Control','GFD')
# patient_ids = list(c(22,35,46,47,81),c(1,72,78),c(7,13,40,53),c(3,4,28,33,43))

# without challenge
patient_groups = c('Control','Active','GFD')
patient_ids = list(c(7,13,40,53),c(22,35,46,47,81),c(3,4,28,33,43))

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

## slightly modify TRJ labeling => make '2 or 1' equivalent to '1 or 2'
all_patients.df[all_patients.df$TRJ=='2 or 1',]$TRJ = '1 or 2'

# set order of patient group
all_patients.df$Group = factor(all_patients.df$Group,levels=c('Control','Active','GFD'))

## finding duplicate chains ##
all_patients.df[duplicated(all_patients.df[,c('Chain','TRV','CDR3','TRJ')]),]

# there are only 4 exact duplicates (matching chain,V,CDR3,J)!
all_patients.df[with(all_patients.df,Chain=='TRD' & TRV=='1' & CDR3=='CALGDQRVPIPWTGGYRHTDKLIF' & TRJ=='1'),]
all_patients.df[with(all_patients.df,Chain=='TRD' & TRV=='1' & CDR3=='CALGEYGRGSWGISHTDKLIF' & TRJ=='1'),]
all_patients.df[with(all_patients.df,Chain=='TRG' & TRV=='2' & CDR3=='CATWDGPNYYKKLF' & TRJ=='2'),]
all_patients.df[with(all_patients.df,Chain=='TRG' & TRV=='9' & CDR3=='CALWEVLYKKLF' & TRJ=='1 or 2'),]

## getting all duplicated CDR3 sequences ##
dup_cdr3_indices = which(duplicated(all_patients.df[,c('Chain','CDR3')])) # 13 clones
dup_cdr3 = list()
for (di in dup_cdr3_indices) {
  dup_cdr3 = c(dup_cdr3,list(all_patients.df[with(all_patients.df,Chain==Chain[di] & CDR3==CDR3[di]),]))
}


require(reshape2)
require(plyr)
require(ggplot2)
require(RColorBrewer)
## sampling summary ##
sample_summary = ddply(all_patients.df,.(Group,ID,Tissue,Chain),summarize,samples=sum(Count))
g = ggplot(sample_summary,aes(factor(ID),samples))
g = g + geom_point(aes(color=Tissue,shape=Chain),size=3)
g = g + facet_grid(~Group,scales='free_x',space='free_x')
g = g + xlab('Patient ID')
g
ggsave(paste(plot_folder,'sample_sizes.png',sep=''),width=5,height=3)

## plotting chain distributions by patient ##
trv_summary = ddply(all_patients.df,.(Group,ID,Tissue,Chain,TRV),summarize,
                   count=sum(Count),freq=sum(Freq))

# create barplot for each combination of tissue and chain
tissues = c('IEL','PBL')
chains = c('TRG','TRD')

# two separate data frames by chain
trv_summary_trd = trv_summary[trv_summary$Chain=='TRD',]
trv_summary_trd$TRV = factor(trv_summary_trd$TRV)
trv_summary_trg = trv_summary[trv_summary$Chain=='TRG',]
trv_summary_trg$TRV = factor(trv_summary_trg$TRV)
trv_summary_split = list(trv_summary_trg,trv_summary_trd)

for (i in 1:length(tissues)) {
  for (j in 1:length(chains)) {
    trv = trv_summary_split[[j]][with(trv_summary_split[[j]],Tissue==tissues[i] & Chain==chains[j]),]
    
    if (j == 1) {
      vcolors = colorRampPalette(brewer.pal(8,'Set2'))(length(levels(trv$TRV)))
    }
    else {
      vcolors = colorRampPalette(brewer.pal(8,'Set1'))(length(levels(trv$TRV)))
    }
    # grouped barplot showing proportions
    g = ggplot(trv,aes(factor(ID),freq,fill=TRV))
    g = g + geom_bar(stat='identity')
    g = g + scale_fill_manual(values=vcolors,limits=levels(trv$TRV))
    g = g + facet_grid(~Group,scales='free_x',space='free_x')
    g = g + xlab('Patient ID')
    g

    ggsave(paste(plot_folder,'TRV_freq_',tissues[i],'_',chains[j],'.png',sep=''),width=5,height=3)
  }
}

## barplot all clones ##
clone_summary = all_patients.df[sample(dim(all_patients.df)[1]),]
clone_summary$Clone = rownames(clone_summary)

# create barplot for each combination of tissue and chain
tissues = c('IEL','PBL')
chains = c('TRG','TRD')

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
    g
    
    ggsave(paste(plot_folder,'clone_freq_',tissues[i],'_',chains[j],'.png',sep=''),width=5,height=3)
  }
}

## calculate shannon diversity from cdr3 frequencies ##
require(vegan)

diversity_summary = ddply(all_patients.df,.(Group,ID,Tissue,Chain),summarize,
                    shannon=diversity(Count,index='shannon'),simpson=diversity(Count,index='simpson'))

# scatter comparison by patient and group
## using shannon
g = ggplot(diversity_summary,aes(factor(ID),shannon))
g = g + geom_point(aes(color=Tissue,shape=Chain),size=3)
g = g + facet_grid(~Group,scales='free_x',space='free_x')
g = g + xlab('Patient ID')
g
ggsave(paste(plot_folder,'shannon_summary.png'),width=5,height=3)

shannon_tissue = dcast(diversity_summary,Group + ID + Chain ~ Tissue,value.var='shannon')
g = ggplot(shannon_tissue,aes(IEL,PBL))
g = g + geom_point(aes(color=Group,shape=Chain),size=3)
g = g + geom_abline(a=1)
g = g + scale_color_manual(values=c('Control'='blue','Active'='red','GFD'='green'))
g
ggsave(paste('plots_nochallenge/shannon_tissue.png',sep=''))

shannon_chain = dcast(diversity_summary,Group + ID + Tissue ~ Chain,value.var='shannon')
g = ggplot(shannon_chain,aes(TRG,TRD))
g = g + geom_point(aes(color=Group,shape=Tissue),size=3)
g = g + geom_abline(a=1)
g
ggsave(paste(plot_folder,'shannon_chain.png',sep=''))

## using simpson
g = ggplot(diversity_summary,aes(factor(ID),simpson))
g = g + geom_point(aes(color=Tissue,shape=Chain),size=3)
g = g + facet_grid(~Group,scales='free_x',space='free_x')
g = g + xlab('Patient ID')
g
ggsave(paste(plot_folder,'simpson_summary.png',sep=''),width=5,height=3)

simpson_tissue = dcast(diversity_summary,Group + ID + Chain ~ Tissue,value.var='simpson')
g = ggplot(simpson_tissue,aes(IEL,PBL))
g = g + geom_point(aes(color=Group,shape=Chain),size=3)
g = g + geom_abline(a=1)
g
ggsave(paste(plot_folder,'simpson_tissue.png',sep=''))

simpson_chain = dcast(diversity_summary,Group + ID + Tissue ~ Chain,value.var='simpson')
g = ggplot(simpson_chain,aes(TRG,TRD))
g = g + geom_point(aes(color=Group,shape=Tissue),size=3)
g = g + geom_abline(a=1)
g
ggsave(paste(plot_folder,'simpson_chain.png',sep=''))

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
      subcounts.df = ddply(subcounts_expand.df,.(Group,ID,Tissue,Chain,TRV,CDR3,TRJ),
                           summarize,Count=length(Count))
      divresults[i] = diversity(subcounts.df$Count,index=divindex)
    }
  }
  return(data.frame(median=median(divresults),lq=quantile(divresults,0.25),
                    uq=quantile(divresults,0.75)))
}

## summarize diversity resampling ##
divshannon_summary = ddply(all_patients.df,.(Group,ID,Tissue,Chain),
                           function(x){return(divsample(x,50,100,'shannon'))})
divsimpson_summary = ddply(all_patients.df,.(Group,ID,Tissue,Chain),
                           function(x){return(divsample(x,50,100,'simpson'))})

# subsampled summary plot
g = ggplot(divshannon_summary,aes(factor(ID),median,color=Tissue,shape=Chain))
g = g + geom_point(size=3)
g = g + geom_errorbar(aes(ymin=lq,ymax=uq))
g = g + facet_grid(~Group,scales='free_x',space='free_x')
g = g + xlab('Patient ID')
g
ggsave(paste(plot_folder,'shannon_subsample_summary.png',sep=''),width=5,height=3)

## tissue and chain subsampled
shannon_tissue = dcast(divshannon_summary,Group + ID + Chain ~ Tissue,value.var='median')
g = ggplot(shannon_tissue,aes(IEL,PBL))
g = g + geom_point(aes(color=Group,shape=Chain),size=3)
g = g + geom_abline(a=1)
g
ggsave(paste(plot_folder,'shannon_subsample_tissue.png',sep=''))

shannon_chain = dcast(divshannon_summary,Group + ID + Tissue ~ Chain,value.var='median')
g = ggplot(shannon_chain,aes(TRG,TRD))
g = g + geom_point(aes(color=Group,shape=Tissue),size=3)
g = g + geom_abline(a=1)
g
ggsave(paste(plot_folder,'shannon_subsample_chain.png',sep=''))

## simpson subsampled
g = ggplot(divsimpson_summary,aes(factor(ID),median,color=Tissue,shape=Chain))
g = g + geom_point(size=3)
g = g + geom_errorbar(aes(ymin=lq,ymax=uq))
g = g + facet_grid(~Group,scales='free_x',space='free_x')
g = g + xlab('Patient ID')
g
ggsave(paste(plot_folder,'simpson_subsample_summary.png',sep=''),width=5,height=3)

## tissue and chain subsampled
simpson_tissue = dcast(divsimpson_summary,Group + ID + Chain ~ Tissue,value.var='median')
g = ggplot(simpson_tissue,aes(IEL,PBL))
g = g + geom_point(aes(color=Group,shape=Chain),size=3)
g = g + geom_abline(a=1)
g
ggsave(paste(plot_folder,'simpson_subsample_tissue.png',sep=''))

simpson_chain = dcast(divsimpson_summary,Group + ID + Tissue ~ Chain,value.var='median')
g = ggplot(simpson_chain,aes(TRG,TRD))
g = g + geom_point(aes(color=Group,shape=Tissue),size=3)
g = g + geom_abline(a=1)
g
ggsave(paste(plot_folder,'simpson_subsample_chain.png',sep=''))


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
  g
  ggsave(paste(plot_folder,'cdr3length_',tissues[i],'.png',sep=''),width=5,height=3)
  
  # hydrophobic
  g = ggplot(cdr3_prop_expand[cdr3_prop_expand$Tissue==tissues[i],],aes(factor(ID),hydrophobic))
  g = g + geom_boxplot(aes(color=Chain))
  g = g + facet_grid(~Group,scales='free_x',space='free_x')
  g = g + xlab('Patient ID')
  g
  ggsave(paste(plot_folder,'hydrophobic_',tissues[i],'.png',sep=''),width=5,height=3)

  # hydrophilic
  g = ggplot(cdr3_prop_expand[cdr3_prop_expand$Tissue==tissues[i],],aes(factor(ID),hydrophilic))
  g = g + geom_boxplot(aes(color=Chain))
  g = g + facet_grid(~Group,scales='free_x',space='free_x')
  g = g + xlab('Patient ID')
  g
  ggsave(paste(plot_folder,'hydrophilic_',tissues[i],'.png',sep=''),width=5,height=3)
  
  # positive
  g = ggplot(cdr3_prop_expand[cdr3_prop_expand$Tissue==tissues[i],],aes(factor(ID),positive))
  g = g + geom_boxplot(aes(color=Chain))
  g = g + facet_grid(~Group,scales='free_x',space='free_x')
  g = g = xlab('Patient ID')
  g
  ggsave(paste(plot_folder,'positive_',tissues[i],'.png',sep=''),width=5,height=3)
  
  # negative
  g = ggplot(cdr3_prop_expand[cdr3_prop_expand$Tissue==tissues[i],],aes(factor(ID),negative))
  g = g + geom_boxplot(aes(color=Chain))
  g = g + facet_grid(~Group,scales='free_x',space='free_x')
  g = g + xlab('Patient ID')
  g
  ggsave(paste(plot_folder,'negative_',tissues[i],'.png',sep=''),width=5,height=3)
  
  # more refined hydrophobic scales
  # Kyte Doolittle (KD)
  g = ggplot(cdr3_prop_expand[cdr3_prop_expand$Tissue==tissues[i],],aes(factor(ID),KD))
  g = g + geom_boxplot(aes(color=Chain))
  g = g + facet_grid(~Group,scales='free_x',space='free_x')
  g = g + xlab('Patient ID')
  g
  ggsave(paste(plot_folder,'hydrophobic_KD_',tissues[i],'.png',sep=''),width=5,height=3)
  
  # Whimley White (WW)
  g = ggplot(cdr3_prop_expand[cdr3_prop_expand$Tissue==tissues[i],],aes(factor(ID),WW))
  g = g + geom_boxplot(aes(color=Chain))
  g = g + facet_grid(~Group,scales='free_x',space='free_x')
  g = g + xlab('Patient ID')
  g
  ggsave(paste(plot_folder,'hydrophobic_WW_',tissues[i],'.png',sep=''),width=5,height=3)
  
  # Hessa von Heigne (HH)
  g = ggplot(cdr3_prop_expand[cdr3_prop_expand$Tissue==tissues[i],],aes(factor(ID),HH))
  g = g + geom_boxplot(aes(color=Chain))
  g = g + facet_grid(~Group,scales='free_x',space='free_x')
  g = g + xlab('Patient ID')
  g
  ggsave(paste(plot_folder,'hydrophobic_HH_',tissues[i],'.png',sep=''),width=5,height=3)
}

#################################
## multiple sequence alignment ##
#################################
# split into gamma and delta chains
all_patients_TRG.df = all_patients.df[all_patients.df$Chain=='TRG',]
all_patients_TRD.df = all_patients.df[all_patients.df$Chain=='TRD',]

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

# save these to file
# save(TRD_scores,TRG_scores,all_patients_TRD.df,all_patients_TRG.df,
#      file='data/pairwise_scores.RData')

# plot distances as heatmap
require(NMF)
quakescale = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), space="Lab")

plot_pairwise = function(scores,patients.df,chain,cwidth){
  ann_colors = list(group=c('red','orange','black','blue'),
                    ID=colorRampPalette(brewer.pal(8,'Set2'))(length(unique(patients.df$ID))),
                    tissue=c('grey','black'),
                    V=colorRampPalette(brewer.pal(8,'Set1'))(length(unique(patients.df$TRV))),
                    freq=c('darkgreen'))
  
  annotation = data.frame(group=patients.df$Group,ID=factor(patients.df$ID),tissue=patients.df$Tissue,
                          V=factor(patients.df$TRV),freq=patients.df$Freq)
  
  return(aheatmap(scores,color=quakescale(100),annCol=annotation,annColors=ann_colors,annRow=annotation,
           filename=paste('plots/pairwise_',chain,'.png',sep=''),cellwidth=cwidth,cellheight=cwidth,
           labCol=NA,labRow=NA,Colv='Rowv'))
}

TRD_ahm = plot_pairwise(TRD_scores,all_patients_TRD.df,'TRD',1.2)
TRG_ahm = plot_pairwise(TRG_scores,all_patients_TRG.df,'TRG',1.5)

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





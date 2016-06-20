
# parsing full vquest output -- excel file with multiple sheets
setwd('~/Documents/Dinner/gd TCR/')

# require(xlsx)
require(XLConnect) # this package lets you import all sheets at once
require(Biostrings) # for nt translation

# helper function to extract matches from regex
get_matches = function(mstring,mresult){
  matches = rep(NA,length(mresult))
  for (i in 1:length(mresult)){
    matches[i] = substr(mstring,mresult[i],mresult[i]+attributes(mresult)$match.length[i]-1)
  }
  return(matches)
}

# helper to parse recorded genes and alleles
get_genes = function(gene_strings){
  match_results = gregexpr("TR\\S*",gene_strings)
  
  genes.df = data.frame(Gene=character(0),Allele=character(0),stringsAsFactors=F)
  alts = list()
  for (i in 1:length(gene_strings)){
    matches = get_matches(gene_strings[i],match_results[[i]])
    matches.lst = strsplit(matches,'\\*')
    genes.df[i,] = matches.lst[[1]]
    if (length(matches) > 1)
      alts[i] = matches.lst[-1]
    else
      alts[i] = NA
  }
  
  return(list(genes.df,alts))
}

## try reading one of the vquest spreadsheets
parse_vquest_gamma = function(vfile) {
  # vfile = 'data/vquest/KL304_vquest 2'
  wb = loadWorkbook(vfile)
  tst.lst = readWorksheet(wb,getSheets(wb))
  
  # this is a list of data frames named by the sheet
  summary.df = tst.lst[[which(names(tst.lst)=='Summary')]]
  
  # unique(summary.df$Functionality)
  # always 'productive','No results', 'unproductive (see comment)', or 'unknown (see comment)'
  # only going to retain productive rearrangements
  prod_filter = summary.df$Functionality == 'productive'

  ## extracting V and J gene
  v_genes = paste(summary.df[prod_filter,]$V.GENE.and.allele,
                  summary.df[prod_filter,]$V.REGION.potential.ins.del,sep=' ')
  j_genes = paste(summary.df[prod_filter,]$J.GENE.and.allele,
                  summary.df[prod_filter,]$J.GENE.and.allele.comment,sep=' ')
  
  # data frames with genes and alleles for each productive sample
  v_genes.df = get_genes(v_genes)[[1]]
  colnames(v_genes.df) = c('VGene','VAllele')
  v_alts = get_genes(v_genes)[[2]]
  j_genes.df = get_genes(j_genes)[[1]]
  colnames(j_genes.df) = c('JGene','JAllele')
  j_alts = get_genes(j_genes)[[2]]
  
  ### individual nucleotide sequences component separated
  nt.df = tst.lst[[which(names(tst.lst)=='Nt-sequences')]][prod_filter,]
  # aa.df = tst.lst[[which(names(tst.lst)=='AA-sequences')]][prod_filter,]

  ## extracting nucleotide sequence components
  # the junction is already separated - don't need to perform end matching
  
  nt_concise.df = nt.df[,c('Sequence.number','Sequence.ID','JUNCTION',
                           'X3.V.REGION','P3.V','N.REGION','P5.J','X5.J.REGION')]
  # replace NA values with empty strings
  nt_concise.df[is.na(nt_concise.df)] = ''
  
  ## need to reshuffle these nt to be codon units
  ## only retaining whole codons in conserved V and J ends
  
  # V end without trailing nucleotides
  nt_concise.df$VEND = sapply(nt_concise.df$X3.V.REGION,
                              function(x) subseq(x,end=-((nchar(x) %% 3)+1)))
  # trailing nucleotides from V end
  nt_concise.df$VNT = sapply(nt_concise.df$X3.V.REGION,
                             function(x) subseq(x,width=nchar(x) %% 3,end=-1))
  # V check
  if (!all(paste(nt_concise.df$VEND,nt_concise.df$VNT,sep='') == nt_concise.df$X3.V.REGION))
    print('V ends not correctly parsed')
  # trailing nucleotides from J end
  nt_concise.df$JNT = sapply(nt_concise.df$X5.J.REGION,
                             function(x) subseq(x,width=nchar(x) %% 3,start=1))
  # J end without trailing nucleotides
  nt_concise.df$JEND = sapply(nt_concise.df$X5.J.REGION,
                             function(x) subseq(x,start=(nchar(x) %% 3)+1))
  # J check
  if (!all(paste(nt_concise.df$JNT,nt_concise.df$JEND,sep='') == nt_concise.df$X5.J.REGION))
    print('J ends not correctly parsed')
  
  # now we just need to concatenate the middle regions
  nt_concise.df$NONGERM = apply(nt_concise.df[,c('VNT','P3.V','N.REGION','P5.J','JNT')],1,
                                  function(x) paste(x,collapse=''))
  # checks
  if (!all(paste(nt_concise.df$VEND,nt_concise.df$NONGERM,nt_concise.df$JEND,sep='') == 
      nt_concise.df$JUNCTION))
      print('Components do not match junction')
  
  # convert nt sequences to AA sequences
  nt_concise.df$JUNCTION_AA = sapply(nt_concise.df$JUNCTION,
                                     function(x) toString(translate(DNAString(x))))
  nt_concise.df$VEND_AA = sapply(nt_concise.df$VEND,
                                 function(x) toString(translate(DNAString(x))))
  nt_concise.df$NONGERM_AA = sapply(nt_concise.df$NONGERM,
                                 function(x) toString(translate(DNAString(x))))
  nt_concise.df$JEND_AA = sapply(nt_concise.df$JEND,
                                 function(x) toString(translate(DNAString(x))))
  if (!all(paste(nt_concise.df$VEND_AA,nt_concise.df$NONGERM_AA,nt_concise.df$JEND_AA,sep='') == 
           nt_concise.df$JUNCTION_AA))
    print('AA Components do not match AA junction')

  parsed_tcr.df = data.frame(TRV=v_genes.df$VGene,TRV_Allele=v_genes.df$VAllele,
                             TRJ=j_genes.df$JGene,TRJ_Allele=j_genes.df$JAllele,
                             CDR3=nt_concise.df$JUNCTION_AA,
                             CDR3_VEND=nt_concise.df$VEND_AA,
                             CDR3_NONGERM=nt_concise.df$NONGERM_AA,
                             CDR3_JEND=nt_concise.df$JEND_AA)
  
  return(parsed_tcr.df)
}

## example single usage
parsed_example.df = parse_vquest_gamma('data/vquest/KL304_vquest 2')

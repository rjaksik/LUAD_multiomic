
########### CNV ##############

CNV_samp= read.table('CNV/gdc_sample_sheet.CNV.2021-12-07.tsv',header=T,sep="\t")
IDs = strsplit2(CNV_samp[,6],', ')[,1]
CNV_samp$IDs = IDs
CNV_samp_nodup = CNV_samp[!duplicated(CNV_samp$IDs),]

CNV_table = data.frame()
for (i in 1:nrow(CNV_samp_nodup)) {
  tab = read.table(paste0('CNV/',CNV_samp_nodup[i,'File.ID'],'/',CNV_samp_nodup[i,'File.Name']),header=T)
  tab$GDC_Aliquot = CNV_samp_nodup$IDs[i]
  CNV_table = rbind(CNV_table,tab)
}
colnames(CNV_table)[1] = 'SampleID'
saveRDS(CNV_table,paste0(DatasetName,'_CNV_v1.RDS'))


########### miRNA ##############
 
miRNA_samp = read.table('miRNA/gdc_sample_sheet.miRNA.2021-12-07.tsv',header=T,sep="\t",stringsAsFactors = F)
miRNA_samp = miRNA_samp[miRNA_samp$Sample.Type=="Primary Tumor",]
miRNA_samp_nodup = miRNA_samp[!duplicated(miRNA_samp$Case.ID),]

miRNA_table = data.frame()
for (i in 1:nrow(miRNA_samp_nodup)) {
  tab = read.table(paste0('miRNA/',miRNA_samp_nodup[i,'File.ID'],'/',miRNA_samp_nodup[i,'File.Name']),header=T)
  colnames(tab)[2] = miRNA_samp_nodup[i,'Case.ID']
  if (i==1) {
    miRNA_table = tab[,1:2]
  } else {
    miRNA_table = merge(miRNA_table,tab[,1:2],by="miRNA_ID")  
  }
}
saveRDS(miRNA_table,paste0(DatasetName,'_miRNA_v1.RDS'))




########### mRNA ##############

mRNA_samp = read.table('mRNA/gdc_sample_sheet.mRNA.2021-12-07.tsv',header=T,sep="\t",stringsAsFactors = F)
mRNA_samp = mRNA_samp[mRNA_samp$Sample.Type=="Primary Tumor",]
mRNA_samp_nodup = mRNA_samp[!duplicated(mRNA_samp$Case.ID),]

mRNA_table = data.frame()
for (i in 1:nrow(mRNA_samp_nodup)) {
  tab = read.table(paste0('mRNA/',mRNA_samp_nodup[i,'File.ID'],'/',mRNA_samp_nodup[i,'File.Name']),header=F)
  colnames(tab) = c('ENSEMBL_ID',mRNA_samp_nodup[i,'Case.ID'])
  if (i==1) {
    mRNA_table = tab
  } else {
    mRNA_table = merge(mRNA_table,tab,by="ENSEMBL_ID")  
  }
}
mRNA_table_flt = mRNA_table[grepl('ENSG',mRNA_table[,1]),]
saveRDS(mRNA_table_flt,paste0(DatasetName,'_mRNA_v1.RDS'))




########### MUT ##############

library('VariantAnnotation')
library('plyr')

MUT_samp = read.table('MUT/gdc_sample_sheet.MUT.2021-12-07.tsv',header=T,sep="\t",stringsAsFactors = F)
IDs = strsplit2(MUT_samp[,6],', ')[,1]
MUT_samp$IDs = IDs
MUT_samp_nodup = MUT_samp[!duplicated(MUT_samp$IDs),]

MUT_table = data.frame()
for (i in 1:nrow(MUT_samp_nodup)) {
  vcf = readVcf(paste0('MUT/',MUT_samp_nodup[i,'File.ID'],'/',MUT_samp_nodup[i,'File.Name']))
  tvarTable = data.frame(rowRanges(vcf))[,c('seqnames','start','end')]
  Annot = ldply(info(vcf)$CSQ,function(x) strsplit2(x[1],'\\|'))
  tvarTableAnnot = cbind(tvarTable,Annot)
  tvarTableAnnot$ID = rownames(info(vcf))
  tvarTableAnnot$SampleID = MUT_samp$IDs[i]
  stvarTableAnnot = tvarTableAnnot[tvarTableAnnot$`3` %in% c('HIGH','MODERATE'),c(71,70,1,2,3,4,26,7,8)]
  colnames(stvarTableAnnot) = c('SampleID','VariantID','Chr','Start','End','Allele','VariantType','GeneSymbol','ENSEMBL_ID')
  if (i==1) {
    MUT_table = stvarTableAnnot
  } else {
    MUT_table = rbind(MUT_table,stvarTableAnnot)  
  }
}
saveRDS(MUT_table,paste0(DatasetName,'_MUT_v1.RDS'))





########### Methylation ##############


Meth_samp = read.table('Meth/gdc_sample_sheet.Meth.2021-12-07.tsv',header=T,sep="\t",stringsAsFactors = F)
Meth_samp = Meth_samp[grepl('HumanMethylation450',Meth_samp$File.Name),]
Meth_samp_nodup = Meth_samp[!duplicated(Meth_samp$Case.ID),]

Meth_table = data.frame()
for (i in 1:nrow(Meth_samp_nodup)) {
  tab = read.table(paste0('Meth/',Meth_samp_nodup[i,'File.ID'],'/',Meth_samp_nodup[i,'File.Name']),header=T,sep="\t",quote="")
  colnames(tab)[1:2] = c('CPG_ID',Meth_samp_nodup[i,'Case.ID'])
  if (i==1) {
    Meth_table = tab[,c('CPG_ID','Chromosome','Start','End')]
    geneIDs = strsplit2(tab$Gene_Symbol,';')
    Meth_table$GeneSymbol = geneIDs[,1]
  } 
  Meth_table = merge(Meth_table,tab[,1:2],by="CPG_ID")  
}
saveRDS(Meth_table,paste0(DatasetName,'_Meth_v1.RDS'))


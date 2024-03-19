
### add gene symbols to the TCGA RNA-seq results
mRNAdata = mRNA_table_sel[mRNA_table_sel$pvalue<0.05,samplesToUse]
mRNAdata$ensembl_gene_id = gsub('mRNA:','',strsplit2(rownames(mRNAdata),'.',fixed=T)[,1])
mRNAdata = mRNAdata[mRNAdata$ensembl_gene_id!="NA" & !is.na(mRNAdata$ensembl_gene_id),]

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters = "ensembl_gene_id", 
                attributes = c("ensembl_gene_id", "hgnc_symbol"),
                values = mRNAdata$ensembl_gene_id, mart = mart)

mRNAdata_annot = merge(mRNAdata,G_list,by='ensembl_gene_id',all.x=T)
mRNAdata_annot$hgnc_symbol[mRNAdata_annot$hgnc_symbol=="NA" | mRNAdata_annot$hgnc_symbol==""] = NA
mRNAdata_annot$hgnc_symbol[is.na(mRNAdata_annot$hgnc_symbol)] = mRNAdata_annot$ensembl_gene_id[is.na(mRNAdata_annot$hgnc_symbol)]
rownames(mRNAdata_annot) = paste0('mRNA:',mRNAdata_annot$hgnc_symbol)

Datasets = rbind(Meth_table_sel[,samplesToUse],
                 mRNAdata_annot[,samplesToUse],
                 miRNA_table_sel[,samplesToUse],
                 CNVraw_table_sel[,samplesToUse],
                 CNV_table_sel[,samplesToUse],
                 GeneCNV_table_sel[,samplesToUse],
                 GeneCNVmsig_table_sel[,samplesToUse],
                 CNVpca_table_sel[,samplesToUse],
                 #CNVpcaAggr_table_sel[,samplesToUse],
                 MutInd_table_sel[,samplesToUse],
                 MutGene_table_sel[,samplesToUse],
                 MutMsig_table_sel[,samplesToUse],
                 third_level[,samplesToUse],
                 last_level[,samplesToUse],
                 MUTpca_table_sel[,samplesToUse])

InputData = data.frame(ID=colnames(Datasets),t(na.omit(Datasets)),check.names = F)
InputDataAnnot = merge(AnnorSel,InputData,by=c('ID'))


library('Boruta')
trData = InputDataAnnot[,!colnames(InputDataAnnot) %in% c('ID')]
trData$Surv2yr = factor(trData$Surv2yr)
boruta_output <- Boruta(Surv2yr~., data=(trData), doTrace=0, maxRuns=1000) 

saveRDS(boruta_output,paste0('4.Classification and var importance/TCGA_LUAD_Multiomic_boruta_output_',Factor,'_v',Version,'.RDS'))
boruta_output = readRDS(paste0('4.Classification and var importance/TCGA_LUAD_Multiomic_boruta_output_',Factor,'_v',Version,'.RDS'))


### #evaluate scale dependance (spoiler: there is none)
#boruta_sel <- Boruta(Surv2yr~., data=(trData[,c('Surv2yr',boruta_signif)]), doTrace=0, maxRuns=1000) 
#imps_sel = attStats(boruta_sel)
#imps_sel = imps_sel[order(-imps_sel$meanImp), ]  # descending sort
#boruta_scale <- Boruta(Surv2yr~., data=data.frame(cbind(scale(trData[,boruta_signif]),Surv2yr=trData$Surv2yr)), doTrace=0, maxRuns=1000) 
#imps_scale = attStats(boruta_scale)
#imps_scale = imps_scale[order(-imps_scale$meanImp), ]  # descending sort





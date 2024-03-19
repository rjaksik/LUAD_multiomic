library('foreach')
library('doParallel')
library('doSNOW')
library('tcltk')
library('glmnet')
library('caret')

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1,outfile="") #not to overload your computer
registerDoSNOW(cl)

set.seed(1)
nFeat=100


Factors = c('Surv2yr') #Surv2yr


LassoRes = LassoFeat = data.frame()
for (fID in 1:length(Factors)) {
    for (data_type in Platforms) {  
      print(paste0('Analyzing data: ',data_type))
      
      tData = Datasets[[data_type]]
      if (data_type %in% c('mRNA','miRNA')) {
          tData$pval = tData$pvalue
      }
      
      #sort by p-value and get top 1k features
      if(data_type %in% c('Meth','mRNA','CNVgeneMsig','MutMsig','CNVgene')) {
        tData = head(tData[order(tData$pval),],1000)
      } 
      
      tData = tData[order(tData$pval),]
      tData$pval=NULL
      tDataFiltr = t(tData[,samplesToUse]) 
      tDataFiltr[is.na(tDataFiltr)] = 0
      tDataFiltr = tDataFiltr[,colSums(abs(tDataFiltr))>0 & MatrixGenerics::colVars(tDataFiltr)>0]
       
      #remove corelated genes      
      correlationMatrix <- cor(tDataFiltr)
      # find attributes that are highly corrected (ideally >0.75)
      highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.9,names =F)
      if (length(highlyCorrelated)>0) {
        tDataFiltr_NoCor = tDataFiltr[,-highlyCorrelated]  
      } else {
        tDataFiltr_NoCor = tDataFiltr
      }
      tDataFiltr_NoCor_Annot = data.frame(ID=rownames(tDataFiltr_NoCor),tDataFiltr_NoCor,stringsAsFactors = F,check.names = F)
      
      #create the input dataset
      InputData = merge(AnnorSel,tDataFiltr_NoCor_Annot,by="ID")
      TestFeatures = colnames(InputData)
      TestFeatures = TestFeatures[3:length(TestFeatures)]
      Nlev = length(unique(AnnorSel[,2]))

      #prepare data for LASSO
      trData = InputData[,!colnames(InputData) %in% c('ID',Factors[fID])]
      factorData = InputData[,Factors[fID]]
      
      #run the analysis
      cv.lasso <- cv.glmnet(data.matrix(trData), as.factor(factorData), family='binomial', nfolds = 10, 
                            alpha=1, parallel=T, standardize=TRUE, type.measure='auc',trace.it=T)
  
      lsRes = data.frame(Dataset = data_type,Factor=Factors[fID],Nfeat=1:length(cv.lasso$cvm), AUC=cv.lasso$cvm, SD=cv.lasso$cvsd)
      tlsFeat = as.matrix(coef(cv.lasso, s=cv.lasso$lambda.min))
      tlsFeat = tlsFeat[rownames(tlsFeat)!='(Intercept)',]
      
      lsFeat = data.frame(Dataset = data_type,Factor=Factors[fID],FeatureID = names(tlsFeat),Importance=as.numeric(tlsFeat))
      lsFeat = lsFeat[order(-abs(lsFeat$Importance)),]
      LassoRes = rbind(LassoRes,lsRes)
      LassoFeat = rbind(LassoFeat,lsFeat)
    }
}


saveRDS(LassoRes, paste0('4.Classification and var importance/TCGA_LUAD_Multiomic_LassoRes_',Factor,'_v2.RDS'))
saveRDS(LassoFeat,paste0('4.Classification and var importance/TCGA_LUAD_Multiomic_LassoFeat_',Factor,'_v2.RDS'))

CPTACT_mRNA = readRDS('CPTAC3_LUAD_Multiomic_mRNA_v1.RDS')
CPTACT_Meth = readRDS('CPTAC3_LUAD_Multiomic_Meth_v1.RDS')

CPTACT_Clin = read.table('clinical.cart.2023-03-30/clinical.tsv',header=T,sep='\t',quote = "")

CPTACT_Clin$days_to_death = as.numeric(CPTACT_Clin$days_to_death)
CPTACT_Clin$days_to_last_follow_up = as.numeric(CPTACT_Clin$days_to_last_follow_up)
Nyears = 2
CPTACT_Clin$Surv2yr = NA
CPTACT_Clin$Surv2yr[CPTACT_Clin$days_to_death>=Nyears*365 | (is.na(CPTACT_Clin$days_to_death) & CPTACT_Clin$days_to_last_follow_up>=Nyears*365) ] = 'Survived'
CPTACT_Clin$Surv2yr[CPTACT_Clin$days_to_death<Nyears*365 & CPTACT_Clin$cause_of_death=="Cancer Related"] = 'Died'
table(CPTACT_Clin$Surv2yr)
rownames(CPTACT_Clin) = CPTACT_Clin$case_id

CPTACT_samplesToUse = table(c(colnames(CPTACT_Meth),colnames(CPTACT_mRNA),CPTACT_Clin$case_submitter_id))
CPTACT_samplesToUse = names(CPTACT_samplesToUse[CPTACT_samplesToUse==3])

CPTACT_ClinSel = CPTACT_Clin[CPTACT_Clin$case_submitter_id %in% CPTACT_samplesToUse,]
rownames(CPTACT_ClinSel)= CPTACT_ClinSel$case_submitter_id

DatasetName = 'LUAD_Multiomic'

SampAnnot = read.table('TCGA_annotation_data/clinical.tsv',header=T,sep="\t",quote='')
table(SampAnnot$occupation_duration_years)

SampAnnot$days_to_death = as.numeric(SampAnnot$days_to_death)
SampAnnot$days_to_last_follow_up = as.numeric(SampAnnot$days_to_last_follow_up)
Nyears = 2
SampAnnot$Surv2yr = NA
SampAnnot$Surv2yr[SampAnnot$days_to_death>=Nyears*365 | (is.na(SampAnnot$days_to_death) & SampAnnot$days_to_last_follow_up>=Nyears*365) ] = 'Survived'
SampAnnot$Surv2yr[SampAnnot$days_to_death<Nyears*365] = 'Died'
table(SampAnnot$Surv2yr)
Nyears = 5
SampAnnot$Surv5yr = NA
SampAnnot$Surv5yr[SampAnnot$days_to_death>=Nyears*365 | (is.na(SampAnnot$days_to_death) & SampAnnot$days_to_last_follow_up>=Nyears*365) ] = 'Survived'
SampAnnot$Surv5yr[SampAnnot$days_to_death<Nyears*365] = 'Died'
table(SampAnnot$Surv5yr)

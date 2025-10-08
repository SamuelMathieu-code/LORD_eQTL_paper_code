# OS

library(survival)
library(survminer)
library(dplyr)

bim <- read.table("full_cohort_genotypes.bim", header = F)
colnames(bim) <- c("chr", "var_id", "cg_morgans", "pos", "allele_other", "allele_effect")
fam <- read.table("full_cohort_genotypes.fam", header = F)
colnames(fam) <- c("FID", "IID", "FIID", "MIID", "sex", "phe")

bed <- read_bed("full_cohort_genotypes.bed", 
         names_loci = bim$var_id, names_ind = fam$IID)
bed <- t(bed)


surv_dat <- read.csv("input_data/survival_pfs_os_vars_curated.csv", na.strings = "na")
staging_dat <- read.csv("input_data/staging_and_diagnosis_lord_filtered_no_neoadj_nsclc_no_rec_or_met.csv")
genetic_pca_dat <- read.table("lord_pca_results.eigenvec", header = F)

colnames(genetic_pca_dat) <- c("Record.ID", "id_bis", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
genetic_pca_dat <- genetic_pca_dat[,c("Record.ID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]
genetic_pca_dat$Record.ID <- as.factor(genetic_pca_dat$Record.ID)

staging_dat$Record.ID <- as.factor(staging_dat$Record.ID)
staging_dat <- staging_dat[!is.na(staging_dat$Histopathological.specification),]
staging_dat$Histopathological.specification <- as.factor(staging_dat$Histopathological.specification)
staging_dat <- staging_dat[!is.na(staging_dat$Histological.type),]
staging_dat$Histological.type <- as.factor(staging_dat$Histological.type)
staging_dat <- staging_dat[!is.na(staging_dat$Pathologic.stade..9th.edition),]
staging_dat$Pathologic.stade..9th.edition <- as.factor(staging_dat$Pathologic.stade..9th.edition)
staging_dat <- staging_dat[,c("Record.ID", "Histopathological.specification", "Histological.type", "Pathologic.stade..9th.edition")]

surv_dat$Record.ID <- as.factor(surv_dat$Record.ID)
surv_dat$Sex.at.birth <- as.factor(surv_dat$Sex.at.birth)
surv_dat$Smoking.status <- as.factor(surv_dat$Smoking.status)

os_dat <- surv_dat[,c("Record.ID", "Smoking.status", "Vital.status", "O.S...2024.", "Sex.at.birth", "age.year")]
os_dat$event <- (os_dat$Vital.status == "Deceased") & (os_dat$O.S...2024. < 5*365)

mymin <- function(myvec, max_value){
    out <- c()
    for(i in 1:length(myvec)){
        out[i] <- min(myvec[i], max_value)
    }
    return(out)
}

os_dat$time_to_event <- mymin(os_dat$O.S...2024., 5*365)

os_dat <- os_dat[,c("Record.ID", "Smoking.status", "event", "time_to_event", "Sex.at.birth", "age.year")]
dat <- inner_join(os_dat, staging_dat, by = "Record.ID")
dat <- inner_join(dat, genetic_pca_dat, by = "Record.ID")
dat <- dat[complete.cases(dat),]

null_model <- SPACox_Null_Model(Surv(time_to_event, event) ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + Smoking.status + Sex.at.birth + age.year + Pathologic.stade..9th.edition + Histological.type, data = dat,
                                pIDs=dat$Record.ID, gIDs=rownames(bed))

results_os <- SPACox(null_model, bed)

write.csv(results_os, "OS_GWAS_V1_results_LORD.csv")

## PFS

## Load clinical data

surv_dat <- read.csv("input_data/survival_pfs_os_vars_curated.csv", na.strings = "na")
staging_dat <- read.csv("input_data/staging_and_diagnosis_lord_filtered_no_neoadj_nsclc_no_rec_or_met.csv")
genetic_pca_dat <- read.table("lord_pca_results.eigenvec", header = F)

colnames(genetic_pca_dat) <- c("Record.ID", "id_bis", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
genetic_pca_dat <- genetic_pca_dat[,c("Record.ID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]
genetic_pca_dat$Record.ID <- as.factor(genetic_pca_dat$Record.ID)

staging_dat$Record.ID <- as.factor(staging_dat$Record.ID)
staging_dat <- staging_dat[!is.na(staging_dat$Histopathological.specification),]
staging_dat$Histopathological.specification <- as.factor(staging_dat$Histopathological.specification)
staging_dat <- staging_dat[!is.na(staging_dat$Histological.type),]
staging_dat$Histological.type <- as.factor(staging_dat$Histological.type)
staging_dat <- staging_dat[!is.na(staging_dat$Pathologic.stade..9th.edition),]
staging_dat$Pathologic.stade..9th.edition <- as.factor(staging_dat$Pathologic.stade..9th.edition)
staging_dat <- staging_dat[,c("Record.ID", "Histopathological.specification", "Histological.type", "Pathologic.stade..9th.edition")]

surv_dat$Record.ID <- as.factor(surv_dat$Record.ID)
surv_dat$Sex.at.birth <- as.factor(surv_dat$Sex.at.birth)
surv_dat$Smoking.status <- as.factor(surv_dat$Smoking.status)

pfs_dat <- surv_dat[,c("Record.ID", "Smoking.status", "Recurrence.progression", "PFS..2024.", "Sex.at.birth", "age.year")]
pfs_dat$event <- pfs_dat$Recurrence.progression != ""
pfs_dat <- pfs_dat[,c("Record.ID", "Smoking.status", "event", "PFS..2024.", "Sex.at.birth", "age.year")]
dat <- inner_join(pfs_dat, staging_dat, by = "Record.ID")
dat <- inner_join(dat, genetic_pca_dat, by = "Record.ID")
dat <- dat[complete.cases(dat),]

## Null model

model <- coxph(Surv(PFS..2024., event) ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + Smoking.status + Sex.at.birth + age.year + Pathologic.stade..9th.edition + Histological.type + bed[dat$Record.ID,"chr1:11171:CCTTG:C"], data = dat)

null_model <- SPACox_Null_Model(Surv(PFS..2024., event) ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + Smoking.status + Sex.at.birth + age.year + Pathologic.stade..9th.edition + Histological.type, data = dat,
                                pIDs=dat$Record.ID, gIDs=rownames(bed))

## Do the GWAS

results <- SPACox(null_model, bed)

write.csv(results, "PFS_GWAS_V1_results_LORD.csv")

library(edgeR)

# 0) Load raw counts, TPM normalized counts and covariates
# counts <- ...
# covariates <- ...
# tpms <- ...

dge<- DGEList(counts)
keep_tpms <- rowSums(tpms>1.2)>=0.5*1030 # keeping tpm>1.2 in >50% of samples
keep_ensg = rownames(tpms)[keep_tpms]

dge <- dge[rownames(dge) %in% keep_ensg, , keep.lib.sizes=F]

# 2) Prepare covariates
group = factor(ifelse(endsWith(colnames(dge$counts), "Sain"), "lung", "tumor"))
temp <- strsplit(colnames(dge), "_")
patients <- factor(apply(sapply(temp, tail, 2), 2, head, 1))

age = rep(0, length(patients))
smoker = rep(0, length(patients))
perctum = rep(0, length(patients))
sex = rep(0, length(patients))
for (i in 1:length(patients)) {
  patient <- patients[i]
  w <- which(rownames(covariates) == patient)
  age[i] <- covariates$age[w]
  sex[i] <- covariates$sex[w]
  smoker[i] <- covariates$smoker[w]
  perctum[i] <- ifelse(group[i] == "tumor", covariates$Percentage.of.Tumor[w], 0)
}

smoker = as.factor(smoker)
sex = as.factor(sex)
perctum[is.na(perctum)] <- mean(perctum, na.rm = TRUE)

# 3) Normalisation (TMM)
dge <- calcNormFactors(dge)

# 4) fit the limmma-voom model
  # prepare model
mm <- model.matrix(~ 0 + patients + group + perctum)
y <- voom(dge, mm, plot = T)

# fit model
fit <- lmFit(y, mm)

contr <- makeContrasts(grouptumor, levels = colnames(coef(fit)))
results <- contrasts.fit(fit, contr)
results <- eBayes(results)
top.table <- topTable(results, sort.by = "P", n = Inf)
top.table$gene = rownames(top.table)
write.csv(top.table, "limma-voom_paired_results_tmm_1.2_50.csv", row.names = F)
# Cis-eQTL analysis

## 1. Expression normalization

TPM normalized readcounts are inverse-normal transformed separately in each tissue compartment (Lung, Tumor, $\Delta$) with R:

```R
invnorm <- function(x){
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
}
for (gene in row.names(expression_df)) {
  expression_df[gene,] <- invnorm(expression_df[gene,])  
}
```

## 2. Peer Factors

See `peertools.sh` script.

## 3. cis-eQTLs

The `qtltools.sh` script implements a parallelized usage of QTLTools which devides the task in 40 smaller tasks. The outpout from each small task is then writen in separate files. They are then concatenated in the `concat.jl` julia script. Multiple_testing is done with `fdr.jl`.

## 4. Fine-mapping

Fine-mapping was performed with `finemaapping.jl`. Ensure to have a R environement with SusieR installed to connect with Julia using RCall.


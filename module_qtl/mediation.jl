using GLM
import Pkg.build
using DataFrames
using CSV
using MultivariateStats
using SnpArrays
using Pipe: @pipe
include("../ensgid_to_hgnc.jl")
ENV["R_HOME"] = "miniconda3/envs/mediationr/lib/R"
build("RCall")
using RCall
R"library(mediation)"


## define functions

"""
from a dataframe containing columns `mods` and `gene`, get a dict of gene for each module. 
"""
function df_to_dic(d::DataFrame, key::Symbol = :mods, values::Symbol = :gene)
    dic_module_genes = Dict()
    for row in eachrow(d)
        if haskey(dic_module_genes, getproperty(row, key))
            push!(dic_module_genes[getproperty(row, key)], getproperty(row, values))
        else
            dic_module_genes[getproperty(row, key)] = Vector([getproperty(row, values)])
        end
    end
    return dic_module_genes
end

"""
plink_bed_encoding_to_genotype(x::UInt8)

    convert plink bed format encoding in number of alt alleles carried.
"""
function plink_bed_encoding_to_genotype(x::UInt8)::Float16
    if x == 0x00
        return 0.0
    elseif x == 0x01
        return NaN
    else
        return x - 1.0
    end
end

"""
calculate first principal component of expression values of genes
"""
function calceigengene(expression::DataFrame, genes)
    pca_dat = Matrix(expression[:, genes])'
    pc_first = @pipe fit(PCA, pca_dat, maxoutdim = 1) |> predict(_, pca_dat) |> vec
    return pc_first
end

## Load data

covariates_l = CSV.read("covariates_lung_60peer.csv", DataFrame)
covariates_t = CSV.read("covariates_tumor_60peer.csv", DataFrame)
genotypes = SnpData("TOPMed_chr_all_LORD_pairedRNA.reordered")
modules_l = CSV.read("modules_lung.tsv", DataFrame)
modules_t = CSV.read("modules_tumor.tsv", DataFrame)
filter!(:mods => !=(0), modules_l)
filter!(:mods => !=(0), modules_t)

expr_l = CSV.read("lung_tpms.tsv", DataFrame, select = ["id"; modules_l.gene])
expr_t = CSV.read("tumor_tpms.tsv", DataFrame, select = ["id"; modules_t.gene])
expr_l[:,2:end] = log.(expr_l[:,2:end] .+ 0.5)
expr_t[:,2:end] = log.(expr_t[:,2:end] .+ 0.5)

modules_l.mods = "mod" .* string.(modules_l.mods)
modules_t.mods = "mod" .* string.(modules_t.mods)
modules_l = df_to_dic(modules_l, :mods, :gene)
modules_t = df_to_dic(modules_t, :mods, :gene)

## CDK15

var_id = "chr2:200876010:C:T"
genes = modules_l["mod5"]
mediators = ["ENSG00000138395.17"]
covars = covariates_l[:,2:end]
eigengene = calceigengene(expr_l, setdiff(genes, mediators))
var_idx = findfirst(==(var_id), genotypes.snp_info.snpid)
g = plink_bed_encoding_to_genotype.(genotypes.snparray[:, var_idx])
oth, eff = genotypes.snp_info.allele1[var_idx], genotypes.snp_info.allele2[var_idx]

fmed = Term(:y) ~ Term(:g) + sum(Term.(Symbol.(names(covars))))
fout = Term(:eigengene) ~ Term(:y) + Term(:g) + sum(Term.(Symbol.(names(covars))))
@rput eigengene g fmed fout covars

out_med_l = DataFrame()
for m in mediators
    y = expr_l[:,m]
    @rput y

    R"""
    df <- cbind (eigengene = eigengene, y = y, g = g, covars)
    mmed <- lm(fmed, df)
    mout <- lm(fout, df)
    mediation_out <- mediate(mmed, mout, treat = "g", mediator = "y", covariates = colnames(covars), sims = 10000)
    med_eff <- mediation_out$d0
    med_ci_low <- mediation_out$d0.ci[1]
    med_ci_high <- mediation_out$d0.ci[2]
    med_p <- mediation_out$d0.p
    dir_eff <- mediation_out$z0
    dir_ci_low <- mediation_out$z0.ci[1]
    dir_ci_high <- mediation_out$z0.ci[2]
    dir_p <- mediation_out$z0.p
    """
    @rget med_eff med_ci_low med_ci_high med_p dir_eff dir_ci_low dir_ci_high dir_p
    push!(out_med_l, (;var = var_id, mod = "mod5", mediator = m, 
                        med_eff = med_eff, med_ci_low = med_ci_low, med_ci_high = med_ci_high, med_p = med_p,
                        dir_eff = dir_eff, dir_ci_low = dir_ci_low, dir_ci_high = dir_ci_high, dir_p = dir_p))
end

out_med_l

CSV.write("CDK15_mediation.csv", out_med_l)

## SCARA5

var_id = "chr8:27329075:T:C"
genes = modules_l["mod18"]
mediators = ["ENSG00000168079.17"]
covars = covariates_l[:,2:end]
eigengene = calceigengene(expr_l, setdiff(genes, mediators))
var_idx = findfirst(==(var_id), genotypes.snp_info.snpid)
g = plink_bed_encoding_to_genotype.(genotypes.snparray[:, var_idx])
oth, eff = genotypes.snp_info.allele1[var_idx], genotypes.snp_info.allele2[var_idx]

fmed = Term(:y) ~ Term(:g) + sum(Term.(Symbol.(names(covars))))
fout = Term(:eigengene) ~ Term(:y) + Term(:g) + sum(Term.(Symbol.(names(covars))))
@rput eigengene g fmed fout covars

out_med_l = DataFrame()
for m in mediators
    y = expr_l[:,m]
    @rput y

    R"""
    df <- cbind (eigengene = eigengene, y = y, g = g, covars)
    mmed <- lm(fmed, df)
    mout <- lm(fout, df)
    mediation_out <- mediate(mmed, mout, treat = "g", mediator = "y", covariates = colnames(covars), sims = 10000)
    med_eff <- mediation_out$d0
    med_ci_low <- mediation_out$d0.ci[1]
    med_ci_high <- mediation_out$d0.ci[2]
    med_p <- mediation_out$d0.p
    dir_eff <- mediation_out$z0
    dir_ci_low <- mediation_out$z0.ci[1]
    dir_ci_high <- mediation_out$z0.ci[2]
    dir_p <- mediation_out$z0.p
    """
    @rget med_eff med_ci_low med_ci_high med_p dir_eff dir_ci_low dir_ci_high dir_p
    push!(out_med_l, (;var = var_id, mod = "mod5", mediator = m, 
                        med_eff = med_eff, med_ci_low = med_ci_low, med_ci_high = med_ci_high, med_p = med_p,
                        dir_eff = dir_eff, dir_ci_low = dir_ci_low, dir_ci_high = dir_ci_high, dir_p = dir_p))
end

out_med_l

CSV.write("SCARA5_mediation.csv", out_med_l)

using GLM
import Pkg.build
using DataFrames
using CSV
using MultivariateStats
using SnpArrays
using JLD2
using Pipe: @pipe
include("../ensgid_to_hgnc.jl")

# 0 : Define usefull functions
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
    oneregression(covars::DataFrame, genotypes::AbstractVector{UInt8}, var::AbstractString, mod::AbstractString, oth::AbstractString, eff::AbstractString, eigengene::Vector{Float64})

        Make regression of one snp on eigengene.

        # Arguments:
        - `covars`: datafrane of the covariates to include
        - `genotypes`: vector of genotypes in plink .bed format (UInt8)
        - `var`: variant id
        - `mod`: module id
        - `oth`: other allele (encoded by 0x00 homozygote)
        - `eff`: effect allele (encoded by 0x03 homozygote)
        - `eigengene`: Vector containing values of the eigengene.
    """
    function oneregression(covars::DataFrame, genotypes::AbstractVector{UInt8}, var::AbstractString, mod::AbstractString, oth::AbstractString, eff::AbstractString, eigengene::AbstractVector{<:AbstractFloat}; sw::Bool = false)::NamedTuple
        data = DataFrame(covars, copycols = false)
        data.y = eigengene
        data.g = plink_bed_encoding_to_genotype.(genotypes)
        f1 = Term(:y) ~ sum(Term.(Symbol.(names(covars)[2:end]))) + Term(:g)
        m1 = lm(f1, data)
        if sw
            swp = residuals(m1) |> ShapiroWilkTest |> pvalue
            return (var = var, mod = mod, oth = oth, eff = eff, β = coef(m1)[end], se = stderror(m1)[end], pval = coeftable(m1).cols[4][end], tval = coeftable(m1).cols[3][end], shapirowilkp = swp)
        else
            return (var = var, mod = mod, oth = oth, eff = eff, β = coef(m1)[end], se = stderror(m1)[end], pval = coeftable(m1).cols[4][end], tval = coeftable(m1).cols[3][end])
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

# 1: load data : covariates (age, sex, smoke, genPCs, 60 peer), expression tpms, genotypes, variant-mod pairs, modules, ciseqtls
    covariates_l = CSV.read("input_modqtls_v6/covariates_lung_60peer.csv", DataFrame)
    covariates_t = CSV.read("input_modqtls_v6/covariates_tumor_60peer.csv", DataFrame)
    genotypes = SnpData("../qtltools/input_data/TOPMed_chr_all_LORD_pairedRNA.reordered")
    var_mod_pairs_l = CSV.read("reprioritized_modqtls_v6_lung_sig.csv", DataFrame, select = [:var_id, :mod])
    var_mod_pairs_t = CSV.read("reprioritized_modqtls_v6_tumor_sig.csv", DataFrame, select = [:var_id, :mod])
    modules_l = CSV.read("out_wgcna_v6/modules_lung.tsv", DataFrame)
    modules_t = CSV.read("out_wgcna_v6/modules_tumor.tsv", DataFrame)
    filter!(:mods => !=(0), modules_l)
    filter!(:mods => !=(0), modules_t)
    
    expr_l = CSV.read("../lung_tpms.tsv", DataFrame, select = ["id"; modules_l.gene])
    expr_t = CSV.read("../tumor_tpms.tsv", DataFrame, select = ["id"; modules_t.gene])
    expr_l[:,2:end] = log.(expr_l[:,2:end] .+ 0.5)
    expr_t[:,2:end] = log.(expr_t[:,2:end] .+ 0.5)

    modules_l.mods = "mod" .* string.(modules_l.mods)
    modules_t.mods = "mod" .* string.(modules_t.mods)
    modules_l = df_to_dic(modules_l, :mods, :gene)
    modules_t = df_to_dic(modules_t, :mods, :gene)
    ciseqtls_l = CSV.read("../qtltools/Lung_eQTLs_fdr5.csv", DataFrame)
    ciseqtls_t = CSV.read("../qtltools/Tumor_eQTLs_fdr5.csv", DataFrame)
    ciseqtls_l = df_to_dic(ciseqtls_l, :var_id, :phe_id)
    ciseqtls_t = df_to_dic(ciseqtls_t, :var_id, :phe_id)

    # verify patient ids are the same 

    covariates_l.id == covariates_t.id == expr_l.id == expr_t.id == (parse.(Int, genotypes.person_info.iid))
        # ↪ is true

# 2: Re-run regression without cis regulated gene

    # --- PSEUDOCODE ---
    # for each var - mod in file do
    #   find cis genes from var
    #   remove genes from gene list of module
    #   Calculate expression eigen gene for the new list
    #   find var in SnpData
    #   make regression

    # Lung
    # output dataframe columns : variant, module, cis genes (HGNC?) (cleared), other allele, effect allele, beta, se, p

    out_lung = DataFrame()
    for row in eachrow(var_mod_pairs_l)
        genes = copy(modules_l[row.mod])
        cis_genes = intersect(genes, ciseqtls_l[row.var_id])
        filter!(!∈(cis_genes), genes)
        cis_genes = reduce((x, y) -> x * ";" * y, convert_to_hgnc(cis_genes))

        eigengene = calceigengene(expr_l, genes)
        var_idx = findfirst(==(row.var_id), genotypes.snp_info.snpid)
        g = genotypes.snparray[:, var_idx]
        oth, eff = genotypes.snp_info.allele1[var_idx], genotypes.snp_info.allele2[var_idx]
        reg = oneregression(covariates_l, g, row.var_id, row.mod, oth, eff, eigengene)
        push!(out_lung, merge(reg, (:cisgenes => cis_genes,)))
    end

    CSV.write("out_sensitivity_modqtls_v6/no_cisgene_modqtls_lung.csv", out_lung)
    
    
    # Tumor
    out_tumor = DataFrame()
    for row in eachrow(var_mod_pairs_t)
        genes = copy(modules_t[row.mod])
        cis_genes = intersect(genes, ciseqtls_t[row.var_id])
        filter!(!∈(cis_genes), genes)
        cis_genes = String(reduce((x, y) -> x * ";" * y, convert_to_hgnc(cis_genes)))

        eigengene = calceigengene(expr_t, genes)
        var_idx = findfirst(==(row.var_id), genotypes.snp_info.snpid)
        g = genotypes.snparray[:, var_idx]
        oth, eff = genotypes.snp_info.allele1[var_idx], genotypes.snp_info.allele2[var_idx]
        reg = oneregression(covariates_t, g, row.var_id, row.mod, oth, eff, eigengene)
        m = merge(reg, (:cisgenes => cis_genes,))
        push!(out_tumor, m)
    end
    
    CSV.write("out_sensitivity_modqtls_v6/no_cisgene_modqtls_tumor.csv", out_tumor)
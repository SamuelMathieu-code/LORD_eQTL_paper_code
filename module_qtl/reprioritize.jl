using DataFrames
using CSV
import Pipe.@pipe
using MultipleTesting: adjust, BenjaminiHochberg, Bonferroni
include("../utils/ensgid_to_hgnc.jl")

# load data

    gene_types = CSV.read("../utils/genes_ensembl_id_type.csv", DataFrame)

    modqtls_l = CSV.read("out_modqtls_v6/out_lung_60peer_noorderquantnorm.csv", DataFrame)
    modqtls_t = CSV.read("out_modqtls_v6/out_tumor_60peer_noorderquantnorm.csv", DataFrame)
    rename!(modqtls_l, :var => :var_id)
    rename!(modqtls_t, :var => :var_id)

    lung_mods_genes_pairs = CSV.read("out_wgcna_v6/modules_lung.tsv", DataFrame)
    tumor_mods_genes_pairs = CSV.read("out_wgcna_v6/modules_tumor.tsv", DataFrame)
    filter!(:mods => !=(0), lung_mods_genes_pairs)
    filter!(:mods => !=(0), tumor_mods_genes_pairs)

# Use all var-module pairs for which var is a cis-eqtl for a gene in module

    dic_module_genes_l = Dict{String7, Set{String31}}()
    dic_module_genes_t = Dict{String7, Set{String31}}()
    for row in eachrow(lung_mods_genes_pairs)
        if haskey(dic_module_genes_l, "mod$(row.mods)")
            push!(dic_module_genes_l["mod$(row.mods)"], row.gene)
        else
            dic_module_genes_l["mod$(row.mods)"] = Set{String31}([row.gene])
        end
    end
    for row in eachrow(tumor_mods_genes_pairs)
        if haskey(dic_module_genes_t, "mod$(row.mods)")
            push!(dic_module_genes_t["mod$(row.mods)"], row.gene)
        else
            dic_module_genes_t["mod$(row.mods)"] = Set{String31}([row.gene])
        end
    end

    function is_valid_pair(mod::AbstractString, var::AbstractString, dic::Dict{String7, Set{String31}}, cis_eqtls::Set{Tuple{T1, T2}}) where T1 <: AbstractString where T2 <: AbstractString
        valid_pairs = [(g, var) for g in dic[mod]]
        return any([vp ∈ cis_eqtls for vp in valid_pairs])
    end

    ciseqtl_l = CSV.read("../qtltools/Lung_eQTLs_fdr5.csv", DataFrame)
    ciseqtl_l = groupby(ciseqtl_l, :phe_id)
    ciseqtl_l = combine(d -> d[argmin(d.nom_pval),:], ciseqtl_l)
    ciseqtl_l = zip(ciseqtl_l.phe_id, ciseqtl_l.var_id) |> Set

    ciseqtl_t = CSV.read("../qtltools/Tumor_eQTLs_fdr5.csv", Dat)
    ciseqtl_t = groupby(ciseqtl_t, :phe_id)
    ciseqtl_t = combine(d -> d[argmin(d.nom_pval),:], ciseqtl_t)
    ciseqtl_t = zip(ciseqtl_t.phe_id, ciseqtl_t.var_id) |> Set

    reprioritized_modqtls_l = filter([:mod, :var_id] => (mod, var) -> is_valid_pair(mod, var, dic_module_genes_l, ciseqtl_l), modqtls_l)
    reprioritized_modqtls_t = filter([:mod, :var_id] => (mod, var) -> is_valid_pair(mod, var, dic_module_genes_t, ciseqtl_t), modqtls_t)

    CSV.write("reprioritized_modqtls_v6_lung.csv", reprioritized_modqtls_l)
    CSV.write("reprioritized_modqtls_v6_tumor.csv", reprioritized_modqtls_t)

    reprioritized_modqtls_l.padj = adjust(reprioritized_modqtls_l.pval, BenjaminiHochberg())
    reprioritized_modqtls_t.padj = adjust(reprioritized_modqtls_t.pval, BenjaminiHochberg())

    sig_modqtls_l = filter(:padj => <(0.05), reprioritized_modqtls_l)
    sig_modqtls_t = filter(:padj => <(0.05), reprioritized_modqtls_t)

    CSV.write("reprioritized_modqtls_v6_lung_sig.csv", sig_modqtls_l)
    CSV.write("reprioritized_modqtls_v6_tumor_sig.csv", sig_modqtls_t)


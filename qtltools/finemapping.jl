println("Loading packages...")

using JLD2
using Random
using DataFrames
using SnpArrays
using CSV
import Pipe.@pipe
using BenchmarkTools
import Pkg
using DelimitedFiles

ENV["R_HOME"] = "/mnt/raid5_hdd/matsam01/miniconda3/envs/r_susier/lib/R"
Pkg.build("RCall")
using RCall
@rlibrary susieR

println("Loading genotypes...")

genotypes = SnpData("../genotypes/plink_files/TOPMed_chr_all_LORD_pairedRNA.reordered")
genotypes_index = Dict(zip(genotypes.snp_info.snpid, 1:size(genotypes.snp_info, 1)))

println("Loading functions...")


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
    make_cs(var_gene::AbstractDataFrame, y::AbstractVector{<:Real})

make credible sets for each <gene, snp> pair in pairs. pairs must contain a single unique gene value.
    y is the expression vector. The user should verify all :var_id in pairs are present in genotypes.snp_info global variable.
    The user should also verify there no missing values in the genotypes plink files. This would lead to incorrect results.
"""
function make_cs(var_gene::AbstractDataFrame, y::AbstractVector{<:AbstractFloat})
    indexes = Vector{Int32}(undef, size(var_gene, 1))
    pip = Vector{Union{Missing, Float64}}(undef, size(var_gene, 1))
    cs_vec = Vector{Union{Missing, Symbol}}(missing, size(var_gene, 1))

    @inbounds @simd for i in 1:lastindex(indexes)
        snp = var_gene.var_id[i]
        indexes[i] = genotypes_index[snp]
    end
    X = @view genotypes.snparray[:, indexes]
    answer = susie(plink_bed_encoding_to_genotype.(X), y)

    converged = rcopy(answer[:converged])
    if converged
        pip .= rcopy(answer[:pip])
        cs_dic = rcopy(answer[:sets][:cs])
        if !isnothing(cs_dic)
            for (k, v) in cs_dic
                if ndims(v) > 0
                    cs_vec[v] .= k
                else
                    cs_vec[v] = k
                end
            end
        end
    else
        pip .= missing
    end

    out = copy(var_gene)
    out.pip = pip
    out.cs = cs_vec
    out.converged .= converged

    return out
end

println("Loading QTLs and expression Lung...")

frd5_pairs = CSV.read("Lung_eQTLs_fdr5.csv", DataFrame)
fdr5_egenes = Set(unique(frd5_pairs.phe_id))
frd5_pairs = nothing
GC.gc()

pairs = CSV.read("Lung_eQTLs_60_peer.tsv", DataFrame, select = [:phe_id, :var_id])
filter!(:phe_id => ∈(fdr5_egenes), pairs)
pairs = groupby(pairs, :phe_id)

header = readline("input_data/Expression_lung_tpm01_20perc_reads6_20perc.csv")
header = @pipe split(header, ',') |> strip.(_, '\"') |> convert.(String, _) |> getindex(_, 2:lastindex(_))
expression = CSV.read("output_lung_boundtol01_vartol1e3/residuals.csv", DataFrame, header = header, transpose = true)
GC.gc()

println("Using susie Lung...")

all_out = Memory{DataFrame}(undef, length(pairs))

@inbounds for (i, k) in enumerate(keys(pairs))
    table_out = make_cs(pairs[k], expression[!, k.phe_id])
    all_out[i] = table_out
end

println("concatenating and writing outputs Lung...")

CSV.write("fine_maping_lung.csv", reduce(vcat, all_out))

println("loading QTLs and expression Tumor...")

frd5_pairs = CSV.read("Tumor_eQTLs_fdr5.jld2", "dataset")
fdr5_egenes = Set(unique(frd5_pairs.phe_id))
frd5_pairs = nothing
GC.gc()

pairs = CSV.read("Tumor_eQTLs_60_peer.tsv", DataFrame, select = [:phe_id, :var_id])
filter!(:phe_id => ∈(fdr5_egenes), pairs)
pairs = groupby(pairs, :phe_id)

header = readline("input_data/Expression_tumor_tpm01_20perc_reads6_20perc.csv")
header = @pipe split(header, ',') |> strip.(_, '\"') |> convert.(String, _) |> getindex(_, 2:lastindex(_))
expression = CSV.read("output_tumor_boundtol01_vartol1e3_v2/residuals.csv", DataFrame, header = header, transpose = true)
GC.gc()

println("Using susie Tumor...")

all_out = Memory{DataFrame}(undef, length(pairs))

@inbounds for (i, k) in enumerate(keys(pairs))
    table_out = make_cs(pairs[k], expression[!, k.phe_id])
    all_out[i] = table_out
end

println("concatenating and writing outputs Tumor...")

CSV.write("fine_maping_tumor.csv", reduce(vcat, all_out))

println("done")
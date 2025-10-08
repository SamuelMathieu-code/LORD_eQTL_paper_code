using JLD2
using DataFrames
using HypothesisTests
using Distributions
using CSV
using DelimitedFiles
import Pipe.@pipe
using MultipleTesting


# background genes
d_qtls_all = CSV.read("Delta_eQTLs_60_peer.tsv", DataFrame, delim = "\t")
all_tested_genes = unique(getindex.(split.(d_qtls_all.phe_id, "."), 1))
writedlm("all_tested_mapped_genes_delta_qtls.txt", all_tested_genes)

# comparing Delta-eQTLs with a gain of effect in the tumor (absolute value) vs a loss of effect (absolute value)
t_qtls = CSV.read("Tumor_eQTLs_fdr5.csv", DataFrame)
l_qtls = CSV.read("Lung_eQTLs_fdr5.csv", DataFrame)
d_qtls = CSV.read("Delta_eQTLs_60peer_fdr5.csv", DataFrame)
t_qtls.unique_id = t_qtls.phe_id .* "_" .* t_qtls.var_id
l_qtls.unique_id = l_qtls.phe_id .* "_" .* l_qtls.var_id
d_qtls.unique_id = d_qtls.phe_id .* "_" .* d_qtls.var_id


delta_ids = Set(d_qtls.unique_id)

context_qtls = outerjoin(t_qtls, l_qtls, on=:unique_id, renamecols = :_t => :_l)
filter!(:unique_id => ∈(delta_ids), context_qtls)

context_augmented = filter([:slope_t, :slope_l] => (t, l) -> ismissing(l) || (!ismissing(t) && abs(t) > abs(l)), context_qtls)
context_diminished = filter([:slope_t, :slope_l] => (t, l) -> ismissing(t) || (!ismissing(l) && abs(l) > abs(t)), context_qtls)

delta_tumor_mapped = @pipe split.(context_augmented.phe_id_t, ".") |> getindex.(_, 1) |> unique
delta_lung_mapped = @pipe split.(context_diminished.phe_id_l, ".") |> getindex.(_, 1) |> unique

all_tested_genes = readdlm("all_tested_mapped_genes_delta_qtls.txt", String15) |> vec
all_types = [CSV.read("HPA_core_sets_ct_specific_lung/".*file, DataFrame, delim = "\t") 
             for file in readdir("HPA_core_sets_ct_specific_lung")]
cell_type_sets = [Set(data.Ensembl) for data in all_types]

successes_pop = [sum([all_tested_genes[i] ∈ set for i in eachindex(all_tested_genes)]) for set in cell_type_sets]
failures_pop = length(all_tested_genes) .- successes_pop

distributions_lung = [Hypergeometric(s, f, length(delta_lung_mapped)) for (s, f) in zip(successes_pop, failures_pop)]
distributions_tumor = [Hypergeometric(s, f, length(delta_tumor_mapped)) for (s, f) in zip(successes_pop, failures_pop)]

successes_lung = [sum([delta_lung_mapped[i] ∈ set for i in eachindex(delta_lung_mapped)]) for set in cell_type_sets]
successes_tumor = [sum([delta_tumor_mapped[i] ∈ set for i in eachindex(delta_tumor_mapped)]) for set in cell_type_sets]

fold_enrichements_lung = (successes_lung./length(delta_lung_mapped))./(successes_pop./length(all_tested_genes)) 
fold_enrichements_tumor = (successes_tumor./length(delta_tumor_mapped))./(successes_pop./length(all_tested_genes))

pvals_lung = [1-cdf(d, x) for (d, x) in zip(distributions_lung, successes_lung)]
pvals_tumor = [1-cdf(d, x) for (d, x) in zip(distributions_tumor, successes_tumor)]

dataframe_lung = DataFrame(cell_type = readdir("HPA_core_sets_ct_specific_lung"), fe = fold_enrichements_lung, pvalue = pvals_lung)
dataframe_tumor = DataFrame(cell_type = readdir("/HPA_core_sets_ct_specific_lung"), fe = fold_enrichements_tumor, pvalue = pvals_tumor)

dataframe_tumor.padj = adjust(dataframe_tumor.pvalue, BenjaminiHochberg())
dataframe_lung.padj = adjust(dataframe_lung.pvalue, BenjaminiHochberg())

CSV.write("cell_type_enrichement_delta_qtls_tumor_augmented_effects_genes.csv", dataframe_tumor)
CSV.write("cell_type_enrichement_delta_qtls_lung_augmented_effects_genes.csv", dataframe_lung)

using DataFrames
using CSV
using GLM
import StatsBase.mean
import StatsBase.std
using SnpArrays
import Pipe.@pipe
import StatsBase.countmap
include("../utils/plink_to_genotype.jl")
using ProgressMeter
import Logging
Logging.disable_logging(Logging.Warn)
import Base.Threads.@threads
using LinearAlgebra
BLAS.set_num_threads(20)


## load data


expr_tumor = CSV.read("Expression_tumor_tpm01_20perc_reads6_20perc.csv", DataFrame)
rename!(expr_tumor, :Column1 => :id)
sort!(expr_tumor, :id)
covars = CSV.read("covariates_tumor_60peer.csv", DataFrame)
covars.age = (covars.age .- mean(covars.age)) / std(covars.age)
covars.smoking = (covars.smoking .- mean(covars.smoking)) / std(covars.age)
covars.sex = (covars.sex .- mean(covars.sex)) / std(covars.sex)
sort!(covars, :id)

# best eQTL
tumor_qtl = CSV.read("Tumor_eQTLs_60_peer.tsv.gz", DataFrame)
filter!(:best_hit => ==(1), tumor_qtl)
select!(tumor_qtl, [:phe_id, :var_id])
rename!(tumor_qtl, :var_id => :top_var)
tumor_best_var_dict = Dict(zip(tumor_qtl.phe_id, tumor_qtl.top_var))

# credible_sets
tumor_cs = CSV.read("fine_maping_tumor.csv", DataFrame)
dropmissing!(tumor_cs)
tumor_cs = @pipe groupby(tumor_cs, [:phe_id, :cs]) |> combine(_, sdf -> sdf[argmax(sdf.pip), :])
tumor_cs = groupby(tumor_cs, :phe_id)
tumor_cs_dict = Dict(zip([g.phe_id for g in keys(tumor_cs)], [Vector{String}(tumor_cs[g].var_id) for g in keys(tumor_cs)]))

# cell type %
proportions = CSV.read("cell_prop_level_2.csv", DataFrame, missingstring = "NA", normalizenames = true)
dropmissing!(proportions)
select!(proportions, Not(:Column1))
filter!(:tissue => ==("Tumeur"), proportions)
rename!(proportions, :Record_ID => :id)
sort!(proportions, :id)
select!(proportions, Not(:tissue))
select!(proportions, Not(:Submucosal_Gland))
for n in names(proportions)[1:end-1]
    proportions[!,n] = (proportions[!,n] .- mean(proportions[!,n]))./std(proportions[!,n])
end

# regressout

X1 = covars[!,[["age", "sex", "smoking"]; ["PC$i" for i in 1:10]]] |> Matrix
X2 = [X1 Matrix(proportions[!,1:end-1])]

function regressout(y, X::Matrix{<:Real})
    y - X * (X \ y)
end

expr_tumor1 = deepcopy(expr_tumor)
@threads for i in 2:lastindex(expr_tumor1, 2)
    expr_tumor1[!,i] = regressout(expr_tumor1[!,i], X1)
end

expr_tumor2 = deepcopy(expr_tumor)
@threads for i in 2:lastindex(expr_tumor2, 2)
    expr_tumor2[!,i] = regressout(expr_tumor2[!,i], X2)
end

# mutations
mutations = CSV.read("data.LORD.n1603.forSam.csv", DataFrame, select = [:record_id, :EGFR, :TP53, :KRAS])
rename!(mutations, :record_id => :id)
filter!(:id => ∈(expr_tumor.id), mutations)
sort!(mutations, :id)

# genotypes
genotypes = SnpData("TOPMed_chr_all_LORD_pairedRNA.reordered")

## Model

r2_driver = Vector{Float64}(undef, size(expr_tumor, 2)-1)
r2_driver_cell_frac_adj = Vector{Float64}(undef, size(expr_tumor, 2)-1)
r2_top_var = Vector{Float64}(undef, size(expr_tumor, 2)-1)
r2_top_var_cell_frac_adj = Vector{Float64}(undef, size(expr_tumor, 2)-1)
r2_cs = Vector{Float64}(undef, size(expr_tumor, 2)-1)
r2_cs_cell_frac_adj = Vector{Float64}(undef, size(expr_tumor, 2)-1)

adjr2_driver = Vector{Float64}(undef, size(expr_tumor, 2)-1)
adjr2_driver_cell_frac_adj = Vector{Float64}(undef, size(expr_tumor, 2)-1)
adjr2_top_var = Vector{Float64}(undef, size(expr_tumor, 2)-1)
adjr2_top_var_cell_frac_adj = Vector{Float64}(undef, size(expr_tumor, 2)-1)
adjr2_cs = Vector{Float64}(undef, size(expr_tumor, 2)-1)
adjr2_cs_cell_frac_adj = Vector{Float64}(undef, size(expr_tumor, 2)-1)

p = Progress(length(r2_driver))
@threads for i in 1:length(r2_driver)
    j = i+1
    model_driver = lm(@formula(phen ~ 1 + EGFR + KRAS + TP53 + KRAS * TP53 + EGFR * TP53), DataFrame(phen = expr_tumor1[!,j], EGFR = mutations.EGFR, KRAS = mutations.KRAS, TP53 = mutations.TP53))
    r2_driver[i] = r2(model_driver)
    adjr2_driver[i] = clamp(adjr2(model_driver), 0.0, 1.0)

    model_driver_cell_frac_adj = lm(@formula(phen ~ 1 + EGFR + KRAS + TP53 + KRAS * TP53 + EGFR * TP53), DataFrame(phen = expr_tumor2[!,j], EGFR = mutations.EGFR, KRAS = mutations.KRAS, TP53 = mutations.TP53))
    r2_driver_cell_frac_adj[i] = r2(model_driver_cell_frac_adj)
    adjr2_driver_cell_frac_adj[i] = clamp(adjr2(model_driver_cell_frac_adj), 0.0, 1.0)

    # top var
    gene = names(expr_tumor)[j]
    if haskey(tumor_best_var_dict, gene)
        gi = findfirst(==(tumor_best_var_dict[gene]), genotypes.snp_info.snpid)
        X = plink_bed_encoding_to_genotype.(genotypes.snparray[:,[gi]])
        X = X .- mean(X, dims = 1)
        model_top_var = lm(X, expr_tumor1[!,j])
        r2_top_var[i] = r2(model_top_var)
        adjr2_top_var[i] = clamp(adjr2(model_top_var), 0.0, 1.0)

        model_top_var = lm(plink_bed_encoding_to_genotype.(genotypes.snparray[:,[gi]]), expr_tumor2[!,j])
        r2_top_var_cell_frac_adj[i] = r2(model_top_var)
        adjr2_top_var_cell_frac_adj[i] = clamp(adjr2(model_top_var), 0.0, 1.0)
    else
        r2_top_var[i] = 0.0
        adjr2_top_var[i] = 0.0
        r2_top_var_cell_frac_adj[i] = 0.0
        adjr2_top_var_cell_frac_adj[i] = 0.0
    end

    # cs
    if haskey(tumor_cs_dict, gene)
        gi = findall(x -> x ∈ tumor_cs_dict[gene], genotypes.snp_info.snpid)
        X = plink_bed_encoding_to_genotype.(genotypes.snparray[:,gi])
        X = X .- mean(X, dims = 1)
        model_cs = lm(X, expr_tumor1[!,j])
        r2_cs[i] = r2(model_cs)
        adjr2_cs[i] = clamp(adjr2(model_cs), 0.0, 1.0)

        model_cs = lm(X, expr_tumor2[!,j])
        r2_cs_cell_frac_adj[i] = r2(model_cs)
        adjr2_cs_cell_frac_adj[i] = clamp(adjr2(model_cs), 0.0, 1.0)
    else
        r2_cs[i] = 0.0
        adjr2_cs[i] = 0.0
        r2_cs_cell_frac_adj[i] = 0.0
        adjr2_cs_cell_frac_adj[i] = 0.0
    end
    next!(p)
end
finish!(p)

r2_df = DataFrame(ensg_id = names(expr_tumor)[2:end], r2_driver = r2_driver, r2_driver_cell_frac_adj = r2_driver_cell_frac_adj,
                  adjr2_driver = adjr2_driver, adjr2_driver_cell_frac_adj = adjr2_driver_cell_frac_adj,
                  r2_cs = r2_cs, r2_cs_cell_frac_adj = r2_cs_cell_frac_adj,
                  adjr2_cs = adjr2_cs, adjr2_cs_cell_frac_adj = adjr2_cs_cell_frac_adj,
                  r2_top_var = r2_top_var, r2_top_var_cell_frac_adj = r2_top_var_cell_frac_adj,
                  adjr2_top_var = adjr2_top_var, adjr2_top_var_cell_frac_adj = adjr2_top_var_cell_frac_adj)

CSV.write("r2_comparison_v2.csv", r2_df)

## plotting

r2_df = CSV.read("r2_comparison_v2.csv", DataFrame)

using CairoMakie
using Colors

df_with_finemapping = filter(:r2_cs_cell_frac_adj => !=(0), r2_df)
f = Figure(size = (250, 350))
ax = Axis(f[1,1], xticks = ([0, 1], ["EGFR, KRAS\n& TP53", "Credible\nsets"]), yticks = ([-6, -5, -4, -3, -2, -1, 0], ["1e-6", "1e-5", "1e-4", "1e-3", "0.01", "0.1", "1"]),
          ylabel = "gene expression r²", ylabelsize = 17)
boxplot!(ax, zeros(size(df_with_finemapping, 1)), log10.(df_with_finemapping.r2_driver_cell_frac_adj), color = RGB(11/255, 0/255, 51/255), alpha = 0.62)
boxplot!(ax, ones(size(df_with_finemapping, 1)), log10.(df_with_finemapping.r2_cs_cell_frac_adj), color = RGB(11/255, 0/255, 51/255), alpha = 0.62)
hidedecorations!(ax, label = false, ticklabels = false, ticks = false)
f
save("r2_comparison.svg", f)

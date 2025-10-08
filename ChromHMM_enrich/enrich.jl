using JLD2
using CSV
using HypothesisTests
using DataFrames
import Pipe.@pipe
using MultipleTesting: adjust, BenjaminiHochberg
using StatsBase


### LUNG

finemapped_positions = @pipe CSV.read("../qtltools/susie_finemapped_fdr5_eqtls_lung.csv", DataFrame, select = [:var_chr, :var_from]) |> unique
rename!(finemapped_positions, [:chr, :pos])
finemapped_positions.pos = finemapped_positions.pos .- 1


all_tested_pos = CSV.read("universal_chromhmm/all_tested_positions_lung_sorted_chromhmm_uni_unique.bed", DataFrame, 
                          delim = '\t', header = false, select = [1, 2, 8])
rename!(all_tested_pos, [:chr, :pos, :category])

finemapped_positions = innerjoin(finemapped_positions, all_tested_pos, on = [:chr, :pos])

background_category_counts = countmap(all_tested_pos.category)
finemapped_category_counts = countmap(finemapped_positions.category)

out = DataFrame()
for k in keys(background_category_counts)
    p = background_category_counts[k] / size(all_tested_pos, 1)
    n = size(finemapped_positions, 1)
    x = haskey(finemapped_category_counts, k) ? finemapped_category_counts[k] : 0
    test = BinomialTest(x, n, p)
    push!(out, (; category = k, pval = pvalue(test), log2fe = log2((x/n)/p), ci_low = log2(confint(test)[1]/p), ci_high = log2(confint(test)[2]/p)))
end
out.padj = adjust(out.pval, BenjaminiHochberg())

CSV.write("lung_enrichments_chromhmm_universal_binomial.csv", out)

### TUMOR

finemapped_positions = @pipe CSV.read("../qtltools/susie_finemapped_fdr5_eqtls_tumor.csv", DataFrame, select = [:var_chr, :var_from]) |> unique
rename!(finemapped_positions, [:chr, :pos])
finemapped_positions.pos = finemapped_positions.pos .- 1


all_tested_pos = CSV.read("universal_chromhmm/all_tested_positions_tumor_sorted_chromhmm_uni_unique.bed", DataFrame, 
                          delim = '\t', header = false, select = [1, 2, 8])
rename!(all_tested_pos, [:chr, :pos, :category])

finemapped_positions = innerjoin(finemapped_positions, all_tested_pos, on = [:chr, :pos])

background_category_counts = countmap(all_tested_pos.category)
finemapped_category_counts = countmap(finemapped_positions.category)

out = DataFrame()
for k in keys(background_category_counts)
    p = background_category_counts[k] / size(all_tested_pos, 1)
    n = size(finemapped_positions, 1)
    x = haskey(finemapped_category_counts, k) ? finemapped_category_counts[k] : 0
    test = BinomialTest(x, n, p)
    push!(out, (; category = k, pval = pvalue(test), log2fe = log2((x/n)/p), ci_low = log2(confint(test)[1]/p), ci_high = log2(confint(test)[2]/p)))
end
out.padj = adjust(out.pval, BenjaminiHochberg())

CSV.write("tumor_enrichments_chromhmm_universal_binomial.csv", out)


### DELTA

credible_lung = @pipe CSV.read("../qtltools/susie_finemapped_fdr5_eqtls_lung.csv", DataFrame, select = [:var_chr, :var_from]) |> unique
credible_tumor = @pipe CSV.read("../qtltools/susie_finemapped_fdr5_eqtls_tumor.csv", DataFrame, select = [:var_chr, :var_from]) |> unique
fdr5_delta = @pipe load("qtltools/Delta_eQTLs_60peer_fdr5.csv", DataFrame) |> select(_, [:var_chr, :var_from]) |> unique
credible_lung_or_tumor = outerjoin(credible_lung, credible_tumor, on = [:var_chr, :var_from])
finemapped_positions = innerjoin(credible_lung_or_tumor, fdr5_delta, on = [:var_chr, :var_from])
rename!(finemapped_positions, [:chr, :pos])
finemapped_positions.pos = finemapped_positions.pos .- 1

all_tested_pos = CSV.read("universal_chromhmm/all_tested_positions_delta_sorted_chromhmm_uni_unique.bed", DataFrame, 
                          delim = '\t', header = false, select = [1, 2, 8])
rename!(all_tested_pos, [:chr, :pos, :category])

finemapped_positions = innerjoin(finemapped_positions, all_tested_pos, on = [:chr, :pos])

background_category_counts = countmap(all_tested_pos.category)
finemapped_category_counts = countmap(finemapped_positions.category)

out = DataFrame()
for k in keys(background_category_counts)
    p = background_category_counts[k] / size(all_tested_pos, 1)
    n = size(finemapped_positions, 1)
    x = haskey(finemapped_category_counts, k) ? finemapped_category_counts[k] : 0
    test = BinomialTest(x, n, p)
    push!(out, (; category = k, pval = pvalue(test), log2fe = log2((x/n)/p), ci_low = log2(confint(test)[1]/p), ci_high = log2(confint(test)[2]/p)))
end
out.padj = adjust(out.pval, BenjaminiHochberg())

CSV.write("delta_enrichments_chromhmm_universal_binomial.csv", out)

out = nothing
finemapped_positions = nothing
all_tested_pos = nothing
GC.gc()

### figures

out_lung = CSV.read("lung_enrichments_chromhmm_universal_binomial.csv", DataFrame)
out_tumor = CSV.read("tumor_enrichments_chromhmm_universal_binomial.csv", DataFrame)
out_delta = CSV.read("delta_enrichments_chromhmm_universal_binomial.csv", DataFrame)
all(out_lung.category .== out_tumor.category .== out_delta.category)
idx = sortperm(out_delta.log2fe, rev = false)
out_delta = out_delta[idx,:]
out_lung = out_lung[idx,:]
out_tumor = out_tumor[idx,:]


i1 = (out_delta.ci_high .< out_lung.ci_low) .&& (out_delta.ci_high .< out_tumor.ci_low)
i2 = (out_delta.ci_low .> out_lung.ci_high) .&& (out_delta.ci_low .> out_tumor.ci_high)
i = i1 .|| i2

out_delta = out_delta[i,:]
out_lung = out_lung[i,:]
out_tumor = out_tumor[i,:]

out_delta = out_delta[end-19:end,:]
out_lung = out_lung[end-19:end,:]
out_tumor = out_tumor[end-19:end,:]

using CairoMakie
using Colors
using LaTeXStrings

n = size(out_delta, 1)

jlc = Colors.JULIA_LOGO_COLORS

function get_bigcategory(s)
    if occursin("Tx", s)
        return "transcribed"
    elseif occursin("Prom", s)
        return "promoter"
    elseif occursin("Enh", s)
        return "enhancer"
    else
        return "zzz"
    end
end

out_delta.bigcategory = get_bigcategory.(out_delta.category)

out = innerjoin(out_delta, out_lung, on = :category, renamecols = "" => "_lung")
out = innerjoin(out, out_tumor, on = :category, renamecols = "" => "_tumor")
sort!(out, [:bigcategory, :log2fe])


ylabels = @pipe split.(out.category, "_") |> getindex.(_, 2)

f = Figure(size = (600, 300))
ax = Axis(f[1, 1], xticks = (2:4:4*n+1, ylabels), 
          ylabel = L"\log_2(FE)", ylabelsize = 18, 
          xticklabelrotation = π/2, xticklabelsize = 13, yticks = 0.5:0.5:2)

errorbars!(ax, 1:4:4*n, out.log2fe, 
          (out.log2fe-out.ci_low), direction = :y,
          color = :grey)
errorbars!(ax, 2:4:4*n+1, out.log2fe_lung, 
          (out.log2fe_lung-out.ci_low_lung), direction = :y,
          color = :grey)
errorbars!(ax, 3:4:4*n+2, out.log2fe_tumor, 
          (out.log2fe_tumor-out.ci_low_tumor), direction = :y,
          color = :grey)
sd = scatter!(ax, 1:4:4*n, out.log2fe , marker = :diamond, 
          markersize = 8, color = :red3)
sl = scatter!(ax, 2:4:4*n+1, out.log2fe_lung, marker = :circle, 
          markersize = 8, color = :black)
st = scatter!(ax, 3:4:4*n+2, out.log2fe_tumor, marker = :utriangle, 
          markersize = 8, color = :black)
vlines!(ax, 4:4:4*n+3, color = :black, alpha = 0.5)
xlims!(ax, 0, 4*n)
hidexdecorations!(ax, ticks = false, ticklabels = false)

Legend(f[1, 2], [sd, sl, st], ["Δ", "Lung", "Tumor"])
f

save("img/top_20_delta_chromhmm_non_overlap_ci.pdf", f)
save("img/top_20_delta_chromhmm_non_overlap_ci.png", f, px_per_unit = 2)
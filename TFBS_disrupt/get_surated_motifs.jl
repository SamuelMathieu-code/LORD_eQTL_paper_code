using CSV
using DataFrames
using Distributions
using JLD2
using DelimitedFiles
import Pipe.@pipe
using MultipleTesting
import StatsBase.countmap

### filter out tfs expressed under 1 tpm (?)

    expression_lung = CSV.read("../lung_tpms.tsv", DataFrame, delim = "\t")
    expression_tumor = CSV.read("../tumor_tpms.tsv", DataFrame, delim = "\t")

    perfectos = CSV.read("merge_output_perfectosAPE.tsv", DataFrame, delim = "\t")

    all_motifs = unique(perfectos.motif)
    all_tfs = @pipe split.(all_motifs, ".") |> getindex.(_, 4)

    motifs_df = DataFrame(motif = all_motifs, tfs = all_tfs)

    symbol_to_ensg = CSV.read("ENSG_id_version_to_gene_symbol", DataFrame) # dataframe with ensembl id and hugo gene name
    symbol_to_ensg = Dict(zip( symbol_to_ensg.hgnc_symbol,  symbol_to_ensg.ensembl_gene_id))
    motifs_df.tfs_ensgs = [[haskey( symbol_to_ensg, gi) ?  symbol_to_ensg[gi] : "." for gi in split(g, ':', keepempty = false)]  for g in motifs_df.tfs]

    countmap(reduce(vcat, motifs_df.tfs_ensgs))["."] # 1 (one non-mapped TF name)

    function keep_tf(ensgs) # verify if gene is expressed at 1 TPM or + in at least Lung OR Tumor
        ans = Vector{Bool}(undef, size(ensgs, 1))
        for i in eachindex(ans)
            j = findfirst(x -> occursin(ensgs[i], x), names(expression_lung))
            if ensgs[i] == "." || isnothing(j)
                ans[i] = false
            else
                if (mean(expression_lung[:,j]) ≥ 1) | (mean(expression_tumor[:,j]) ≥ 1)
                    ans[i] = true
                else
                    ans[i] = false
                end
            end
            if isnothing(j)
                ans[i] = true
            end
        end

        return reduce(&, ans)
    end

    motifs_df.keep = keep_tf.(motifs_df.tfs_ensgs)
    
    filter!(:keep => identity, motifs_df)
    CSV.write("expressed_motifs.csv", select(motifs_df, [:motif, :tfs]))
    
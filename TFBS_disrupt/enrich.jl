using CSV
using DataFrames
using Distributions
using DelimitedFiles
using MultipleTesting: adjust, BenjaminiHochberg, PValues
using HypothesisTests
import Pipe.@pipe

# load data

    # delta
    delta = CSV.read("../qtltools/Delta_eQTLs_60peer_fdr5.jld2", DataFrame)

    # all test lung and tumor
    lung_all = CSV.read("../qtltools/Lung_eQTLs_60_peer.tsv", DataFrame, delim = '\t', select = [:phe_id, :var_id, :slope])
    tumor_all = CSV.read("../qtltools/Tumor_eQTLs_60_peer.tsv", DataFrame, delim = '\t', select = [:phe_id, :var_id, :slope])

    # credible sets.
    lung_cs = CSV.read("../qtltools/susie_finemapped_fdr5_eqtls_lung.csv", DataFrame)
    tumor_cs = CSV.read("../qtltools/susie_finemapped_fdr5_eqtls_tumor.csv", DataFrame)
    lung_cs = groupby(lung_cs, [:phe_id, :cs])
    tumor_cs = groupby(tumor_cs, [:phe_id, :cs])

# equalize credible set sizes

    function equalize_cs_size(d, p::Real)
        idx = partialsortperm(d.pip, 1:max(1, Int(round(p*size(d, 1)))), rev = true)
        return d[idx, :]
    end

    mean_cs_size_lung = [size(d, 1) for d in lung_cs] |> mean
    mean_cs_size_tumor = [size(d, 1) for d in tumor_cs] |> mean

    # mean size of lung is 1.48 times lower.
    # As there are more variants in tumor cs, the odds of at least one variant in cs disrupting a TF motif are much higher for tumor.
    # To account for this bias, we equalize the credible sets sizes by keeping ∼2/3 of the best variants in each tumor credible set.

    p = mean_cs_size_lung/mean_cs_size_tumor

    lung_cs = combine(lung_cs, x -> equalize_cs_size(x, 1.0))
    tumor_cs = combine(tumor_cs, x -> equalize_cs_size(x, p))

    [size(d, 1) for d in groupby(lung_cs, [:phe_id, :cs])] |> mean
    # ↪ 30.02
    [size(d, 1) for d in groupby(tumor_cs, [:phe_id, :cs])] |> mean
    # ↪ 30.02

    # Perfectos
    perfectos = CSV.read("merge_output_perfectosAPE.tsv", DataFrame, delim = "\t")
    filter!(:fold_change => x -> (x<1/5) || (x>5), perfectos)
    filter!([:pvalue_1, :pvalue_2] => (x, y) -> (x<1e-5) || (y<1e-5), perfectos)

# Delta-eQTLs enrichment

    delta_vars_set = Set(unique(delta.var_id))

    function is_lead_delta(nom_pval, var_id)
        i = argmin(nom_pval)
        return var_id[i] ∈ delta_vars_set
    end

    tfbs_disrupt_events_lung = combine(lung_cs_g, :motif_id => is_motif_disrupted => AsTable, [:nom_pval, :var_id] => is_lead_delta => :is_delta)
    tfbs_disrupt_events_tumor = combine(tumor_cs_g, :motif_id => is_motif_disrupted => AsTable, [:nom_pval, :var_id] => is_lead_delta => :is_delta)

    N_delta = count(tfbs_disrupt_events_lung.is_delta) + count(tfbs_disrupt_events_tumor.is_delta)

    x_successes_delta = Vector{Int}(undef, size(clusters, 1))
    log_fold_enrich_delta = Vector{Float64}(undef, size(clusters, 1))
    pvalue_delta = Vector{Float64}(undef, size(clusters, 1))
    
    for (i, tf) in enumerate(clusters.motif_id)
        x_successes_delta[i] = count(tfbs_disrupt_events_lung[!,tf] .&& tfbs_disrupt_events_lung.is_delta) + 
                                count(tfbs_disrupt_events_tumor[!,tf] .&& tfbs_disrupt_events_tumor.is_delta)
        d_delta = Hypergeometric(x_successes_lung[i]+x_successes_tumor[i], N_lung+N_tumor-x_successes_lung[i]-x_successes_tumor[i], N_delta)
        p_obs_delta = x_successes_delta[i]/N_delta
        p = (x_successes_lung[i]+x_successes_tumor[i])/(N_lung+N_tumor)
        pvalue_delta[i] = 1 - cdf(d_delta, x_successes_delta[i]-1)
        log_fold_enrich_delta[i] = log2(p_obs_delta/p)
    end

    out_delta = DataFrame(motif_id = clusters.motif_id,
                        x_successes_delta = x_successes_delta,
                        N_delta = N_delta,
                        log_fold_enrich = log_fold_enrich_delta,
                        pval = pvalue_delta, 
                        padj = adjust(PValues(pvalue_delta), BenjaminiHochberg()))
    out_delta = leftjoin(clusters, out_delta, on = :motif_id)
    sort!(out_delta, :pval)

    CSV.write("output_enrich_delta_all_v6.csv", out_delta)

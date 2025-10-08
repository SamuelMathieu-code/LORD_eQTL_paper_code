using CSV
using DataFrames
using Folds
import Pipe.@pipe

files = @pipe [CSV.read("out_delta_60/out_$i.txt", DataFrame, delim = " ", header = 
        [:phe_id, :phe_chr, :phe_from, :phe_to, :phe_strd, :n_var_in_cis, 
        :dist_phe_var, :var_id, :var_chr, :var_from, :var_to, :nom_pval, 
        :r_squared, :slope, :slope_se, :best_hit]) for i in 1:40] |> Folds.reduce(vcat, _, init = DataFrame())

CSV.write("Delta_eQTLs_60_peer.tsv", files, delim = "\t")

files = @pipe [CSV.read("out_lung_60/out_$i.txt", DataFrame, delim = " ", header = 
        [:phe_id, :phe_chr, :phe_from, :phe_to, :phe_strd, :n_var_in_cis, 
        :dist_phe_var, :var_id, :var_chr, :var_from, :var_to, :nom_pval, 
        :r_squared, :slope, :slope_se, :best_hit]) for i in 1:40] |> Folds.reduce(vcat, _, init = DataFrame())

CSV.write("Lung_eQTLs_60_peer.tsv", files, delim = "\t")

files = @pipe [CSV.read("out_tumor_60/out_$i.txt", DataFrame, delim = " ", header = 
        [:phe_id, :phe_chr, :phe_from, :phe_to, :phe_strd, :n_var_in_cis, 
        :dist_phe_var, :var_id, :var_chr, :var_from, :var_to, :nom_pval, 
        :r_squared, :slope, :slope_se, :best_hit]) for i in 1:40] |> Folds.reduce(vcat, _, init = DataFrame())

CSV.write("Tumor_eQTLs_60_peer.tsv", files, delim = "\t")

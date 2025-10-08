using CSV
using DataFrames
using JLD2
using StatsBase

pos = CSV.read("../qtltools/Tumor_eQTLs_60_peer.tsv", DataFrame, delim = "\t")
select!(pos, [:var_chr, :var_from, :var_to, :phe_id])
pos.var_from .-= 1
pos.var_to = pos.var_from .+ 1
CSV.write("all_tested_positions_tumor.bed", pos, delim = "\t", header = false)
GC.gc()

pos = CSV.read("../qtltools/Lung_eQTLs_60_peer.tsv", DataFrame, delim = "\t")
select!(pos, [:var_chr, :var_from, :var_to, :phe_id])
pos.var_from .-= 1
pos.var_to = pos.var_from .+ 1
CSV.write("all_tested_positions_lung.bed", pos, delim = "\t", header = false)
GC.gc()

pos = CSV.read("../qtltools/Delta_eQTLs_60_peer.tsv", DataFrame, delim = "\t")
select!(pos, [:var_chr, :var_from, :var_to, :phe_id])
pos.var_from .-= 1
pos.var_to = pos.var_from .+ 1
CSV.write("all_tested_positions_delta.bed", pos, delim = "\t", header = false)

using DataFrames
using CSV
using MultipleTesting


eqtls_t = CSV.read("Tumor_eQTLs_60_peer.tsv", DataFrame, delim="\t", 
                   types = Dict(:var_id => String127, :var_chr => String7, :var_from => Int32),
                   pool = Dict(:phe_id => true, :var_chr => true))

eqtls_s = CSV.read("Lung_eQTLs_60_peer.tsv", DataFrame, delim="\t", 
                   types = Dict(:var_id => String127, :var_chr => String7, :var_from => Int32),
                   pool = Dict(:phe_id => true, :var_chr => true))

eqtls_d6 = CSV.read("Delta_eQTLs_60_peer.tsv", DataFrame, delim="\t", 
                    types = Dict(:var_id => String127, :var_chr => String7, :var_from => Int32),
                    pool = Dict(:phe_id => true, :var_chr => true))

eqtls_s.padj = adjust(PValues(eqtls_s.nom_pval), BenjaminiHochberg())
eqtls_t.padj = adjust(PValues(eqtls_t.nom_pval), BenjaminiHochberg())
eqtls_d6.padj = adjust(PValues(eqtls_d6.nom_pval), BenjaminiHochberg())
filter!(:padj => <(5e-2), eqtls_s)
filter!(:padj => <(5e-2), eqtls_t)
filter!(:padj => <(5e-2), eqtls_d6)

CSV.write("Tumor_eQTLs_fdr5.csv", eqtls_t)
CSV.write("Lung_eQTLs_fdr5.csv", eqtls_s)
CSV.write("Delta_eQTLs_60peer_fdr5.csv", eqtls_d4)

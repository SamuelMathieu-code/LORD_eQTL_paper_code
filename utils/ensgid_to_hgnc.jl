using DataFrames
using CSV
import Pipe.@pipe

function convert_to_hgnc(x; isversion = true)
   dat = CSV.read("../utils/ENSG_id_version_to_gene_symbol", DataFrame, delim = " ")
   if isversion
      x2 = @pipe split.(x, ".") |> getindex.(_, 1) 
   else
      x2 = x
   end
   dat = Dict(zip(dat.ensembl_gene_id, dat.hgnc_symbol))
   return [(haskey(dat, g)) ? dat[g] : g for g in x2]
end
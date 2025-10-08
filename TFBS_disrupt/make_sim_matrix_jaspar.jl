using Distributed
addprocs(20, lazy = false)
using Dagger
using Pkg
using DataFrames
using CSV
import Pipe.@pipe
using LinearAlgebra
using DelimitedFiles
import IterTools.subsets
ENV["R_HOME"] = "miniconda3/envs/pwmenrich_r/lib/R"
Pkg.build("RCall")


@everywhere begin
    using RCall
end

@everywhere begin
    R"library(PWMEnrich)"
end

@everywhere function get_similarity(m1::Matrix{Int}, m2::Matrix{Int})::Float64
    @rput m1 m2
    R"rownames(m1) <- c(\"A\", \"C\", \"G\", \"T\")"
    R"rownames(m2) <- c(\"A\", \"C\", \"G\", \"T\")"
    R"s <- motifSimilarity(m1, m2)"
    return @rget s
end

tfs_df = CSV.read("expressed_motifs.csv", DataFrame)
output_mat_tasks = Matrix{Dagger.DTask}(undef, size(tfs_df, 1), size(tfs_df, 1))
motifs = @pipe split.(tfs_df.motif, " ") |> getindex.(_, 1)

n = binomial(length(motifs), 2)

k = 1
for (i, j) in subsets(eachindex(motifs), 2)
    print("$(k)/$(n)")
    m1 = @pipe DelimitedFiles.readdlm("motif_repo/motif_$(motifs[i]).pfm",
                                        ' ', Float64, '\n', skipstart=1).*100 |> round.(_) |> Int.(_)
    m2 = @pipe DelimitedFiles.readdlm("motif_repo/motif_$(motifs[j]).pfm",
                                        ' ', Float64, '\n', skipstart=1).*100 |> round.(_) |> Int.(_)

    @inbounds output_mat_tasks[i, j] = output_mat_tasks[j, i] = @Dagger.spawn get_similarity(m1, m2)

    print("\r")
    k+=1
end
println()
println("done")

for i in eachindex(motifs)
    output_mat_tasks[i, i] = Dagger.@spawn identity(1.0)
end

output_mat = fetch.(output_mat_tasks)

CSV.write("expressed_tfs_motif_similrity_matrix.csv", Tables.table(output_mat, header = motifs))

include("../plink_to_genotype.jl")
using CSV
using DataFrames
using Survival
using SnpArrays
using StatsModels
using StatsBase
import Base.Threads.@threads
using ProgressMeter

############### PFS ###############
surv_dat = CSV.read("input_data/PFS_NSCLC_and_covariates_clean_data.csv", DataFrame; types = Dict("Record_ID" => String), normalizenames = true)
genotypes = SnpData("full_cohort_genotypes")
index_persons = findfirst.( .==(surv_dat[!, "Record_ID"]), Ref(genotypes.person_info.iid))
surv_dat.event_time = EventTime.(surv_dat.PFS_2024_, surv_dat.event)

g = DataFrame(g = plink_bed_encoding_to_genotype.(genotypes.snparray[index_persons, 1]))
f2 = @formula(event_time ~  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + Smoking_status + Sex_at_birth + age_year + Pathologic_stade_9th_edition + Histological_type + g)
f2 = apply_schema(f2, schema(f2, hcat(surv_dat, g)))

function make_cox(g_index::Int64; tol = 1e-4)
    g = DataFrame(g = plink_bed_encoding_to_genotype.(genotypes.snparray[index_persons, g_index]))
    model = try
        coxph(f2, hcat(surv_dat, g), tol = tol)
    catch e
        nothing
    end
    if !isnothing(model)
        b, se = model.model.β[end], sqrt(model.model.vcov[end, end])
    else
        b, se = NaN, NaN
    end
    return b, se
end

β = Vector{Float64}(undef, size(genotypes.snp_info, 1))
se = Vector{Float64}(undef, size(genotypes.snp_info, 1))

make_cox(1; tol = 1e-2) # precompile

p = Progress(size(genotypes.snp_info, 1))
@threads for i in 1:lastindex(genotypes.snp_info, 1)
    β[i], se[i] = make_cox(i; tol = 1e-2)
    next!(p)
end

CSV.write("beta_se_est_PFS.csv", DataFrame(snpid = genotypes.snp_info.snpid, other_allele = genotypes.snp_info.allele1, effect_allele = genotypes.snp_info.allele2, beta = β, se = se))

############### OS ###############
surv_dat = CSV.read("input_data/OS_NSCLC_and_covariates_clean_data.csv", DataFrame; types = Dict("Record_ID" => String), normalizenames = true)
surv_dat.event_time = EventTime.(surv_dat.time_to_event, surv_dat.event)

g = DataFrame(g = plink_bed_encoding_to_genotype.(genotypes.snparray[index_persons, 1]))
f2 = @formula(event_time ~  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + Smoking_status + Sex_at_birth + age_year + Pathologic_stade_9th_edition + Histological_type + g)
f2 = apply_schema(f2, schema(f2, hcat(surv_dat, g)))

function make_cox(g_index::Int64; tol = 1e-4)
    g = DataFrame(g = plink_bed_encoding_to_genotype.(genotypes.snparray[index_persons, g_index]))
    model = try
        coxph(f2, hcat(surv_dat, g), tol = tol)
    catch e
        nothing
    end
    if !isnothing(model)
        b, se = model.model.β[end], sqrt(model.model.vcov[end, end])
    else
        b, se = NaN, NaN
    end
    return b, se
end

β = Vector{Float64}(undef, size(genotypes.snp_info, 1))
se = Vector{Float64}(undef, size(genotypes.snp_info, 1))

make_cox(1; tol = 1e-2) # precompile

p = Progress(size(genotypes.snp_info, 1))
@threads for i in 1:lastindex(genotypes.snp_info, 1)
    β[i], se[i] = make_cox(i; tol = 1e-2)
    next!(p)
end

CSV.write("beta_se_est_OS.csv", DataFrame(snpid = genotypes.snp_info.snpid, other_allele = genotypes.snp_info.allele1, effect_allele = genotypes.snp_info.allele2, beta = β, se = se))


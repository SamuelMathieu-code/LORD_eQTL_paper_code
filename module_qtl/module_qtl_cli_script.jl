using CSV
using DataFrames
using DelimitedFiles
using GLM
using Dagger
using MultivariateStats
using SnpArrays
import IterTools.product
using ArgParse
using Distributions
using HypothesisTests
import Pipe.@pipe

s = ArgParseSettings()
@add_arg_table s begin
    "--covariates", "-c"
        help = "csv/tsv file containing covariates. With a header (covars names) and first colum sample ids named \"id\"."
        required = true
    "--snplist", "-s"
        help = "text file containing list of snps to be tested"
        required = true
    "--modules", "-m"
        help = "csv/tsv file containing gene ids in first column and module number in the second (must have exactly two columns)."

    "--expression", "-e"
        help = "csv/tsv file containing expression tpms. With a header (gene ids) and first column sample ids named \"id\"."

    "--eigengenes", "-y"
        help = "csv/tsv file containing eigen genes. With header (module ids) and first column sample ids named \"id\""

    "--genotypes", "-g"
        help = "path to .bed .bim .fam files with prefix."
        required = true
    "--noorderquantnorm"
        help = "Flag : do not perform order quantile normalization over eigengene values before regressions."
        action = :store_true
    "--outputfile"
        help = "file to write results table in. (Will create or overwrite.)"
        required = true
    "--logfile"
        help = "file to log informations about the analysis in real time and potential errors"

    "--writeeigengenes"
        help = "file to write eigen_genes (before ordered quantile normalization). (Will create or overwrite.)"
end
global args = parse_args(ARGS, s, as_symbols = false)


"""
    calceigengenes(expression::DataFrame, modules::DataFrame)

calculate eigengenes for each module
"""
function calceigengenes(expression::DataFrame, modules::DataFrame)::DataFrame
    n_modules = length(unique(modules.modules))
    modules_g = groupby(modules, :modules)
    eigengenes = Matrix{Float64}(undef, size(expression, 1), n_modules)
    moduleno = Vector{String31}(undef, n_modules)
    @inbounds for (i, (k, g)) in collect(enumerate(pairs(modules_g)))
        expr_module = expression[:, g.genes]
        pca_dat = Matrix(expr_module)'
        pc_first = @pipe fit(PCA, pca_dat, maxoutdim = 1) |> predict(_, pca_dat) |> vec
        eigengenes[:, i] = pc_first
        moduleno[i] = "mod"*string(k.modules)
    end
    out = DataFrame(eigengenes, moduleno)
    return insertcols(out, 1, "id" => expression.id)
end


function orderquantilenorm(v)
    quantile(Normal(), (sortperm(v) .- 0.5) / length(v))
end


"""
    plink_bed_encoding_to_genotype(x::UInt8)

convert plink bed format encoding in number of alt alleles carried.
"""
function plink_bed_encoding_to_genotype(x::UInt8)::Float16
    if x == 0x00
        return 0.0
    elseif x == 0x01
        return NaN
    else
        return x - 1.0
    end
end


"""
    oneregression(covars::DataFrame, genotypes::AbstractVector{UInt8}, var::AbstractString, mod::AbstractString, oth::AbstractString, eff::AbstractString, eigengene::Vector{Float64})

Make regression of one snp on eigengene.

# Arguments:
- `covars`: datafrane of the covariates to include
- `genotypes`: vector of genotypes in plink .bed format (UInt8)
- `var`: variant id
- `mod`: module id
- `oth`: other allele (encoded by 0x00 homozygote)
- `eff`: effect allele (encoded by 0x03 homozygote)
- `eigengene`: Vector containing values of the eigengene.
"""
function oneregression(covars::DataFrame, genotypes::AbstractVector{UInt8}, var::AbstractString, mod::AbstractString, oth::AbstractString, eff::AbstractString, eigengene::AbstractVector{<:AbstractFloat}; sw::Bool = false)::NamedTuple
    data = DataFrame(covars, copycols = false)
    data.y = eigengene
    data.g = plink_bed_encoding_to_genotype.(genotypes)
    f1 = Term(:y) ~ sum(Term.(Symbol.(names(covars)[2:end]))) + Term(:g)
    m1 = lm(f1, data)
    if sw
        swp = residuals(m1) |> ShapiroWilkTest |> pvalue
        return (var = var, mod = mod, oth = oth, eff = eff, β = coef(m1)[end], se = stderror(m1)[end], pval = coeftable(m1).cols[4][end], tval = coeftable(m1).cols[3][end], shapirowilkp = swp)
    else
        return (var = var, mod = mod, oth = oth, eff = eff, β = coef(m1)[end], se = stderror(m1)[end], pval = coeftable(m1).cols[4][end], tval = coeftable(m1).cols[3][end])
    end
end


"""
    moduleregressions(genotypes::SnpData, variants::Vector{<:AbstractString}, covariates::DataFrame, eigengenes::DataFrame)

make regression for all pairs of snp - module

# Arguments:
- `genotypes`: SnpData for plink files
- `variants`: Vector of variant ids to use
- `covariates`: DataFrame of covariates (samples in rows)
- `eigengenes`: DataFrame of eigengenes (samples in rows)
"""
function moduleregressions(genotypes::SnpData, variants::Vector{<:AbstractString}, covariates::DataFrame, eigengenes::DataFrame; sw::Bool = false)::DataFrame
    @info "searching variants in .bim file..."
    varidxs = indexin(variants, genotypes.snp_info.snpid)
    @info "$(count(isnothing, varidxs)) missing variants in .bim file"

    @info "making regressions..."

    allpairs = product(names(eigengenes)[2:end], zip(variants, varidxs))
    nmods, nvars = size(allpairs)
    ntests = nmods * nvars
    outtasks = Vector{Dagger.DTask}(undef, ntests)

    orderpersons = [findfirst(==(i), genotypes.person_info.iid) for i in covariates.id]

    if any(isnothing(orderpersons))
        throw(ArgumentError(), "Some sample are missing record in plink .fam file")
    end
    
    for (i, (mod, (var, idx))) in enumerate(allpairs)
        oth, eff = genotypes.snp_info.allele1[idx], genotypes.snp_info.allele2[idx]
        outtasks[i] = Dagger.@spawn oneregression(covariates, @view(genotypes.snparray[orderpersons, idx]), var, mod, oth, eff, @view(eigengenes[:,mod]), sw = sw)
    end
    out = DataFrame()
    push!.(Ref(out), fetch.(outtasks))
    return out
end


function modqtls()
    @info "loading data..."
    
    #covariates
    covariates = CSV.read(args["covariates"], DataFrame, types = Dict(1 => String))
    rename!(covariates, 1 => :id)
    sort!(covariates, :id)

    #modules
    if args["modules"] !== nothing
        modules = CSV.read(args["modules"], DataFrame)
        rename!(modules, [:genes, :modules])
        filter!(:modules => !=(0), modules)
    end

    # SNP list to try
    snplist = readdlm(args["snplist"], String) |> vec

    #expression
    if args["expression"] !== nothing
        expressiondf = CSV.read(args["expression"], DataFrame, select = ["id"; modules.genes], types = Dict(1 => String))
        rename!(expressiondf, 1 => :id)
        sort!(expressiondf, :id)
        expressiondf[:,2:end] = log.(expressiondf[:,2:end] .+ 0.5)
    end

    #genotypes
    genotypes = SnpData(args["genotypes"])

    # get eigen gene values (read or calculate)
    if args["eigengenes"] === nothing
        @info "calculating eigen genes..."
        eigengenes = calceigengenes(expressiondf, modules)
        if !isnothing(args["writeeigengenes"])
            @info "writing eigen genes..."
            CSV.write(args["writeeigengenes"], eigengenes)
        end
    else
        eigengenes = CSV.read(args["eigengenes"], DataFrame, types = Dict(1 => String))
        sort!(eigengenes, :id)
    end

    if eigengenes.id != covariates.id
        throw(ArgumentError(), "sample ids are not the same in expression and covariates data")
    end

    # ordered quantile normalization
    if !args["noorderquantnorm"]
        for col in names(eigengenes)[2:end]
            eigengenes[:, col] = orderquantilenorm(eigengenes[:, col])
        end
    end

    outtable = moduleregressions(genotypes, snplist, covariates, eigengenes, sw = args["noorderquantnorm"])
    
    @info "writing output table..."
    CSV.write(args["outputfile"], outtable)

    @info "done."
end


function main()
    if ((args["expression"] === nothing) | (args["modules"] === nothing)) & (args["eigengenes"] ===  nothing)
        throw(ArgParseError("Expected --expression and --modules or --eigengenes."))
    end
    if isnothing(args["logfile"])
        modqtls()
    else
        open(args["logfile"], "w") do io
            redirect_stderr(modqtls, io)
        end
    end
end

main()
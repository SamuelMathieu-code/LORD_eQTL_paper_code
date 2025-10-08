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
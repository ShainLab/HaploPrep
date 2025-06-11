# Variant Caller Testing

We tested our application with HaplotypeCaller as an alternative to the default FreeBayes variant caller. 
Simple code modifications were required to accommodate differences in VCF headers and INFO structure between the two callers.

## Files

- `FB-InformativeSNPs.vcf` - Output using FreeBayes variant caller
- `HC-InformativeSNPs.vcf` - Output using HaplotypeCaller variant caller  
- `variantcallercomparison.pdf` - Summary of metrics comparing both callers

## Notes

- Both VCF files contain the same input dataset processed through these different variant callers
- Final phasing step was not performed as it would not affect the outcome of this analysis

## Results Summary

FreeBayes: 15,384 final informative SNPs
HaplotypeCaller: 15,520 final informative SNPs
Overlapping variants: 15,176 SNPs (97.8% concordance)
Variants unique to FreeBayes: 191 SNPs (<1.24% of total variants)
Variants unique to HaplotypeCaller: 344 SNPs (<2.2% of total variants)

See `variantcallercomparison.pdf` for detailed metrics.

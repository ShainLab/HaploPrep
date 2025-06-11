# Variant Caller Testing

To verify that the user could use alternative approaches to come up with an initial candidate list of SNPs, we tested our pipeline using HaplotypeCaller as an alternative to the default FreeBayes variant caller on a representative sample.

## Files

- `FB-InformativeSNPs.vcf` - Output using FreeBayes variant caller
- `HC-InformativeSNPs.vcf` - Output using HaplotypeCaller variant caller  
- `variantcallercomparison.pdf` - Summary of metrics comparing both callers

## Notes

- Both VCF files contain the same input dataset processed through these different variant callers
- Simple code modifications were required to accommodate differences in VCF headers and INFO structure between the two callers tested
- Final phasing step was not performed as it would not affect the outcome of this analysis

## Results Summary

FreeBayes: 15,384 final informative SNPs  
HaplotypeCaller: 15,520 final informative SNPs  
Overlapping variants: 15,176 SNPs **(97.8% concordance)**  
Variants unique to FreeBayes: 191 SNPs **(<1.24% of total variants)**  
Variants unique to HaplotypeCaller: 344 SNPs **(<2.2% of total variants)**  

See `variantcallercomparison.pdf` for detailed metrics.

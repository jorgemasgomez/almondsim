# Sequential genomic OPR

This replaces the previous coancestry-endpoint optimization design.

| GS-2 early1500 scenario | Maximum annual diversity loss |
|---|---:|
| OPR Diversity | 5% |
| OPR Balanced | 10% |
| OPR Gain | 20% |

All three use 32 parents, 20 crosses of 75 offspring, and genomic seedling selection from 1,500 to 500. Parental candidates are HPT2 individuals, age five, with two phenotyping years and target aggregate h²=0.2.

## Annual sequence

1. Predict GEBVs for current parents and candidates using the same genomic model.
2. Measure initial parental diversity: mean expected heterozygosity `mean(2*p*(1-p))` across all SNP-chip markers. Frequencies p come from the 32 current parents. Keep monomorphic SNPs in the denominator.
3. Restrict the candidate list to the top 50 distinct new individuals from the 500 eligible candidates, ranked by descending GEBV. Never proceed to candidate 51 after a rejection.
4. Test each candidate sequentially. Temporarily merge it with the current parents and retain the 32 highest-GEBV individuals. For a single candidate this means replacing the lowest-GEBV current parent, provided the candidate is strictly better. Ties retain the existing parent.
5. Measure diversity of that proposed group. Accept only if it is at least 95%, 90%, or 80% of the START-OF-YEAR diversity, depending on the scenario. After acceptance, the updated parent group is used for the next candidate, but the diversity reference does not change during that year.
6. Continue until all 32 initial parents have been replaced or the top-50 list is exhausted. There is no eight-parent quota. A rejected or discarded candidate never counts as an incorporation.
7. Generate the next cohort using the existing random compatible-cross procedure.

`opr_max_diversity_loss` is 0.05, 0.10, or 0.20. `opr_candidate_limit` is 50. Renewal is free and conditional, from zero to 32 individuals; other existing scenarios retain their exact renewal rules.

Example: initial diversity He=0.30 gives minimum He=0.285 for Diversity, 0.27 for Balanced, and 0.24 for Gain. The reference resets next year. These are annual group-diversity limits, not fixed multi-year retention targets or inbreeding percentages.

The diversity criterion uses observed SNP dosages, not simulated true genetic variance. It preserves SNP heterozygosity subject to the threshold; it does not guarantee preservation of every rare allele or the character's additive variance, and introduces no external alleles.

## Recorded outputs

`oprMaxDiversityLoss`, `opr_beforeHe`, `opr_afterHe`, `opr_minHe`, `opr_lossPercent`, `opr_expectedBV`, `opr_tested`, `opr_rejectedDiversity`, `newParents`, and `realizedRenewalPercent`. `renewalRule` is `OPR_free`. Negative loss means diversity increased.

Validation passed: GEBV order; strict top-50 cutoff even when candidate 51 would be acceptable; annual cumulative diversity limits; acceptance and rejection; free renewal beyond eight and up to the full parental pool; no-op and zero-quota cases; existing scheme regressions; and two future years per OPR variant through a reduced AlphaSimR pipeline with early selection and compatible crossing. No full OPR simulations have been launched.

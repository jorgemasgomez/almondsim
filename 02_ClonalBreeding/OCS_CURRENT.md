# Current OCS implementation (19 October 2026)

Diversity, Balanced and Gain use Ne control parameters 30, 10 and 5. These are constraint proxies, not estimates of realized effective population size.

Native optiSel max.BV/default solver uses one variable per genotype. Mutually incompatible S-genotype groups have combined contribution <=50%; sharing only one S allele permits crossing. Self-compatible carriers are treated separately.

GS6 retains the individual combined contribution cap of 10%. Speed breeding uses a 3.2% adult cap (age >=4), a 10% young cap, and at most 50% contribution from pollen-only candidates. This guarantees at least 16 adult contributors, verified after rounding/mating; it does not set a maximum number of parents. Mothers must be >=4, pollen donors >=2 in speed and >=4 otherwise. SI, no selfing and 1500 offspring followed by early retention of 500 remain hard constraints.

Continuous and rounded coancestry excess up to 10% relative to the bound plus 1e-7 is audited and permitted. Native offspring counts remain fixed during mating allocation. No repair, refinement or automatic retry. Failed inputs and native fits are captured for diagnosis.

The 10-replicate completion run reuses GS6 replica1 from v11 and speed replica1 from v12 without relabelling their designVersion. Replicas2..10 use the exact respective frozen source. No old v10 trajectories or failed speed-v11 outcomes are pooled.

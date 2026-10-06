# GS-2 + OCS: Diversity

Reference: GS_2_Early1500 (1,500 offspring -> 500 seedlings selected by GS). Eligible parents are current
breeders plus HPT2, age five, after two phenotyping years (aggregate h2=0.20).
OCS chooses identities and realized contributions freely, allowing complete
renewal and at most 32 active parents. No obligatory replacement quota.

20 crossing events x 75 offspring = 1,500, identical to reference. Repeated
parental pairs are permitted, as in random crossing. Selfing and completely
S-incompatible crosses are forbidden. Each family has exactly 75 survivors.

Bound: K_min + 0.05*(K_gain-K_min). The endpoints are recalculated from the
same eligible pool annually: minimum coancestry found by local search and the
best compatible pair by EBV. This is an adaptive trade-off restriction, NOT a
fixed annual inbreeding rate or a target Ne. Lower alpha favors diversity.

K is half the centered genomic relationship matrix with SNP frequencies fixed
at the founders. This score can be negative for pairwise relationships and is
not an absolute IBD probability. Realized plan coancestry is checked.

The optimizer is a multistart integer local-search approximation, not a proof
of the global OCS optimum. It considers every eligible candidate in moves,
keeps the at-most-32 constraint and verifies the bound after all assignments.
Mate allocation minimizes kinship further while preserving contributions.

Provisional calibration uses existing almond genotypes; repeat on the new
32-founder burn-in when bank150 is complete. Run 00RUNME.R alone or main.R
for all 22 schemes. Twelve replicates by default; none launched here.

Early culling uses GS predictions before phenotyping. OCS acts on parental
selection and contributions at HPT2; it does not replace seedling GS selection.

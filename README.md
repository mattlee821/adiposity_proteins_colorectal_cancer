[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7636326.svg)](https://doi.org/10.5281/zenodo.7636326)

# Adiposity, proteins, and colorectal cancer risk: Mendelian randomization analyses

We performed a series of two-sample summary level univariable (UV) and
multivariable (MV) Mendelian randomization (MR) analyses to investigate
whether adiposity-associated proteins mediate the relationship between
adiposity and colorectal cancer (CRC).

All scripts were run sequentially and are numbered accordingly. The CRC
data was not publicly available and so all data/results which contain
data from these summary statistics are not available here (e.g., outcome
data files are not made available for the association between adiposity
and CRC).

## Data

### Adiposity

We used three measures of adiposity: body mass index (BMI), waist hip
ratio (WHR), and WHR adjusted for BMI (WHRadjBMI). Summary statistics
were available from [Pulit et al.,
(2020)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6298238/).
Sex-combined and sex-specific summary statistics were used for all three
measures.

### Proteins

In the main analysis, summary statistics for 4,907 proteins (Somascan)
were available from [Ferkingstad et al.,
(2019)](https://www.nature.com/articles/s41588-021-00978-w) - only
sex-combined data were available. In the additional analysis, summary
statistics for 1,463 proteins (Olink) were available from [Sun et al.,
(2022)](https://www.biorxiv.org/content/10.1101/2022.06.17.496443v1) -
only sex-combined data were available.

### Colorectal cancer

Summary statistics for 16 colorectal cancer outcomes were available from
[Huyghe et al.,
(2019)](https://www.nature.com/articles/s41588-018-0286-6), this
included: overall, overall measured in Europeans only, overall male,
overall female, overall colon, overall proximal, overall distal, overall
rectal, female colon, male colon, female proximal, male proximal, female
distal, male distal, female rectal, and male rectal. In addition,
summary statistics for early onset colorectal cancer were available from
XXX et al., (). These data are nto publicly available and can instead be
requested from the study authors.

## Methods

This study has four main analyses that were performed sequentially to
estimate: (I) the association between adiposity measures and CRC, (II)
the association between adiposity measures and proteins, (III) the
association between adiposity-associated proteins (identified in step
II) and CRC, and (IV) the potential mediating effect of
adiposity-associated proteins in the adiposity CRC association
(identified in step II and III). The reverse association was tested for
each step. Where proteins were used as the exposure, cis and trans
variants were used in the main analysis. For all proteins included in
step IV, the following additional analyses were performed where data
were available: (I) cis-only MR of the protein-CRC pair, (II)
colocalization of the protein-CRC pair, and (III) UVMR using cis and
trans variants from a second protein dataset.

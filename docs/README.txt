Data/Output backup for amr_analysis project.
code available here: https://github.com/henrichung/amr_analysis

Contents:

Data:
+ antibiotic classes.csv - reference table that classifies specific antibiotics into broad types.
+ CLSIbreakpoint_refTAblePalantir.csv - reference table listing the MIC breakpoint values for SIR classification with different ABs and animals.
+ identifier_edited.csv - reference table that associates gene identifiers with specific genes.

+ CARD_lookups - list of genes to manually look up in CARD database.
+ CARD_lookups_post - list of gene names returned by AMRfinder/Abricate with their CARD entry equivalents
+ gene_resistance - list of CARD Antimicrobials and labels whether they are a drug_class or antibiotic.

+ Ecoli_auto_msss.csv - phenotype data for approximately 600 samples from various livestock.
+ ecoli_ppy*_astdata.csv - phenotype data for approxmiately 300 samples from various livestock.

+ reusdaamrdata/ecoli_pp*_abricate|amrfinder.csv - lists gene identifiers found approximately 800 samples from various livestock.


Outputs
+ significant_gene_animal_phenotypes.pdf - plot that shows the AM/animal combinations where there is a signicant difference in the number of AMR genes across SIR phenotypes by animal.
+ significant_gene_animal_phenotypes.csv - table which shows the number of samples where there was a significant difference
in the number of AMR genes across SIR phenotypes by animal.
+ significant_chisquare_gene_animal.csv - results for chisquare analysis where there is a significant interaction between presence/absence of a specific gene and the number of SIR samples per AM.
+ gene_presence_per_animal - table which shows the prevalence of specific gene in animal samples.

+ heatmaps - shows presence/absence of AMR genes clustered in multiple ways.

+ sample_histogram_animal - histogram showing composition of samples by animal.
+ sample_histogram_animal_species - histogram showing composition of samples by animal species name.

Docs
+ *.ppt - powerpoints for project meetings.
+ *.pdf - draft of writeup.
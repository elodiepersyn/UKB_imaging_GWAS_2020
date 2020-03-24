# Genome-wide association study of MRI markers of cerebral small vessel disease in 42,310 participants

This directory contains the code that was used for the analyses in "Genome-wide association study of MRI markers of cerebral small vessel disease in 42,310 participants" (Persyn et al.).

## Content
- The `workflows` directory contains different `snakemake` pipelines:

    - `GWAS_WMH`: meta-analysis GWAS for white matter hyperintensities volume (WMH)
    - `GWAS_FA_MD`: GWAS for fractional anisotropy (FA) and mean diffusivity (MD)
    - `Association_brain_regions`: FA and MD association analyses for a subset of SNPs in the 48 brain regions
    - `GWASsumstats_analyses`: a pipeline for post-GWAS analyses such as:
        - SNP heritability estimation with [`ldsc`](https://github.com/bulik/ldsc)
        - Partionned heritability with [`ldsc`](https://github.com/bulik/ldsc)
        - Genetic correlation estimation with [`ldsc`](https://github.com/bulik/ldsc)
        - Transcriptome-wide association studies and colocalization analyses with [`FUSION`](http://gusevlab.org/projects/fusion/) software.
        - Gene Ontology (GO) gene-set enrichment analysis from GWAS results with [`MAGMA`](https://ctg.cncr.nl/software/magma)
    
    To run these pipelines, it is necessary to modifiy the configuration files in the `config` directory to specify the location of datasets, results and programs. Please read the documentation of [`snakemake`](https://snakemake.readthedocs.io/en/stable/#) tool.


- The `other` directory contains the following scripts :
    
    - `PhenoScanner.R`: annotation of significant SNPs with [`PhenoScanner`](http://www.phenoscanner.medschl.cam.ac.uk/) database to explore the association with other traits
    - `TWAS_GSEA.sh`: Gene Ontology (GO) gene-set enrichment analysis from TWAS results with [`TWAS-GSEA`](https://github.com/opain/TWAS-GSEA)
    - `hyprcoloc_analysis.R`: Colocalization analysis for multiple traits with [`hyprcoloc`](https://github.com/jrs95/hyprcoloc)

## Programs to install

Here is the list of necessary programs to run the analyses:

- [`snakemake`](https://snakemake.readthedocs.io/en/stable/#)
- [`R`](https://www.r-project.org/)
- [`python`](https://www.python.org/) (the version and the dependencies might change in function of other programs)
- [`plink v2.0`](https://www.cog-genomics.org/plink/2.0/)
- [`qctool v2.0`](https://www.well.ox.ac.uk/~gav/qctool/)
- [`plink v1.9`](https://www.cog-genomics.org/plink/1.9/)
- [`metal`](https://genome.sph.umich.edu/wiki/METAL_Documentation)
- [`ldsc`](https://github.com/bulik/ldsc)
- [`FUSION`](http://gusevlab.org/projects/fusion/)
- [`MAGMA`](https://ctg.cncr.nl/software/magma)
- [`TWAS-GSEA`](https://github.com/opain/TWAS-GSEA)

`R` packages to install are :

- `circlize`
- `coloc`
- `data.table`
- `dplyr`
- `FactoMineR`
- `ggplot2`
- [`hyprcoloc`](https://github.com/jrs95/hyprcoloc)
- `mice`
- `ukbtools`
- `yaml`




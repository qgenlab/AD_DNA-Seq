
# Installation
Please refer to `environment.yml` for detail. For your quick reference, you need the packages below to run our scripts
```
  - bcftools=1.20
  - bedtools=2.31.1
  - bioconductor-limma=3.58.1
  - biopython=1.85
  - htslib=1.20
  - ipykernel=6.29.5
  - jupyter=1.0.0
  - make=4.4.1
  - matplotlib=3.9.2
  - matplotlib-base=3.9.2
  - minimap2=2.1.1
  - numpy=1.26.4
  - numpy-base=1.26.4
  - pandas=2.2.2
  - parasail-python=1.3.4
  - polars=1.6.0
  - rpy2=3.5.11
  - scipy=1.13.1
  - pysam=0.22.1
  - python=3.12.5
  - r-base=4.3.3
  - r-statmod=1.5.0
  - seaborn=0.13.2
  - statsmodels=0.14.2
```

You can easily install the dependent packages using `annoconda` or `miniconda`.
```
git clone https://github.com/qgenlab/AD_DNA-Seq
cd AD_DNA-Seq
conda env create -f environment.yml
source activate AD_DNA_Seq
```

# Reproducibility script

All Python and bash scripts are available under the directory of `scripts`. The testing examples were provided here. 



# Help

Please refer to the [AD_DNA-Seq issue page](https://github.com/qgenlab/AD_DNA-Seq/issues) for posting your issues. We will respond your questions quickly. Your comments/suggestions are always appreciated and will benefit other users.

# Citation

Please cite the publication below for script and data usage and/or other purposes

Q. Liu, R. Doshi, H. Derbel, A. Bavero, C. Harris, Z. Zhao, P. Schulz. Multi-faceted genomic and epigenomic profiling in Alzheimer's disease using Oxford Nanopore sequencing. Submitted.


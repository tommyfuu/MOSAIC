# Data sources and licenses

The MOSAIC code is licensed under GPL-3.0 (see `../LICENSE`). The data files in this folder are derived from
public databases and keep the licenses of their sources. Please cite the database and the original studies when
you use them.

## MicrobiomeHD (16S rRNA gene sequencing, genus-level counts)

Folders: `autism_2_microbiomeHD/`, `cdi_3_microbiomeHD/`, `intermediate_*_microbiomeHD/`, `pruned_autism_2_microbiomeHD/`,
`pruned_cdi_3_microbiomeHD/`, and `cleaned_data/autism_2_microbiomeHD/`, `cleaned_data/cdi_3_microbiomeHD/`,
`cleaned_data/cdi_3_microbiomeHD_clean/`.

- License: Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0). These files may be shared
  and adapted with attribution, for non-commercial use only. They are not covered by the GPL-3.0 license of the code.
- Source: Duvallet C, Gibbons SM, Gurry T, Irizarry RA, Alm EJ. Meta-analysis of gut microbiome studies identifies
  disease-specific and shared responses. Nat Commun 8, 1784 (2017). Data: https://zenodo.org/records/1146764
- Original studies:
  - autism: Son JS et al., PLoS One 2015; Kang DW et al., PLoS One 2013
  - CDI: Schubert AM et al., mBio 2014; Vincent C et al., Microbiome 2013; Youngster I et al., Clin Infect Dis 2014

## curatedMetagenomicData (shotgun metagenomics, species-level relative abundances)

Folders: `pruned_ibd_3_CMD/`, `pruned_crc_8_CMD/`, `pruned_crc_11_CMD/`, `T2D_10_CMD/`, and
`cleaned_data/ibd_3_CMD/`, `cleaned_data/crc_8_CMD/`, `cleaned_data/crc_11_CMD/`. `cmd_phenotype_survey.csv` lists
sample counts per disease in the database.

- License: the curatedMetagenomicData Bioconductor package is released under Artistic-2.0.
- Source: Pasolli E et al. Accessible, curated metagenomic data through ExperimentHub. Nat Methods 14, 1023-1024
  (2017). doi:10.18129/B9.bioc.curatedMetagenomicData
- Original studies:
  - IBD: Lloyd-Price J et al., Nature 2019; Schirmer M et al., Nat Microbiol 2018; Li J et al., Nat Biotechnol 2014;
    Nielsen HB et al., Nat Biotechnol 2014
  - CRC (11 studies): Feng Q et al., Nat Commun 2015; Gupta A et al., mSystems 2019; Dhakan DB et al., GigaScience 2019;
    Hannigan GD et al., mBio 2018; Thomas AM et al., Nat Med 2019; Vogtmann E et al., PLoS One 2016; Wirbel J et al.,
    Nat Med 2019; Yachida S et al., Nat Med 2019; Yu J et al., Gut 2017; Zeller G et al., Mol Syst Biol 2014

## Simulation template

`../benchmark/ibd_150.Rdata` is a genus-level subset of the IBD data from the Integrative Human Microbiome Project
(IBDMDB), obtained through the HMP2Data Bioconductor package (Artistic-2.0). Cite: Lloyd-Price J et al., Nature 2019;
Integrative HMP (iHMP) Research Network Consortium, Nature 2019.

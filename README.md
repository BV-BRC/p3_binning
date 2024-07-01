# MetaGenome Binning

## Overview

Metagnome binning is the process of extracting individual genomes from
a metagenomic sample.  If the sample is in the form of reads, they are
assembled into contigs first.  Binning contigs involves separating them
into sets (called _bins_) that represent contigs believed to be from
a single genome.

This is a supervised binning process:  it will only return bins that match
genomes in BV-BRC known to be of high quality.  You can specify binning of
bacterial/archael genomes, virus genomes, or both.  All bins will be processed
by the [Genome Annotation Service](https://bv-brc.org/app/Annotation) to
produce annotations and determine the quality of the bin.  The resulting
genomes will be indexed in the main BV-BRC database as private data and available in 
searches as well as for use in other BV-BRC services.

## About this module

This module is a component of the BV-BRC build system. It is designed to fit into the
`dev_container` infrastructure which manages development and production deployment of
the components of the BV-BRC. More documentation is available [here](https://github.com/BV-BRC/dev_container/tree/master/README.md).

The application service specification is given in [MetaGenomeBinning](app_specs/MetagenomeBinning.md)the BV-BRC Comprehensive Genome Annotation service.

The code for the metagenome service itself, along with the code for genome evaluation, 
is found [here](https://github.com/BV-BRC/p3_code).


## References

Aziz, R. K. et al. The RAST Server: rapid annotations using subsystems technology. BMC genomics 9, 75 (2008).

VIGOR4, https://github.com/JCVenterInstitute/VIGOR4.

Parrello B, Butler R, Chlenski P, Olson R, Overbeek J, Pusch G, Vonstein V, Overbeek R. (2019) A machine learning-based service for estimating quality of genomes using PATRIC. BMC Bioinformatics 20, 486 (2019).

Parrello B, Butler R, Chlenski P, Pusch G, Overbeek R. (2021) Supervised extraction of near-complete genomes from metagenomic samples: A new service in PATRIC. PLoS ONE 16(4):e0250092. April 2021
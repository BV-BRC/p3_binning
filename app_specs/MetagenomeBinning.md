
# Application specification: MetagenomeBinning

This is the application specification for service with identifier MetagenomeBinning.

The backend script implementing the application is [App-MetagenomeBinning.pl](../service-scripts/App-MetagenomeBinning.pl).

The raw JSON file for this specification is [MetagenomeBinning.json](MetagenomeBinning.json).

This service performs the following task:   Assemble, bin, and annotate metagenomic sample data

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | :------: | ------------ |
| paired_end_libs |  | group  |  |  |
| single_end_libs |  | group  |  |  |
| srr_ids | SRR ID | string  |  |  |
| contigs | Contig file | WS: Contigs  |  |  |
| genome_group | Output Genome Group | string  |  |  |
| skip_indexing | Don't index bins | bool  |  | 0 |
| recipe | Annotation recipe | string  |  |  |
| viral_recipe | Viral nnotation recipe | string  |  |  |
| output_path | Output Folder | folder  | :heavy_check_mark: |  |
| output_file | File Basename | wsid  | :heavy_check_mark: |  |
| force_local_assembly | Force local assembly | bool  | :heavy_check_mark: | 0 |
| force_inline_annotation | Force inline annotation | bool  |  | 1 |
| perform_bacterial_binning | Perform bacterial binning | bool  |  | 1 |
| perform_viral_binning | Perform viral binning | bool  |  | 0 |
| perform_viral_annotation | Perform viral annotation | bool  |  | 0 |
| perform_bacterial_annotation | Perform bacterial annotation | bool  |  | 1 |
| assembler | Assembler to use | string  |  |  |
| danglen | Dangling element length | string  |  | 50 |
| min_contig_len | Minimal output contig length | int  |  | 400 |
| min_contig_cov | Minimal output contig coverage | float  |  | 4 |


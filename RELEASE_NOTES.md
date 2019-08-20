### Version 2.1.2
- annotate_metagenome added in order to creat KBaseMetagenome.AnnotatedMetagenomeAssembly output objects.
  the function removes a limit on the size of the incoming fasta/assembly.

### Version 1.0.4
- Features IDs in the annotated genomes are now unique IDs generated
  by Prokka, not gene names (the latter are now stored as feature aliases).

### Version 1.0.3
- Check didn't work on for older contigs that were missing the metadata field

### Version 1.0.2
- Added check and error message for large contig inputs


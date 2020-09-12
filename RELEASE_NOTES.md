### Version 3.2
- Per Chris Henry's request deactiving SEED Subsystem Ontologies
- The EC numbers are no longer made into aliases and instead are made into ontologies
- The EC is currently going from a file in the Repo, down the line we need to move that to an ontology object
- TODO Chris would like COGs captured as well as ontologgies from PROKKA output.
- Note this needs to be done carefully because another ontology event needs to be made and kept track of appropriately

### Version 3.1
- Updated the backend Prokka to 1.14.5

### Version 2.1.5
- Adding in RE taxonomy dropdown for scientific_name field.
- Removing output viewer for AnnotatedMetagenomeAssembly

### Version 2.1.4
- Bug fix for setting ontologies_present to the correct type in a certain execution path

### Version 2.1.3
- Fixing output methods for annotate_metagenome to no longer create an extra Assembly object.
- Bug fixes, all tests passing.

### Version 2.1.2
- annotate_metagenome added in order to create KBaseMetagenome.AnnotatedMetagenomeAssembly output objects.
  the function removes a limit on the size of the incoming fasta/assembly.

### Version 1.0.4
- Features IDs in the annotated genomes are now unique IDs generated
  by Prokka, not gene names (the latter are now stored as feature aliases).

### Version 1.0.3
- Check didn't work on for older contigs that were missing the metadata field

### Version 1.0.2
- Added check and error message for large contig inputs


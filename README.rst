[![Build Status](https://travis-ci.org/rsutormin/ProkkaAnnotation.svg?branch=master)](https://travis-ci.org/rsutormin/ProkkaAnnotation)


Prokka Annotate Assembly
^^^^^^^^^^^^^^^^^^^^^^^^
1) Retrieve the seed system ontology from "KBaseOntology/seed_subsystem_ontology" and create a table of SSO Numbers with EC numbers that map to them.

2) Inspect assembly to see if it has too many contigs (>30,000). This might mean it is not a metagenome or is non prokaryotic.

3) Retrieve the FASTA from the assembly_ref using AssemblyUtils.getAssembly_as_fasta, and rename the contigs, but keep a record of the old names

5) Run Prokka with the modified, renamed fasta file

6) Retrieve the results, check for FAA, FFN, and GFF files in output

7) Create an EVIDENCE EVENT for adding to the SSO_REF.

8) Parse the output GFF, but map the contigs back to the old IDs

9) For each record in the GFF get a list of features, and for each feature, get the following fields

feature = {"id": fid, "location": location, "type": "gene",
           "aliases": aliases, "md5": md5, "dna_sequence": dna,
           "dna_sequence_length": len(dna)}


10) Count the number of non hypothetical proteins and number of genes with SSO Refs

11) If the discovered EC number in the list of known ecs with an SSO,, create an SSO Ref and add to Ontology_terms for the feature

12) For feature, if a CDS or MRNA is discovered, add that to the list of CDS and MRNAS

13) Create a genome with the features, cds, mrnas and various parameters input during runtime and then save the genome and create a report

Prokka Re-Annotate Genome
^^^^^^^^^^^^^^^^^^^^^^^^^
1) Retrieve the seed system ontology, as in step 1 above

2) Using the genome api, download a non downgraded genome

3) For each feature, extract the fasta dna sequence and write to a output features.fasta file and run prokka on the features.fasta file

4) Search the output GFF file for its features. For each feature, extract the product and eC_Number

5) Create evidence field for any SSO annotations. If the ec number is found in the seed system ontology, create an SSO Event. Save the following from the gff: 'id','function->product' and 'ontology_terms->eC_number'

6) Annotate the existing genome features with the new function. The function is overwritten with the new function. The functions array is overwritten to be an array with one item, the function.

7) If it is a new genome, append a SSO number to the ontology event, as well as an ontology event index, in case of multiple annotations of the same SSO Number

8) If it is an old genome, append the actual SSO event in the style as in the Prokka assembly directly into the genome. The evidence will automatically be set correctly.

13) Save the genome and create a report, updating only the output_genome name and the features. No other fields will be updated.
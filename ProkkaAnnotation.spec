/*
A KBase module: ProkkaAnnotation
*/

module ProkkaAnnotation {

    /*
        Reference to an Assembly object in the workspace
        @id ws KBaseGenomeAnnotations.Assembly
    */
    typedef string assembly_ref;

    /*
        Reference to an Genome object in the workspace
        @id ws KBaseGenomes.Genome
    */
    typedef string genome_ref;

    typedef structure {
        assembly_ref assembly_ref;
        string output_workspace;
        string output_genome_name;
    } AnnotateContigsParams;

    typedef structure {
        genome_ref output_genome_ref;
    } AnnotateContigsOutput;

    funcdef annotate_contigs(AnnotateContigsParams params) returns (AnnotateContigsOutput) authentication required;
};

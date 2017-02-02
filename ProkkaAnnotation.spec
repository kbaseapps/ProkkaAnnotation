/*
A KBase module: ProkkaAnnotation
*/

module ProkkaAnnotation {

    /* A boolean. 0 = false, anything else = true. */
    typedef int boolean;

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

    /*
        Required parameters:
            assembly_ref - reference to Assembly object,
            output_workspace - output workspace name,
            output_genome_name - output object name,
        Optional parameters (correspond to PROKKA command line arguments):
          --kingdom [X]     Annotation mode: Archaea|Bacteria|Mitochondria|Viruses (default 'Bacteria')
          --genus [X]       Genus name (triggers to use --usegenus)
          --gcode [N]       Genetic code / Translation table (set if --kingdom is set) (default '0')
          --gram [X]        Gram: -/neg +/pos (default '')
          --metagenome      Improve gene predictions for highly fragmented genomes (default OFF)
          --rawproduct      Do not clean up /product annotation (default OFF)
          --fast            Fast mode - skip CDS /product searching (default OFF)
          --mincontiglen [N] Minimum contig size [NCBI needs 200] (default '1')
          --evalue [n.n]    Similarity e-value cut-off (default '1e-06')
          --rfam            Enable searching for ncRNAs with Infernal+Rfam (SLOW!) (default '0')
          --norrna          Don't run rRNA search (default OFF)
          --notrna          Don't run tRNA search (default OFF)
    */
    typedef structure {
        assembly_ref assembly_ref;
        string output_workspace;
        string output_genome_name;
        string kingdom;
        string genus;
        int gcode;
        string gram;
        boolean metagenome;
        boolean rawproduct;
        boolean fast;
        int mincontiglen;
        string evalue;
        boolean rfam;
        boolean norrna;
        boolean notrna;
    } AnnotateContigsParams;

    typedef structure {
        genome_ref output_genome_ref;
        string report_name;
        string report_ref;
    } AnnotateContigsOutput;

    funcdef annotate_contigs(AnnotateContigsParams params) returns (AnnotateContigsOutput) authentication required;
};

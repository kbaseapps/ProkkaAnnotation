/*
A KBase module: ProkkaAnnotation
*/

module ProkkaAnnotation {
   /* A boolean. 0 = false, anything else = true. */
    typedef int boolean;

    /*
        Reference to an Assembly or Genome object in the workspace
        @id ws KBaseGenomeAnnotations.Assembly
        @id ws KBaseGenomes.Genome
    */
    typedef string data_obj_ref;

    /*
        Reference to an Genome object in the workspace
        @id ws KBaseGenomes.Genome
    */
    typedef string genome_ref;

    /*
        Required parameters:
            object_ref - reference to Assembly or Genome object,
            output_workspace - output workspace name,
            output_genome_name - output object name,
        Optional parameters (correspond to PROKKA command line arguments):
          --scientific_name Genome scientific name (default 'Unknown')
          --kingdom [X]     Annotation mode: Archaea|Bacteria|Mitochondria|Viruses (default 'Bacteria')
          --genus [X]       Genus name (triggers to use --usegenus)
          --gcode [N]       Genetic code / Translation table (set if --kingdom is set) (default '11')
          --metagenome      Improve gene predictions for highly fragmented genomes (default OFF)
          --rawproduct      Do not clean up /product annotation (default OFF)
          --fast            Fast mode - skip CDS /product searching (default OFF)
          --mincontiglen [N] Minimum contig size [NCBI needs 200] (default '1')
          --evalue [n.n]    Similarity e-value cut-off (default '1e-06')
          --rfam            Enable searching for ncRNAs with Infernal+Rfam (SLOW!) (default OFF)
          --norrna          Don't run rRNA search (default OFF)
          --notrna          Don't run tRNA search (default OFF)
    */
    typedef structure {
        data_obj_ref object_ref;
        string output_workspace;
        string output_genome_name;
        string scientific_name;
        string kingdom;
        string genus;
        int gcode;
        boolean metagenome;
        boolean rawproduct;
        boolean fast;
        int mincontiglen;
        string evalue;
        boolean rfam;
        boolean norrna;
        boolean notrna;
    } AnnotateParams;

    typedef structure {
        genome_ref output_genome_ref;
        string report_name;
        string report_ref;
    } AnnotateOutput;

    funcdef annotate(AnnotateParams params) returns (AnnotateOutput) authentication required;
};

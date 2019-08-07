
# mgv-data 

Scripts for building the server-side data to serve the Multiple Genome Viewer (MGV).

The data backend for MGV is organized by genome, and comprises gene models and genome assemblies.
The default top-level script (build.sh) builds a data set for MGV comprising the 19 sequenced and 
annotated inbred strains.
Most of the data comes from Ensembl, although specific pieces and certain transformations are MGI specific.
The scripts are designed to be reconfigurable with little effort for building
from other genomes available at Ensembl.

## Mouse-specific bits
For each of the mouse strains, the build script retrieves the genome assembly (.fa) and genome annotation 
(.gff3) files from Ensembl.
One exception is the genome annotation for the C57BL/6J strain, which comes from MGI. The MGI file
combines gene model predictions from multiple sources, including Ensembl.
In addition, all genes from all strains are tagged with their equivalent ("canonical") MGI id.

## Files

### ./bin:
* build.sh Top level build script
* genomes.tsv List of genomes the drives a build. Read by build.sh. Table has 3 columns: genomePath, genomeName,
and taxonId. genomePath is the name of the organism as it appears in paths at Ensembl, e.g., "mus_musculus_aj".
genomeName is a label to use for the genome, e.g., "A/J". TaxonId is the NCBI taxon id for the organism, e.g., "10090"
* getGenome.sh Gets genome data from Ensembl and imports into MGV data backend.
* importGff3.py Imports a genome annotation (GFF3) file for one genome, splits it up into chunks, and stores it 
in the right directory structure.
* importFasta.py Imports a genome assembly (Fasta) file for one genome, and writes it into the format requires for MGV.
* sequenceHound.py A CGI that extracts specified sequences from the assembly files. Called by the SequenceCart download
function in MGV.

### ./downloads
Where compressed GFF3 and Fasta files go when downloaded from Ensembl/MGI.

### ./output
Where the generated output goes, i.e., the directory that will be served to MGV. Each genome to be served is a separate subdirectory under this one. 
* ./output/index.json A json-formatted list of the subdirectories (genomes) to be served. MGV does not attempt to determine automatically what genomes are available. Instead, it reads this file to find out. The contents is a single list of subdirectory names followed by "/", e.g. `[ "mus_musculus_aj/", "mus_musculus_dba2j/", ...]`
* ./output/fetch.cgi A CGI wrapper script that invoked `./bin/sequenceHound.py`
* ./output/mus_musculus_aj/ (et. al.) The data for each genome is stored in its own subdirectory
** ./output/mus_musculus_aj/index.json A json formatted object describing the genome
** ./output/mus_musculus_aj/genes Contains a single file, named '0' which is a GFF3 file containing just the top level gene features for the genome.
** ./output/mus_musculus_aj/transcripts Contains all the transcripts, exons, and CDSs. These are divided into subdirectories by chromosome, and the data for each chromosome is divided into 4MB chunks, each named simply by its chunk number. So for example, the file ./output/mus_musculus_aj/transcripts/12/8 refers to the 8th chunk (which covers 32-36 Mbp) of chromosome 12 of the A/J strain.
** ./output/mus_musculus_aj/sequences Contains the genome assembly, one file per chromosome, named by that chromsome. E.g., ./output/mus_musculus_aj/sequences/10 is the file containing the sequence for chromosome 10 of A/J. The sequence files are not
Fasta, but simply the base sequence. (No header line, no line breaks.)

Gene models are stored as modified GFF3 as follows:
* all genes (top level features) for the genome are in a single file
* transcripts, exons, and CDSs are stored in specialized GFF3 files:
 * organized by chromosome and divided into 4MB chunks
 * 




# mgv-data 

Scripts for building the server-side data directories to serve the Multiple Genome Viewer (MGV).

The data backend for MGV is organized by genome, and comprises gene models and genome assemblies.
The default top-level script (build.sh) builds a data set for MGV comprising the 19 sequenced and 
annotated inbred strains.
Most of the data comes from Ensembl, although specific pieces and certain transformations are MGI specific.
The scripts are designed to be reconfigurable with little effort for building
from other genomes available at Ensembl.

## MGI-specific bits

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

### ./downloads
Where compressed GFF3 and Fasta files go when downloaded from Ensembl/MGI.

### ./output

Gene models are stored as modified GFF3 as follows:
* all genes (top level features) for the genome are in a single file
* transcripts, exons, and CDSs are stored in specialized GFF3 files:
 * organized by chromosome and divided into 4MB chunks
 * 

built from plain files and one CGI script.


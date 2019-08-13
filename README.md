
# mgv-data 

Scripts for building and serving data for the Multiple Genome Viewer (MGV).

## Building the default data set
This repo comes preconfigured to build the data set being served by MGV at MGI,
which consists of the genome assemblies and gene model annotations for 19 inbred strains of mice.
Most of the data comes from Ensembl, except for C57BL/6J gene models and MGI canonical ids and nomenclature,
which come from MGI. To build this data set:
1. Install the repo and cd to the bin directory. 
```bash
git clone git@github.com:JoelRichardson/mgv-data.git
cd mgv-data/bin
```
2. Edit config.sh for your environment. You'll probably want to set data directories for where files are downloaded (DDIR), built into (ODIR) and served from (WDIR).
3. Ensure things are executable:
 ```bash
chmod u+x build.sh getGenome deploy.sh
```
4. Run the build: 
 ```bash
./build.sh
```

## Serving data to MGV

The data that serves MGV is a directory hierarchy of static files, organized by genome.
The data comprise genome features, stored as (modified) GFF3 and served as static files, and genome assemblies, stored as plain strings and served by a CGI script.

1. Run the deployment script. The deployment directory is configured in config.sh. 
 ```bash
./deploy.sh
```
Rsync is used to copy new/changed files from the build area to the deployment area. 
You can save space and time by building directly into the deployment area (set both ODIR and WDIR to the same value).

The CGI is a Python script, fetch.py, which is invoked by a shell wrapper, fetch.cgi. The deploy script copies both pieces to the deployment directory and makes the wrapper executable. You may have to adjust this step, depending on your local environment. 

## Customizing a build

If you're using data from Ensembl:
1. Edit build.sh to call getGenome.sh for the organisms you want. Optionally change the default version number (in config.sh).
2. If you need to do custom translations on the GFF3, supply/write the appropriate Python modules, and specify them to getGenome.sh (see -m command line option).

The existing build script (build.sh) and custom translations (pg_MGI.py, pg_Ensembl.py, and pg_tagEnsemblWithMgi.py) should provide sufficient examples to go by.

Using data NOT from Ensembl:

The data files need not come from Ensembl, provided they are in the same format as those provided by Ensembl.
1. Use the --gff-url argument to directly specify the location of the compressed GFF3 (.gff3.gz) file containing the genome annotations.
2. Use the --fasta-url argument to directly specify the location of the compressed Fasta (.fa.gz) file containing the genome assembly.

## Directory and file structure

### ./bin:
* build.sh Top level build script
* getGenome.sh Gets genome data from Ensembl and imports into MGV data backend.
* importGff3.py Imports a genome annotation (GFF3) file for one genome, splits it up into chunks, and stores it 
in the right directory structure.
* importFasta.py Imports a genome assembly (Fasta) file for one genome, and writes it into the format requires for MGV.
* fetch.py A CGI that extracts specified sequences from the assembly files. Called by the SequenceCart download
function in MGV.
* fetch.cgi A shell wrapper that invokes fetch.py
* pg_MGI.py Custom translator appied to the MGI.gff3 file. Sets the cID attribute.
* pg_ensembl.py Custom translator for Ensembl. Sets column 3 to "protein_coding_gene" when col 3 says "gene" and col 9 biotype says "protein_coding".
* pg_tagEnsemblWithMgi.py Custom translator applied Ensembl gene records. Finds corresponding MGI gene and tags record with the MGI id and symbol.

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



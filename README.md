
# mgv-data 

Scripts for building and serving data for the Multiple Genome Viewer (MGV).

## Building the default data set
This repo comes preconfigured to build the data set being served by MGV at MGI,
which consists of the genome assemblies and gene model annotations for 19 inbred strains of mice.
Most of the data comes from Ensembl, except for C57BL/6J gene models and MGI canonical ids and nomenclature,
which come from MGI. To build this data set:
1. Install the repo and cd to the bin directory. 
```bash
git clone git@github.com:mgijax/mgv-data.git
(or git clone https://github.com/mgijax/mgv-data.git)
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

Note that these file are NOT compressed.
Instead, the web server should be configured to use http compression.
For example, to configure Apache, you could add a .htaccess file with the following:
```
AddType application/json .json
AddType text/tab-separated-values .gff3
AddOutputFilterByType DEFLATE application/json text/tab-separated-values
```

The CGI is a Python script, fetch.py, which is invoked by a shell wrapper, fetch.cgi. The deploy script copies both pieces to the deployment directory and makes the wrapper executable. You may have to adjust this step, depending on your local environment. 

## Customizing a build

If you're using data from Ensembl:
1. Edit build.sh to call getGenome.sh for the organisms you want. Optionally change the default version number (in config.sh).
2. If you need to do custom translations on the GFF3, supply/write the appropriate Python modules, and specify them to getGenome.sh (see -m command line option).
3. One primary purpose of a custom translation is setting the cID ("canonical" or "class" ID) attribute on the genes. Generally, this is the official/MOD ID of the gene (eg, MGI id for mouse genes, HGNC id for human, FlyBase for drosophila, etc. The cID is the determinant of equivalence between genomes in the same species. It is also the basis for homology relationships across species. For a custom build, you'll need to figure out how to tag gene features with appropriate cID values, and write appropriate translators to do it.

The existing build script (build.sh) and custom translations should provide sufficient examples to go by.

Using data NOT from Ensembl:

The data files need not come from Ensembl, provided they are in the same format as those provided by Ensembl.
1. Use the --gff-url argument to directly specify the location of the compressed GFF3 (.gff3.gz) file containing the genome annotations.
2. Use the --fasta-url argument to directly specify the location of the compressed Fasta (.fa.gz) file containing the genome assembly.

## Scripts

* build.sh Top level build script
* getGenome.sh Gets genome data from Ensembl and imports into MGV data backend.
* importGff3.py Imports a genome annotation (GFF3) file for one genome, splits it up into chunks, and stores it 
in the right directory structure. 
* getHomologies.sh Downloads homology data from the Alliance of Genome Resources. This consists of gene-to-gene orthology assertions.
* importFasta.py Imports a genome assembly (Fasta) file for one genome, and writes it into the format requires for MGV.
* fetch.py A CGI that extracts specified sequences from the assembly files. Called by the SequenceCart download
function in MGV.
* fetch.cgi A shell wrapper that invokes fetch.py
* MGI.py Custom translator appied to the MGI.gff3 file. Sets the cID attribute.
* ensembl.py Custom translator for Ensembl. Sets column 3 to "protein_coding_gene" when col 3 says "gene" and col 9 biotype says "protein_coding".
* tagEnsemblFish.py Tags D.rerio genes with their ZFIN ids.
* tagEnsemblFly.py Tags  D.melonogater with their FlyBase ids.
* tagEnsemblHuman.py Tags H.sapiens genes with their HGNC ids.
* tagEnsemblWorm.py Tags C.elegans genes with their WormBase ids.
* tagEnsemblYeast.py Tags S.cerevisiae genes with their SGD ids.
* tagEnsemblRat.py Tags R.norvegicus genes with their RGD ids.

## Directory and file structure

Where the generated output goes, i.e., the directory that will be served to MGV. Each genome to be served is a separate subdirectory under this one. 
* ./index.json A json-formatted list of the subdirectories (genomes) to be served. MGV does not attempt to determine automatically what genomes are available. Instead, it reads this file to find out. The contents is a single list of subdirectory names followed by "/", e.g. `[ "mus_musculus_aj/", "mus_musculus_dba2j/", ...]`
* ./fetch.cgi A CGI wrapper script that invokes `./fetch.py`
* ./fetch.py The actual CGI. Its job is to extract subsequences from the genome assemblies, possibly doing reverse complementation or protein translation, and return the results as Fasta.
* ./mus_musculus_aj/ (et. al.) The data for each genome is stored in its own subdirectory
** ./mus_musculus_aj/index.json A json formatted object describing the genome
** ./mus_musculus_aj/genes/ Contains a single file, named '0.gff3' which is a GFF3 file containing just the top level gene features for the genome.
** ./mus_musculus_aj/transcripts/ Contains all the transcripts, exons, and CDSs. These are divided into subdirectories by chromosome, and the data for each chromosome is divided into 4MB chunks, each named simply by its chunk number. So for example, the file ./mus_musculus_aj/transcripts/12/8.gff3 refers to the 8th chunk (which covers 32-36 Mbp) of chromosome 12 of the A/J strain.
** ./mus_musculus_aj/sequences/ Contains the genome assembly, one file per chromosome, named by that chromsome. E.g., ./mus_musculus_aj/sequences/10.txt is the file containing the sequence for chromosome 10 of A/J. The sequence files are not
Fasta, but simply the linearized sequence. (No header line, no line breaks.)

Gene models are stored as modified GFF3 as follows:
* all genes (top level features) for the genome are in a single file
* transcripts, exons, and CDSs are stored in specialized GFF3 files:
 * organized by chromosome and divided into 4MB chunks
 * 



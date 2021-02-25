
# mgv_data 

Scripts for building and serving data for the [Multiple Genome Viewer (MGV)](https://github.com/mgijax/mgv).

## Licensing
See the [LICENSE](LICENSE.txt) file for license rights and limitations (MIT expat).

## Overview
The Multiple Genome Viewer (MGV) displays gene models, genome sequences, and homology relationships for one or more configured genomes.
The purpose of this repo, mgv_data, is to build the server-side data repository which comprises
GFF3, plain text, and JSON files, organized in a directory hierarchy by genome and type.
A "build" generally involves (1) downloading data files from various providers, (2) transforming data to the
formats and file organization used internally by MGV, and (3) deploying data to a web-accessible location.
At runtime, most requests by MGV are for static files in this hierarchy.
There is one CGI script that serves requests for sequences.

## Building the default data set
This repository comes preconfigured to build the data set being served at MGI (www.informatics.jax.org/mgv),
which consists of the genome assemblies and gene model annotations for 19 inbred strains of mice,
the genome assemblies and gene model annotations for human, rat, and several other species, and
homology relationships among genomes.
The data come from Ensembl, NCBI, the Alliance of Genome Resources, and MGI.
To build this data set:
1. Install the repo and cd to the main directory. 
```bash
$ git clone git@github.com:mgijax/mgv_data.git
(or git clone https://github.com/mgijax/mgv_data.git)
$ cd mgv_data
```
2. First verify that you can run the wrapper script and the underlying Python script:
```bash
$ ./build -h
```
You should see a help message showing all the paramaters you can specify.
If you get an error from Python, make sure the PYTHON variable defined in wrapperConfig.sh points to Python 3 for your environment.

3. Test a small piece of the build by downloading the genome annotation (GFF3) data for the A/J strain:
```bash
$ ./build -g mus_musculus_aj -t models -p download
```
The results should be saved under ./mgv_data_files:
```bash
$ ls -R mgv_data_files/
mgv_data_files/:
downloads

mgv_data_files/downloads:
mus_musculus_aj

mgv_data_files/downloads/mus_musculus_aj:
Mus_musculus_aj.A_J_v1.101.gff3.gz
```
4. At this point, you should be able to run the entire build. However, you may want to change some configuration settings,
(e.g., output directories) and now would be a good time to do that. Use your favarite editor to edit wrapperConfig.sh:
```bash
$ vim wrapperConfig.sh
```
See file comments for details.

5. A build comprises a series of phases (download/import/deploy) that are run for a set of data types (assemblies, models, and orthologies) for a list of genomes (multiple mouse strains, plus human, rat, etc). As indicated above, you can selectively run specific phases, types, etc.
E.g., to completely refresh the gene model data (run all phases) for all the mouse strains:
```bash
$ ./build -g "mus_.*" -t models
```
To deploy the assemblies of all genomes:
```bash
$ ./build -p deploy -t assembly
```
Note that the download and import phases must already have been run for this command to work. (The build script is not that smart.)

To run a full build, invoke the script with no arguments
```bash
$ ./build
```

2. The top level directory contains:
  * build: a wrapper script for running builds.
  * wrapperConfig.sh: used by the build script, mainly defines various directories 
  By default, everything will be output to a subdirectory of there mgv_data is installed.
  This is fine for testing purposes, but you'll probably want to change that for a production setting.
  * buildConfig.json: defines the actual build details, data sources, etc.
  To customize a build, you will need to edit this file.
  Alternatively, you can create another config file and change the BCONFIG variable in wrapperConfig.sh to point to it.
  See below for details on configuing builds.


## Serving data to MGV

The data that serves MGV is a directory hierarchy of static files, organized by genome.
The data comprise genome features, stored as GFF3 and served as static files and genome assemblies, 
stored as plain strings and served by a CGI script.

Note that these file are NOT compressed.

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



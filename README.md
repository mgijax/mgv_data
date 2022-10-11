
# mgv_data 

Scripts for building and serving data for the [Multiple Genome Viewer (MGV)](https://github.com/mgijax/mgv).

## Licensing
See the [LICENSE](LICENSE.txt) file for license rights and limitations (MIT expat).

## Overview
The Multiple Genome Viewer (MGV) displays gene models, genome sequences, and homology relationships for one or 
more configured genomes.
The purpose of this repo, mgv_data, is to build the server-side data repository which comprises
indexed GFF3, FASTA, and JSON files, organized in a directory hierarchy by genome.
A "build" generally involves (1) downloading data files from various providers, (2) transforming data to the
formats and file organization used internally by MGV, and (3) deploying data to a web-accessible location.
At runtime, all data requests go through a fetch CGI.


## Installation (SECTION IN PROGRESS)

Prequisites: Samtools
The back end relies on three tools in the Samtools suite: tabix, for indexing and retrieval of GFF data; 
faidx, for indexing/rerieval of FASTA sequences; and the compression tool, bgzip.
Samtools is available from: http://www.htslib.org/

1. Clone this repo
```bash
$ git clone git@github.com:mgijax/mgv_data.git
(or git clone https://github.com/mgijax/mgv_data.git)
$ cd mgv_data
```
2. Compile the build script. 
```
$ cd mgv_data
$ bin/compile config.yaml build.sh
```
3. Run the build
```
$ ./build.sh
```




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

To run a build with a custom config file:
```bash
$ ./build -b myConfig.yaml
```

2. The top level directory contains:
  * build: a wrapper script for running builds.
  * wrapperConfig.sh: used by the build script, mainly defines various directories 
  By default, everything will be output to a subdirectory of there mgv_data is installed.
  This is fine for testing purposes, but you'll probably want to change that for a production setting.
  * buildConfig.yaml: defines the actual build details, data sources, etc.
  To customize a build, you will need to edit this file or create your own.
  To run the build using your config file, use the -b option.
  Alternatively, you can change the BCONFIG variable in wrapperConfig.sh.
  See below for details on configuing builds.

## Serving data to MGV

The data that serves MGV is a directory hierarchy of static files, organized by genome.
The data comprise genome features, stored as GFF3 and served as static files and genome assemblies, 
stored as plain strings and served by a CGI script.

Note that these files are stored uncompressed, and we use server compression instead.
The following Apache .htaccess file is included in the www directory and is copied to the data root directory during deployment.
```
AddType application/json .json
AddType text/tab-separated-values .gff3
AddOutputFilterByType DEFLATE application/json text/tab-separated-values
```
The CGI is a Python 3 script, fetch.py, which is invoked by a shell wrapper, fetch.cgi. Deployment copies both pieces to the deployment directory and makes the wrapper executable. You may have to adjust this step, depending on your local environment. 

## Customizing a build

If you're using data from Ensembl:

Using data NOT from Ensembl:

# Internal data format
## File structure.
* Genomes are served to MGV from a static file directory.
  * TO specify this directory, set WDIR in wrapperConfig.sh
  * Gene models are transferred by direct file requests.
  * Sequences are served by a CGI that reads genome sequence files.
* Each genome is a subdirectory of a common root data directory.
  * The root directory contains an info file (index.json) naming the available genomes.
  * Each genome directory also contains an index.json, providing metadata on that genome.
  * Gene models and genome assemblies hang off the genomes' root directories.

## General flow of requests from MGV to the back end (helps explain the file structures):
* When viewer loads:
  * reads the index.json file at the configured location (see mgv/public/runtimeConfig.json)
  * reads the index.json file for each genome configured in the root index.json
* When a genome is first selected for display, reads the genes gff file (see below) for that genome
* As user navigates to different regions, reads transcript chunk files (see below)
* When the user zooms in far enoughs and when the user downloads a sequences, invokes the fetch CGI script (see below)

## Root info file (index.json):
* simply a list of subdirectory names (genomes) to serve
  * all must end with "/"
* subdirectories not listed here are not visible to MGV
* example: ["mus_caroli/", "mus_musculus_129s1svimj/" ]

## Genome info file (<genome>/index.json):
* Metadata for this genome contained in a json object. Fields:
  * name - the name of the genome (displayed to the user)
  * shortname - (deprecated) allows a shorter version of the name in URLs 
  * pathname - the name of the subdirectory containing this genome
  * timestamp - when this genome's data was last updated. The browser caches data on a per genome basis, and uses the timestamp to know when to flush the cache.
  * chromosomes - list of chromosomes for this genome.
    * each chromosome has a name and a length
    * the name is the same as that used for naming chromosome files and directories
    * the order of the list is the order the browser presents them
  * tracks - configures the data types and access for this genome. (At this point, the tracks concept is more vision that reality, and the configuration reflects this!) Each track descriptor has a name and a type, and may have additional parameters
    * name: one of: genes, transcripts, sequences
    * type: one of: PlainSequence, ChunkedGff3
    * chunkSize (if type = ChunkedGFF3) - size of one chunk. 
  * linkouts - used by viewer to display links in a feature's popup. Each linkout spec has:
    * attr - which attribute to use for linking
    * url - link url prefix. The linking attr is appended to form the complete URL
    * text - the link text to display
  * metadata - the metadata values to be displayed for a genome (Nothing listed above is displayed automatically. If you want it visible, include it here.)
* Example:
```javascript
        {
          "name": "C57BL/6J (GRCm38)",
          "shortname": null,
          "pathname": "mus_musculus_grcm38",
          "timestamp": "Mon Apr  5 14:54:13 2021",
          "chromosomes": [
            {
              "name": "1",
              "length": 195471971
            },

            {
              "name": "Y",
              "length": 91744698
            } 
          ],
          "tracks": [
            { 
              "name": "genes",
              "type": "ChunkedGff3",
              "chunkSize": 0
            },
            { 
              "name": "transcripts",
              "type": "ChunkedGff3",
              "chunkSize": "4000000"
            },
            {
              "name": "sequences",
              "type": "PlainSequence"
            } 
          ],
          "linkouts": [
            { 
              "text": "MGI",
              "url": "http://www.informatics.jax.org/accession/",
              "attr": "cID"
            }
          ],
          "metadata": {
            "name": "C57BL/6J (GRCm38)",
            "pathname": "mus_musculus_grcm38",
            "taxonid": "10090",
            "assemblyBuild": "GRCm38",
            "annotationSource": "mgi",
            "annotationRelease": "4 Jan 2021"
          }
        }
```

## Genome assemblies.
* In a genome's assembly directory, there is one file per chromosome, named <chr>.txt, e.g., 11.txt for chromosome 11.
* The file contains the forward stand of that chromosome as a single string with no header and no line breaks. 
* The "PlainSequence" type in the track descriptor referes to this format

## Gene models. Gene models are partitioned in the internal format. 
* Top-level features are split off to a separate file
* Transcripts features hierarchies are compressed and partitioned into fixed-size chunks
* The fame format ("ChunkedGFF3") is used for both genes and transcripts. For the gene track, chunkSize=0, which means everything going in a single file; for transcripts, the chunkSize is 4Mb. For genomes with higher or lower feature densities, may want to adjust this number.

### Genes (et al). 
* All the top level features for a genome (and only the top level features) are stored in the file <genome>/models/genes/0.gff3
* The file is in GFF3 format.
* It is sorted by chromosome and start position.
* The file should not contain very large features (like chromosomes). 
* The ID of each feature must be universally unique (e.g., Pax6 in genome A must have a different ID than Pax6 in genome B).
* The canonical id for each feature (for a mouse gene, this would be the MGI id) should be attached via a "cID" attribute in column 9.
* The only column 9 attributes used by the viewer are ID, cID, and name. Other attributes are fine, but currently ignored.
* Example:
```
1       MGI     protein_coding_gene     3269956 3741733 .       -       .       ID=MGI_C57BL6J_3528744_GRCm39;Name=Xkr4;gene_id=MGI:3528744;curie=MGI:3528744;Dbxref=ENSEMBL:ENSMUSG00000051951,NCBI_Gene:497097;mgi_type=protein coding gene;so_term_name=protein_coding_gene;cID=MGI:3528744;long_name=X-linked Kx blood group related 4
```

### Transcripts. 
* Transcripts are stored files with names of the form: <genome>/models/transcripts/<chr>/<chunk_num>.gff3
  * there is a chunk size defined for the genome
  * if a transcript crosses a chunk boundary, it appears in both chunk files
* Files are GFF-like:
  * same line format 
  * but the Parent for a transcript refers to a top level feature in the genes file.
  * subfeatures (exons/CDSs) encoded into column 9 of the transcript
* Exon coodinates are encoded as an "exons" attribute in colummn 9 of the transcript feature.
  * formatted as a (comma-separated) list of (underscore-separated) integer pairs specifying delta (from the start of the transcript) and length of each exon. 
  * for example: exons=0_422,7218_150,9607_95
* CDS coodinates are encoded as a "cds" attribute in column 9 of the transcript feature
  * formatted as a pipe separated list of four values: ID of CDS from imported file, protein id, start of first piece of CDS, end of last piece of CDS (these coordinates are absolute, not deltas)
  * for example: cds=MGI_C57BL6J_1915252_cds_2|RefSeq:NP_080779.2|23996251|24005600
* Example:
```
1       NCBI    mRNA    3269956 3741733 .       -       .       ID=MGI_C57BL6J_3528744_transcript_7;Name=XM_006495550.5;Parent=MGI_C57BL6J_3528744_GRCm39;transcript_id=RefSeq:XM_006495550.5;gene_id=MGI:3528744;exons=0_7585,13706_3530,221969_200,470819_959;cds=MGI_C57BL6J_3528744_cds_5|RefSeq:XP_006495613.1|3286245|3741571
```

## CGI script.
* Requests for sequences are served by a CGI script:
  * fetch.cgi is a wrapper shell script; fetch.py does the heavy lifting
  * Both scripts live in the data root directory by default. You can change this by setting CDIR in wrapperConfig.sh
* Request format is documented in fetch.py
* Return format is FASTA

## Homologies
* Homology data live in a "homologies" directory under the data root.
  * Orthology data stored in files named homologies/orthology/<taxonid>.json
  * (eventually) Paralogy data will go under homologies/paralogy
* An orthology data file is a list of ID pairs containing all orthology assertions for one taxon.
  * Each list item has the form: [id1, taxon1, td2, taxon2, flags]. (flags is not currently used by the viewer)
  * Each list item is an assertion that id2/taxon2 is an orthlog of id1/taxon1
  * No inverse relationship is assumed; the inverse must be explicitly provided in taxon2's othology file.
  * Example: the first few lines of homologies/orthology/10090.json
```
[["MGI:1917606", "10090", "FB:FBgn0260484", "7227", "YY"]
,["MGI:1923224", "10090", "FB:FBgn0260486", "7227", "YN"]
,["MGI:106321", "10090", "FB:FBgn0260486", "7227", "YN"]
...
```
  * The file homologies/orthology/7227.json has the inverse relationships, e.g. ["FB:FBgn0260484", "7227", "MGI:1917606", "10090", "YY"]
  * This model derives directly from the Alliance of Genome Resources, strict orthology set.

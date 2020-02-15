#!/usr/bin/env bash
#
# build.sh
#
# Top level build script, specific for building the annotated mouse strains.
# To build everything:
#	$ ./build.sh

# Any command line parameters are passed along to getGenome.
# Examples:
#
# Run the whole thing:
#	$ ./build.sh
#
# Only build for the A/J strain:
#	$ ./build.sh -G mus_musculus_aj
# 
# Only run the GFF3 import phase:
#	$ ./build.sh -t models -p import
#

source config.sh
source utils.sh

logit "mgv_data: starting build..."

logit "mgv_data: getting homolgy data..."
./getHomologies.sh 2>> ${LOGFILE}

# default list of filter/transformation steps to run for mouse genomes
FILTERS="ensembl,stripPrefix,tagEnsemblMgi"
# non-B6 mouse strains
./getGenome.sh -x 10090 -m ${FILTERS} -g mus_musculus_129s1svimj -n 129S1/SvImJ $* >${LOGFILE} 2>&1
./getGenome.sh -x 10090 -m ${FILTERS} -g mus_musculus_aj         -n A/J         $* >${LOGFILE} 2>&1
./getGenome.sh -x 10090 -m ${FILTERS} -g mus_musculus_akrj       -n AKR/J       $* >${LOGFILE} 2>&1
./getGenome.sh -x 10090 -m ${FILTERS} -g mus_musculus_balbcj     -n BALB/cJ     $* >${LOGFILE} 2>&1
./getGenome.sh -x 10090 -m ${FILTERS} -g mus_musculus_c3hhej     -n C3H/HeJ     $* >${LOGFILE} 2>&1
./getGenome.sh -x 10090 -m ${FILTERS} -g mus_musculus_c57bl6nj   -n C57BL/6NJ   $* >${LOGFILE} 2>&1
./getGenome.sh -x 10089 -m ${FILTERS} -g mus_caroli              -n CAROLI/EiJ  $* >${LOGFILE} 2>&1
./getGenome.sh -x 10090 -m ${FILTERS} -g mus_musculus_casteij    -n CAST/EiJ    $* >${LOGFILE} 2>&1
./getGenome.sh -x 10090 -m ${FILTERS} -g mus_musculus_cbaj       -n CBA/J       $* >${LOGFILE} 2>&1
./getGenome.sh -x 10090 -m ${FILTERS} -g mus_musculus_dba2j      -n DBA/2J      $* >${LOGFILE} 2>&1
./getGenome.sh -x 10090 -m ${FILTERS} -g mus_musculus_fvbnj      -n FVB/NJ      $* >${LOGFILE} 2>&1
./getGenome.sh -x 10090 -m ${FILTERS} -g mus_musculus_lpj        -n LP/J        $* >${LOGFILE} 2>&1
./getGenome.sh -x 10090 -m ${FILTERS} -g mus_musculus_nodshiltj  -n NOD/ShiLtJ  $* >${LOGFILE} 2>&1
./getGenome.sh -x 10090 -m ${FILTERS} -g mus_musculus_nzohlltj   -n NZO/HlLtJ   $* >${LOGFILE} 2>&1
./getGenome.sh -x 10093 -m ${FILTERS} -g mus_pahari              -n PAHARI/EiJ  $* >${LOGFILE} 2>&1
./getGenome.sh -x 10090 -m ${FILTERS} -g mus_musculus_pwkphj     -n PWK/PhJ     $* >${LOGFILE} 2>&1
./getGenome.sh -x 10096 -m ${FILTERS} -g mus_spretus             -n SPRET/EiJ   $* >${LOGFILE} 2>&1
./getGenome.sh -x 10090 -m ${FILTERS} -g mus_musculus_wsbeij     -n WSB/EiJ     $* >${LOGFILE} 2>&1

# C57BL/6J is a little different.
# Get the GFF3 from MGI, and use different preprocessing.
# Assembly still comes from Ensembl.
./getGenome.sh -x 10090 -m MGI        -g mus_musculus            -n C57BL/6J --gff-url "${MGI_URL}" $* >${LOGFILE} 2>&1

# Alliance organisms
./getGenome.sh -x 9606   -m ensembl,stripPrefix,tagEnsemblHuman  -g homo_sapiens             -n H.sapiens       $* >${LOGFILE} 2>&1
./getGenome.sh -x 10116  -m ensembl,stripPrefix,tagEnsemblRat    -g rattus_norvegicus        -n R.norvegicus    $* >${LOGFILE} 2>&1
./getGenome.sh -x 7955   -m ensembl,stripPrefix,tagEnsemblFish   -g danio_rerio              -n D.rerio         $* >${LOGFILE} 2>&1
./getGenome.sh -x 6239   -m ensembl,stripPrefix,tagEnsemblWorm   -g caenorhabditis_elegans   -n C.elegans       $* >${LOGFILE} 2>&1
./getGenome.sh -x 7227   -m ensembl,stripPrefix,tagEnsemblFly    -g drosophila_melanogaster  -n D.melanogaster  $* >${LOGFILE} 2>&1
./getGenome.sh -x 559292 -m ensembl,stripPrefix,tagEnsemblYeast  -g saccharomyces_cerevisiae -n "S.cerevisiae (S288C)" $* >${LOGFILE} 2>&1


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

./getHomologies.sh

MODULES="ensembl,stripPrefix,tagEnsemblWithMgi"
./getGenome.sh -x 10090 -m ${MODULES} -g mus_musculus_129s1svimj -n 129S1/SvImJ $*
./getGenome.sh -x 10090 -m ${MODULES} -g mus_musculus_aj         -n A/J $*
./getGenome.sh -x 10090 -m ${MODULES} -g mus_musculus_akrj       -n AKR/J $*
./getGenome.sh -x 10090 -m ${MODULES} -g mus_musculus_balbcj     -n BALB/cJ $*
./getGenome.sh -x 10090 -m ${MODULES} -g mus_musculus_c3hhej     -n C3H/HeJ $*
# C57BL/6J is a little different.
# Get the GFF3 from MGI, and use different preprocessing.
# Assembly still comes from Ensembl.
./getGenome.sh -x 10090 -m MGI        -g mus_musculus            -n C57BL/6J --gff-url "${MGI_URL}" $*
./getGenome.sh -x 10090 -m ${MODULES} -g mus_musculus_c57bl6nj   -n C57BL/6NJ $*
./getGenome.sh -x 10089 -m ${MODULES} -g mus_caroli              -n CAROLI/EiJ $*
./getGenome.sh -x 10090 -m ${MODULES} -g mus_musculus_casteij    -n CAST/EiJ $*
./getGenome.sh -x 10090 -m ${MODULES} -g mus_musculus_cbaj       -n CBA/J $*
./getGenome.sh -x 10090 -m ${MODULES} -g mus_musculus_dba2j      -n DBA/2J $*
./getGenome.sh -x 10090 -m ${MODULES} -g mus_musculus_fvbnj      -n FVB/NJ $*
./getGenome.sh -x 10090 -m ${MODULES} -g mus_musculus_lpj        -n LP/J $*
./getGenome.sh -x 10090 -m ${MODULES} -g mus_musculus_nodshiltj  -n NOD/ShiLtJ $*
./getGenome.sh -x 10090 -m ${MODULES} -g mus_musculus_nzohlltj   -n NZO/HlLtJ $*
./getGenome.sh -x 10093 -m ${MODULES} -g mus_pahari              -n PAHARI/EiJ $*
./getGenome.sh -x 10090 -m ${MODULES} -g mus_musculus_pwkphj     -n PWK/PhJ $*
./getGenome.sh -x 10096 -m ${MODULES} -g mus_spretus             -n SPRET/EiJ $*
./getGenome.sh -x 10090 -m ${MODULES} -g mus_musculus_wsbeij     -n WSB/EiJ $*

# Alliance organisms
./getGenome.sh -x 9606  -m ensembl,stripPrefix,tagEnsemblHuman -g homo_sapiens     -n H.sapiens $*
./getGenome.sh -x 10116 -m ensembl,stripPrefix,tagEnsemblRat  -g rattus_norvegicus -n R.norvegicus $*
./getGenome.sh -x 7955  -m ensembl,stripPrefix,tagEnsemblFish   -g danio_rerio -n D.rerio $*
./getGenome.sh -x 6239  -m ensembl,stripPrefix,tagEnsemblWorm   -g caenorhabditis_elegans -n C.elegans $*
./getGenome.sh -x 7227 -m ensembl,stripPrefix,tagEnsemblFly   -g drosophila_melanogaster -n D.melanogaster $*
./getGenome.sh -x 559292 -m ensembl,stripPrefix,tagEnsemblYeast  -g saccharomyces_cerevisiae -n "S.cerevisiae (S288C)" $*


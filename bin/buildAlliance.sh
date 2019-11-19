#!/usr/bin/env bash
#
# buildAlliance.sh
#

source config.sh
source utils.sh
#
./getHomologies.sh -x 10090 -x 10116 -x 9606 -x 7955 -x 6239 -x 7227 -x 559292 > "${DDIR}/homologies.txt"
#
./getGenome.sh -x 10090 -m MGI,tagEnsemblWithHid -g mus_musculus -n M.musculus --gff-url "${MGI_URL}" $*
./getGenome.sh -x 9096  -m ensembl,stripPrefix,tagEnsemblHuman,tagEnsemblWithHid -g homo_sapiens     -n H.sapiens $*
./getGenome.sh -x 10116 -m ensembl,stripPrefix,tagEnsemblRat,tagEnsemblWithHid  -g rattus_norvegicus -n R.norvegicus $*
./getGenome.sh -x 7955  -m ensembl,stripPrefix,tagEnsemblWithHid  -g danio_rerio -n D.rerio $*
./getGenome.sh -x 6239  -m ensembl,stripPrefix,tagEnsemblWithHid  -g caenorhabditis_elegans -n C.elegans $*
./getGenome.sh -x 7227 -m ensembl,stripPrefix,tagEnsemblWithHid  -g drosophila_melanogaster -n D.melanogaster $*

# fixme: yeast name should include "S288C"
./getGenome.sh -x 559292 -m ensembl,stripPrefix,tagEnsemblWithHid  -g saccharomyces_cerevisiae -n S.cerevisiae $*


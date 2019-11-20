#!/usr/bin/env bash
#
# buildAlliance.sh
#

source config.sh
source utils.sh
#
./getHomologies.sh
#
./getGenome.sh -x 10090 -m MGI -g mus_musculus -n M.musculus --gff-url "${MGI_URL}" $*
./getGenome.sh -x 9606  -m ensembl,stripPrefix,tagEnsemblHuman -g homo_sapiens     -n H.sapiens $*
./getGenome.sh -x 10116 -m ensembl,stripPrefix,tagEnsemblRat  -g rattus_norvegicus -n R.norvegicus $*
./getGenome.sh -x 7955  -m ensembl,stripPrefix,tagEnsemblFish   -g danio_rerio -n D.rerio $*
./getGenome.sh -x 6239  -m ensembl,stripPrefix,tagEnsemblWorm   -g caenorhabditis_elegans -n C.elegans $*
./getGenome.sh -x 7227 -m ensembl,stripPrefix,tagEnsemblFly   -g drosophila_melanogaster -n D.melanogaster $*

# fixme: yeast name should include "S288C"
./getGenome.sh -x 559292 -m ensembl,stripPrefix,tagEnsemblYeast  -g saccharomyces_cerevisiae -n S.cerevisiae $*


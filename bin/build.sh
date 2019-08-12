#!/usr/bin/env bash

source config.sh

RELEASE="97"
MODULES="-m pg_ensembl,pg_tagEnsemblWithMgi"
# 
MGI_URL="http://www.informatics.jax.org/downloads/mgigff3/MGI.gff3.gz"

./getGenome.sh -r ${RELEASE} -x 10090 ${MODULES} -g mus_musculus_129s1svimj -n 129S1/SvImJ $*
./getGenome.sh -r ${RELEASE} -x 10090 ${MODULES} -g mus_musculus_aj -n A/J $*
./getGenome.sh -r ${RELEASE} -x 10090 ${MODULES} -g mus_musculus_akrj -n AKR/J $*
./getGenome.sh -r ${RELEASE} -x 10090 ${MODULES} -g mus_musculus_balbcj -n BALB/cJ $*
./getGenome.sh -r ${RELEASE} -x 10090 ${MODULES} -g mus_musculus_c3hhej -n C3H/HeJ $*
./getGenome.sh -r ${RELEASE} -x 10090 ${MODULES} -g mus_musculus_c57bl6nj -n C57BL/6NJ $*
./getGenome.sh -r ${RELEASE} -x 10089 ${MODULES} -g mus_caroli -n CAROLI/EiJ $*
./getGenome.sh -r ${RELEASE} -x 10090 ${MODULES} -g mus_musculus_casteij -n CAST/EiJ $*
./getGenome.sh -r ${RELEASE} -x 10090 ${MODULES} -g mus_musculus_cbaj -n CBA/J $*
./getGenome.sh -r ${RELEASE} -x 10090 ${MODULES} -g mus_musculus_dba2j -n DBA/2J $*
./getGenome.sh -r ${RELEASE} -x 10090 ${MODULES} -g mus_musculus_fvbnj -n FVB/NJ $*
./getGenome.sh -r ${RELEASE} -x 10090 ${MODULES} -g mus_musculus_lpj -n LP/J $*
./getGenome.sh -r ${RELEASE} -x 10090 ${MODULES} -g mus_musculus_nodshiltj -n NOD/ShiLtJ $*
./getGenome.sh -r ${RELEASE} -x 10090 ${MODULES} -g mus_musculus_nzohlltj -n NZO/HlLtJ $*
./getGenome.sh -r ${RELEASE} -x 10093 ${MODULES} -g mus_pahari -n PAHARI/EiJ $*
./getGenome.sh -r ${RELEASE} -x 10090 ${MODULES} -g mus_musculus_pwkphj -n PWK/PhJ $*
./getGenome.sh -r ${RELEASE} -x 10096 ${MODULES} -g mus_spretus -n SPRET/EiJ $*
./getGenome.sh -r ${RELEASE} -x 10090 ${MODULES} -g mus_musculus_wsbeij -n WSB/EiJ $*
#
# C57BL/6J is a little different. Get the GFF3 from MGI
curl "${MGI_URL}" > ${DDIR}/mus_musculus.gff3.gz
gunzip -c ${DDIR}/mus_musculus.gff3.gz | \
    ${PYTHON} prepGff3.py -x "${EXCLUDE_TYPES}" -c "${CHR_REGEX}" -m pg_MGI | \
    ${PYTHON} importGff3.py -p mus_musculus -g C57BL/6J -x 10090 -k ${K} -d ${ODIR}
# get the assembly from Ensembl
./getGenome.sh -r ${RELEASE} -x 10090 ${MODULES} -g mus_musculus -n C57BL/6J -t assembly


#!/usr/bin/bash

export SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

source ${SCRIPT_DIR}/buildcore.sh

activateVenv
python ${SCRIPT_DIR}/compileConfig.py $1 > $2
chmod +x $2



#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
python2.7 ${DIR}/fetch.py --cgi --dir ${DIR}

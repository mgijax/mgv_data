#
# Set-and-forget parameters for an installation. 
# Mainly directories. Edit as needed. 
#

# BCONFIG: Build config file. Configures the actual build steps, data sources, etc.
export BCONFIG=${SCRIPT_DIR}/buildConfig.yaml

# DDIR: Downloads directory. 
# Where files downloaded from external sources will go.
# Default is under the installation directory.
export DDIR=${SCRIPT_DIR}/mgv_data_files/downloads

# TDIR: Temp directory / work area
export TDIR=/tmp

# ODIR: Output directory. This is where the transformed/internalized results of the import phase go.
# Default is under the installation directory, next to the downloads directory.
export ODIR=${SCRIPT_DIR}/mgv_data_files/output

# WDIR: Web data directory. Web-accessible directory where the files are actually served from.
# Files from the output directory are copied here during the deployment phase.
# The copy is skipped if the two directories are the same, which is the default.
export WDIR=${ODIR}

# CDIR: CGI directory. A web-accessible directory that is allowed to serve CGI scripts.
# Default is the web data directory.
export CDIR=${WDIR}

# Python 3 executable
export PYTHON="python3"

# MouseMine URL - used by some filters
export MOUSEMINE_URL="https://www.mousemine.org"

#
export SAMTOOLSBIN="/home/jer/samtools/build/bin"
export SAMTOOLS="${SAMTOOLSBIN}/samtools"
export TABIX="${SAMTOOLSBIN}/tabix"
export FAIDX="${SAMTOOLS} faidx"
export BGZIP="${SAMTOOLSBIN}/bgzip"
export CURL="curl"
export GUNZIP="gunzip"
export PYTHONPATH="${SCRIPT_DIR}/bin/lib"


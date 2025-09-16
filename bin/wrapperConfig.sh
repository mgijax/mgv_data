#
# Set-and-forget parameters for an installation. 
# Mainly directories. Edit as needed. 
#

# set to true generate development deployments.
# Development deployments disable CORS origin restrictions.
export DEVMODE="false"

# DROOTDIR: Data root directory 
# Default is under the installation directory.
export OROOTDIR=${INSTALL_DIR}/mgv_data_files

# DDIR: Downloads directory. 
# Where files downloaded from external sources will go.
export DDIR=${OROOTDIR}/downloads

# TDIR: Temp directory / work area
export TDIR=${OROOTDIR}/work

# ODIR: Output directory. This is where the transformed/internalized results of the import phase go.
export ODIR=${OROOTDIR}/output

# WDIR: Web data directory. Web-accessible directory where the files are actually served from.
# Files from the output directory are copied here during the deployment phase.
# The copy is skipped if the two directories are the same, which is the default.
export WDIR=${ODIR}

# CDIR: CGI directory. A web-accessible directory that is allowed to serve CGI scripts.
# Default is the web data directory.
export CDIR=${WDIR}

# Python 3 executable
export PYTHON="/opt/python/bin/python3"
export PYTHONPATH="${SCRIPT_DIR}/bin/lib"

# MouseMine URL - used by some filters
export MOUSEMINE_URL="https://www.mousemine.org"

# Tabix/Samtools executables
export TABIX="/usr/local/tabix/bin/tabix"
export SAMTOOLS="/usr/local/tabix/bin/samtools"
export BGZIP="/usr/local/tabix/bin/bgzip"

# Other stuff
export CURL="curl"
export GUNZIP="gunzip"


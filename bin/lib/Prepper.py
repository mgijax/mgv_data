#
# Prepper.py
#
# Filters/munges downloaded data into the form needed by the Importer.
#

from .Filter import filterNameMap

class Prepper:
    def __init__(self, config):
        self.config = config


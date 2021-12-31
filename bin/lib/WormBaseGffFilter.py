from .GffFilter import GffFilter
class WormBaseGffFilter (AllianceGff) :
    def processFeature (self, f) :
        AllianceGff.processFeature(self, f)
        attrs = f[8]
        if 'ID' in attrs and attrs['ID'].startswith('Transcript:'):
            attrs['transcript_id'] = 'WB:' + attrs['ID'][11:]
        return f


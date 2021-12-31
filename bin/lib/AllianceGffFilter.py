from .GffFilter import GffFilter
class AllianceGffFilter (GffFilter) : 
    def processFeature(self, f) :
        attrs = f[8]
        attrs.pop("description", None)
        if "so_term_name" in attrs:
            f[2] = attrs.pop("so_term_name")
        elif "Ontology_term" in attrs:
            soid = None
            soids = list(filter(lambda i: i.startswith("SO:"), attrs["Ontology_term"]))
            if len(soids):
                soid = soids[0]
            if soid == "SO:0001217":
                f[2] = "protein_coding_gene"
            elif soid == "SO:0000336":
                f[2] = "pseudogene"
            elif soid == "SO:0001263":
                f[2] = "ncRNA_gene"
        return f


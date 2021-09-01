import sys
import gzip
from typing import List
from dataclasses import dataclass, field
from xml.etree.ElementTree import iterparse

#path = 'data/demo.xml'
path = sys.argv[1]
if path[-2:] == 'gz':
    handle = gzip.open(path, 'r')
else:
    handle = open(path, 'r')


@dataclass
class Gene:
    """Holds metadata representing a gene"""

    accession: str = ""
    sequence: str = ""
    existence: str = ""
    citations: int = 0
    taxid: int = 0
    names: List[str] = field(default_factory=list)
    ipr: List[str] = field(default_factory=list)
    pfam: List[str] = field(default_factory=list)
    go: List[str] = field(default_factory=list)    
    
    def __repr__(self):
        """Prints the gene as a row"""
        fields = [
            self.accession,
            self.existence,
            self.citations,
            self.taxid,
            "|".join(self.names),
            "|".join(self.ipr),
            "|".join(self.pfam),
            "|".join(self.go),
            self.sequence,
        ]
        fields = [str(fi) for fi in fields]
        row = "\t".join(fields)
        return row

# Write header
print("id\tevidence\tcitations\ttaxid\tname\tIPR\tPFAM\tGO\tseq")
gene = Gene()
for _, elem in iterparse(handle, events=('end',)):
    tag = elem.tag
    # Tags all have a URL prefix between {}, so we trim it
    if '}' in tag:
        tag = tag.split('}', 1)[1]
    if tag == 'entry':
        print(gene)
        gene = Gene()
    elif tag == 'accession':
        gene.accession = elem.text
    elif tag == 'sequence':
        gene.sequence = elem.text
    elif tag == 'proteinExistence':
        gene.existence = elem.attrib['type']
    elif tag == 'citation':
        gene.citations += 1
    elif tag == 'fullName':
        gene.names.append(elem.text)
    elif tag == 'dbReference':
        rtype = elem.attrib['type']
        rid = elem.attrib['id']
        if rtype  == 'NCBI Taxonomy':
            gene.taxid = rid
        elif rtype == 'GO':
            gene.go.append(rid)
        elif rtype == 'InterPro':
            gene.ipr.append(rid)
        elif rtype == 'Pfam':
            gene.pfam.append(rid)
    elem.clear()

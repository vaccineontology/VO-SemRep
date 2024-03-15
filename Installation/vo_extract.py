"""
 /rhome/wjrogers/studio/python/rdf/efo_extract.py, Mon Jan  9 14:05:05 2012,
 edit by Will Rogers

Extract concepts and synomyms from EFO_inferred_v2.18.owl and generate
UMLS format tables for use by MetaMap's Data File Builder.


Original Author: Willie Rogers, 09jan2012
Updated on 03sep2014 for 2014 version of EFO ontology (not inferred.)
"""
import argparse
from rdflib import Namespace, ConjunctiveGraph
from readrdf import readrdf
import re
import sys

from mwi_utilities import normalize_ast_string

from pprint import pprint
import umls_references as umls
from vo_umls_refs import vo_umls_refs


datafile = 'brucellosis_merged.owl' # '/efo.owl'
# gwa = '/net/lhcdevfiler/vol/cgsb5/ind/II_Group_WorkArea'
# collections = gwa + '/wjrogers/Collections'
# efo_datafile = collections + '/EFO/efo.owl'
# EFO = Namespace("http://www.ebi.ac.uk/efo/")
EFO = Namespace("http://purl.obolibrary.org/obo/VO#")
RDFS = Namespace("http://www.w3.org/2000/01/rdf-schema#")
RDF = Namespace('http://www.w3.org/1999/02/22-rdf-syntax-ns#')
OWL = Namespace("http://www.w3.org/2002/07/owl#")
BFO = Namespace("http://www.ifomis.org/bfo/1.1/snap#MaterialEntity")

mmencoding = 'utf-8'


"""
For IDOBRU:-   
    BFO, IDO, OGMS, CHEBI, CLO (Cell line), CL (cell type), 
    GO, IAO, NCBITaxon, OBI, PRO, RO, UBERON, VO, OGG, PR
"""
prefixdict = {
    'http://purl.obolibrary.org/obo/VO#': 'VO',
    'http://purl.org/dc/elements/1.1/': 'dc',
    'http://www.w3.org/2000/01/rdf-schema#': 'rdfs',
    'http://purl.org/obo/owl/PATO#': 'PATO',
    'http://purl.org/obo/owl/OBO_REL#': 'OBO_REL',
    'http://www.orpha.net/ORDO/': 'ORDO',
    'http://xmlns.com/foaf/0.1/': 'foaf',
    'http://purl.obolibrary.org/obo/bto#': 'bto',
    'http://usefulinc.com/ns/doap#': 'doap',
    'http://purl.obolibrary.org/obo#': 'obo2',
    'http://www.co-ode.org/patterns#': 'patterns',
    'http://purl.obolibrary.org/obo/ncbitaxon#': 'ncbitaxon',
    'http://semanticscience.org/resource/': 'resource',
    'http://www.obofoundry.org/ro/ro.owl#': 'ro',
    'http://purl.obolibrary.org/obo/uberon#': 'uberon',
    'http://www.ebi.ac.uk/efo/': 'efo',
    'http://www.geneontology.org/formats/oboInOwl#': 'oboInOwl',
    'http://purl.obolibrary.org/obo/go#': 'go',
    'http://purl.obolibrary.org/obo/uberon/phenoscape-anatomy#':
    'phenoscape-anatomy',
    'http://purl.org/dc/terms/': 'terms',
    'http://www.orphanet.org/rdfns#': 'rdfns',
    'http://www.w3.org/2002/07/owl#': 'owl',
    'http://purl.obolibrary.org/obo/uberon/core#': 'core',
    'http://www.ifomis.org/bfo/1.1/snap#': 'snap',
    'http://www.ifomis.org/bfo/1.1/span#': 'span',
    'http://www.w3.org/1999/02/22-rdf-syntax-ns#': 'rdf',
    'http://www.w3.org/2001/XMLSchema#': 'xsd',
    'http://purl.obolibrary.org/obo/': 'obo',
    'http://purl.obolibrary.org/obo/BFO#': 'bfo',
    'http://www.w3.org/2003/11/swrl#': 'swrl',
    'http://www.w3.org/2003/11/swrlb#': 'swrlb',
    'http://www.w3.org/2006/12/owl2-xml#': 'owl2xml',
    'http://purl.obolibrary.org/obo/ido/': 'ido',
    'http://www.obofoundry.org/ro/ro.owl#': 'ro',
    'http://purl.obolibrary.org/obo/ubprop#': 'ubprop',
    'http://protege.stanford.edu/plugins/owl/protege#': 'protege'
}

prefixlist = [
    'http://purl.obolibrary.org/obo/VO#',
    'http://purl.org/dc/elements/1.1/',
    'http://www.w3.org/2000/01/rdf-schema#',
    'http://purl.org/obo/owl/PATO#',
    'http://purl.org/obo/owl/OBO_REL#',
    'http://www.orpha.net/ORDO/',
    'http://xmlns.com/foaf/0.1/',
    'http://purl.obolibrary.org/obo/bto#',
    'http://usefulinc.com/ns/doap#',
    'http://purl.obolibrary.org/obo#',
    'http://www.co-ode.org/patterns#',
    'http://purl.obolibrary.org/obo/ncbitaxon#',
    'http://semanticscience.org/resource/',
    'http://www.obofoundry.org/ro/ro.owl#',
    'http://purl.obolibrary.org/obo/uberon#',
    'http://www.ebi.ac.uk/efo/',
    'http://www.geneontology.org/formats/oboInOwl#',
    'http://purl.obolibrary.org/obo/go#',
    'http://purl.obolibrary.org/obo/uberon/phenoscape-anatomy#',
    'http://purl.org/dc/terms/',
    'http://www.orphanet.org/rdfns#',
    'http://www.w3.org/2002/07/owl#',
    'http://purl.obolibrary.org/obo/uberon/core#',
    'http://www.ifomis.org/bfo/1.1/snap#',
    'http://www.ifomis.org/bfo/1.1/span#',
    'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
    'http://www.w3.org/2001/XMLSchema#',
    'http://purl.obolibrary.org/obo/',
    'http://purl.obolibrary.org/obo/BFO#',
    'http://www.w3.org/2003/11/swrl#',
    'http://www.w3.org/2003/11/swrlb#',
    'http://www.w3.org/2006/12/owl2-xml#',
    'http://purl.obolibrary.org/obo/ido/',
    'http://www.obofoundry.org/ro/ro.owl#',
    'http://purl.obolibrary.org/obo/ubprop#',
    'http://protege.stanford.edu/plugins/owl/protege#'
    
    ]
srclist = ['VO', 'dc', 'rdfs', 'PATO', 'OBO_REL', 'ORDO', 'foaf', 'bto',
           'doap', 'obo2', 'patterns', 'ncbitaxon', 'resource', 'ro',
           'uberon', 'efo', 'oboInOwl', 'go', 'phenoscape-anatomy', 'terms',
           'rdfns', 'owl', 'core', 'snap', 'span', 'rdf', 'xsd', 'obo', 'bfo',
           'swrl', 'swrlb', 'owl2xml', 'ido', 'ro', 'ubprop', 'protege'
           ]

def list_id_and_labels(graph):
    for s, p, o in graph:
        if p.__str__ == 'http://www.w3.org/2000/01/rdf-schema#label':
            print(s, p, o)


def query_labels(graph):
    ns = dict(efo=EFO, rdfs=RDFS)
    return graph.query('SELECT ?aname ?bname WHERE { ?a rdfs:label ?b }',
                       initNs=ns)


def query_efo_type(graph, typename):
    ns = dict(efo=EFO, rdfs=RDFS)
    return graph.query(
        'SELECT ?aname ?bname WHERE { ?a efo:%s ?b }' % typename,
        initNs=ns)


def query_efo_concepts(graph):
    ns = dict(owl=OWL, rdf=RDF)
    return graph.query('SELECT ?s { ?s rdf:type owl:Class }', initNs=ns)


def query_efo_concept_labels(graph):
    """ return all labels (concepts) from graph (EFO_inferred_v2.18.owl)"""
    ns = dict(owl=OWL, rdf=RDF, rdfs=RDFS)
    return graph.query(
        'SELECT ?s ?o { ?s rdf:type owl:Class . ?s rdfs:label ?o }',
        initNs=ns)


def query_efo_concept_synonyms(graph):
    """return all alternative_terms (synonyms) from graph
    (EFO_inferred_v2.18.owl)"""
    ns = dict(owl=OWL, rdf=RDF, efo=EFO)
    return graph.query(
        'SELECT ?s ?o { ?s rdf:type owl:Class . ?s efo:alternative_term ?o }',
        initNs=ns)


def escape_re_chars(prefix):
    """ escape regular expression special characters in url"""
    return prefix.replace('?', r'\?').replace('#', r'\#')


def abbrev_uri_original(uri):
    " remove prefix from uri leaving unique part of identifier "
    for prefix in prefixlist:
        m = re.match(r"(%s)(.*)" % escape_re_chars(prefix), uri)
        if m is not None:
            return m.group(2).replace(':', '_')
    return uri


def abbrev_uri(urichars):
    " remove prefix from uri leaving unique part of identifier "
    if isinstance(urichars, bytes):
        uri = urichars.decode(mmencoding)
    else:
        uri = urichars
    for prefix in prefixlist:
        # print('prefix = ''%s''' % prefix)
        m = uri.find(prefix)
        if m == 0:
            newuri = uri[len(prefix):].replace(':', '_')
            if newuri.find(':') >= 0:
                print("problem with abbreviated uri: %s" % newuri)
            return newuri
    return uri


def get_source_name_original(uri):
    " derive source name from uri "
    for prefix in prefixdict.keys():
        m = re.match(r"(%s)(.*)" % escape_re_chars(prefix), uri)
        if m is not None:
            return prefixdict[m.group(1)]
    return uri


def get_source_name(uristring):
    " derive source name from uri "
    if isinstance(uristring, bytes):
        uri = uristring.decode(mmencoding)
    else:
        uri = uristring
    for prefix in prefixdict.keys():
        m = uri.find(prefix)
        if m == 0:
            return prefixdict[uri[0:len(prefix)]]
    return uri


def collect_concepts(graph):
    """
    Return dictionaries (maps) of concepts and synonyms
    (alternative_terms) from results of SPARQL queries
    """
    cdict = {}
    conceptresult = query_efo_concept_labels(graph)
    # serialnum = 1
    for row in conceptresult:
        key = tuple(row)[0].__str__()
        if key in cdict:
            cdict[key].append(row)
        else:
            cdict[key] = [row]
    syndict = {}
    synonymresult = query_efo_concept_synonyms(graph)
    for row in synonymresult:
        key = tuple(row)[0].__str__()
        if key in syndict:
            syndict[key].append(row)
        else:
            syndict[key] = [row]
    return cdict, syndict


def is_valid_cui(cui):
    return re.match(r"[A-Za-z]+[\_]*[0-9]+", cui)


def gen_mrcon_original(graph, filename):
    """
    Generate UMLS format MRCON table.

    return rows of the form:
    EFO_0003549|ENG|P|L0000001|PF|S0000001|caudal tuberculum|0|
    EFO_0003549|ENG|S|L0000002|SY|S0000002|posterior tubercle|0|
    EFO_0003549|ENG|S|L0000003|SY|S0000003|posterior tuberculum|0|
    """
    conceptresult = query_efo_concept_labels(graph)
    fp = open(filename, 'w', encoding='utf-8')
    serialnum = 1
    for row in conceptresult:
        if is_valid_cui(abbrev_uri(tuple(row)[0])):
            fp.write("%s|ENG|P|L%07d|PF|S%07d|%s|0|\n" % (
                abbrev_uri(tuple(row)[0]),
                serialnum, serialnum,
                tuple(row)[1]))
            serialnum = serialnum + 1
    fp.close()


def gen_mrcon(filename, cdict={}, syndict={}, strdict={}, luidict={}):
    """
    Generate UMLS format MRCON (concepts) table.

    return rows of the form:
    EFO_0003549|ENG|P|L0000001|PF|S0000001|caudal tuberculum|0|
    EFO_0003549|ENG|S|L0000002|SY|S0000002|posterior tubercle|0|
    EFO_0003549|ENG|S|L0000003|SY|S0000003|posterior tuberculum|0|
    """
    fp = open(filename, 'w', encoding='utf-8')
    serialnum = 1
    for key in cdict.keys():
        for row in cdict[key]:
            if is_valid_cui(abbrev_uri(tuple(row)[0])):
                term = tuple(row)[1].strip()
                lui = get_lui(luidict, term)
                if lui is None:
                    sys.stderr.write(
                        'gen_mrso: LUI missing for prefname %s\n' % (term))
                sui = strdict.get(term)
                if sui is None:
                    sys.stderr.write(
                        'gen_mrcon:SUI missing for preferred name %s\n' %
                        (term))
                fp.write(
                    "%s|ENG|P|%s|PF|%s|%s|0|\n" % (abbrev_uri(tuple(row)[0]),
                                                   lui, sui, term))
                serialnum = serialnum + 1
        if key in syndict:
            if is_valid_cui(abbrev_uri(tuple(row)[0])):
                for row in syndict[key]:
                    term = tuple(row)[1].strip()
                    lui = get_lui(luidict, term)
                    if lui is None:
                        sys.stderr.write(
                            'gen_mrso: LUI missing for synonym %s\n' % (term))
                    sui = strdict.get(term)
                    if sui is None:
                        sys.stderr.write(
                            'gen_mrcon:SUI missing for synonym %s\n' % (term))
                    fp.write("%s|ENG|S|%s|SY|%s|%s|0|\n" %
                             (abbrev_uri(tuple(row)[0]),
                              lui, sui, term))
                    serialnum = serialnum + 1
    fp.close()


def gen_mrso(filename, cdict={}, syndict={}, strdict={}, luidict={}):
    """
    Generate UMLS format MRSO (sources) table.

   EFO_0003549|L0000001|S0000001|EFO|PT|http://www.ebi.ac.uk/efo/EFO_0003549|0|
   EFO_0003549|L0000002|S0000002|EFO|SY|http://www.ebi.ac.uk/efo/EFO_0003549|0|
   EFO_0003549|L0000003|S0000003|EFO|SY|http://www.ebi.ac.uk/efo/EFO_0003549|0|
   """
    fp = open(filename, 'w', encoding='utf-8')
    serialnum = 1
    for key in cdict.keys():
        for row in cdict[key]:
            if is_valid_cui(abbrev_uri(tuple(row)[0])):
                term = tuple(row)[1].strip()
                lui = get_lui(luidict, term)
                if lui is None:
                    sys.stderr.write(
                        'gen_mrso: LUI missing for prefname %s\n' % (term))
                sui = strdict.get(term)
                if sui is None:
                    sys.stderr.write(
                        'gen_mrso: SUI missing for prefname %s\n' % (term))
                fp.write("%s|%s|%s|%s|PT|%s|0|\n" % (
                    abbrev_uri(tuple(row)[0]),
                    lui, sui,
                    get_source_name(tuple(row)[0]),
                    tuple(row)[0]))
                serialnum = serialnum + 1
        if key in syndict:
            for row in syndict[key]:
                if is_valid_cui(abbrev_uri(tuple(row)[0])):
                    term = tuple(row)[1].strip()
                    lui = get_lui(luidict, term)
                    if lui is None:
                        sys.stderr.write(
                            'gen_mrso: LUI missing for synonym %s\n' % (term))
                    sui = strdict.get(term)
                    if sui is None:
                        sys.stderr.write(
                            'SUI missing for synonym %s\n' % (term))
                    fp.write("%s|%s|%s|%s|SY|%s|0|\n" % (
                        abbrev_uri(tuple(row)[0]),
                        lui, sui,
                        get_source_name(tuple(row)[0]),
                        tuple(row)[0]))
                    serialnum = serialnum + 1
    fp.close()


def gen_mrconso(filename, cdict={}, syndict={}, auidict={},
                strdict={}, luidict={}):
    """
    Generate UMLS RRF format MRCONSO (concept+sources) table

    cui|lat|ts|lui|stt|sui|ispref|aui|saui|scui|sdui|sab|tty|code|str|srl|suppress|cvf

    EFO0003549|ENG|P|L0000001|PF|S0000001|Y|A0003549||||EFO|PT|http://www.ebi.ac.uk/efo/EFO_0003549|caudal tuberculum|0|?|N||
    """
    fp = open(filename, 'w', encoding='utf-8')
    serialnum = 1
    for key in cdict.keys():
        for row in cdict[key]:
            uri = tuple(row)[0].strip()
            if is_valid_cui(abbrev_uri(uri)):
                sab = get_source_name(uri)
                term = tuple(row)[1].strip()
                aui = auidict.get((term, sab))
                if aui is None:
                    sys.stderr.write(
                        'gen_mrconso:AUI missing for preferred name %s, %s\n' %
                        (term, sab))
                lui = get_lui(luidict, term)
                if lui is None:
                    sys.stderr.write(
                        'gen_mrso: LUI missing for prefname %s\n' % term)
                sui = strdict.get(term)
                if sui is None:
                    sys.stderr.write(
                        'gen_mrconso:SUI missing for preferred name %s\n' %
                        term)
                fp.write("%s|ENG|P|%s|PF|%s|Y|%s||||%s|PT|%s|%s|0|N||\n" %
                         (abbrev_uri(uri), lui, sui, aui, sab, uri, term))
                serialnum = serialnum + 1
        if key in syndict:
            uri = tuple(row)[0].strip()
            if is_valid_cui(abbrev_uri(uri)):
                for row in syndict[key]:
                    sab = get_source_name(uri)
                    term = tuple(row)[1].strip()
                    aui = auidict.get((term, sab))
                    if aui is None:
                        sys.stderr.write(
                            'gen_mrconso:AUI missing for synonym %s,%s\n' %
                            (term, sab))
                    lui = get_lui(luidict, term)
                    if lui is None:
                        sys.stderr.write(
                            'gen_mrso: LUI missing for synonym %s\n' % (term))
                    sui = strdict.get(term)
                    if sui is None:
                        sys.stderr.write(
                            'gen_mrconso:SUI missing for synonym %s\n' %
                            (term))
                    fp.write("%s|ENG|S|%s|SY|%s|Y|%s||||%s|PT|%s|%s|0|N||\n" %
                             (abbrev_uri(uri), lui, sui,
                              aui, sab, uri, term))
                    serialnum = serialnum + 1

    fp.close()

# semtypes = [('T045', 'Genetic Function', 'genf'),
#             ('T028', 'Gene or Genome', 'gegm'),
#             ('T116', 'Amino Acid, Peptide, or Protein', 'aapp'),
#             ('T121', 'Pharmacologic Substance', 'phsu')]
    
"""
For IDOBRU:-   
    BFO /, IDO /, OGMS /, CHEBI /, CLO (Cell line) //, CL (cell type) /, 
    GO /, IAO /, NCBITaxon /, OBI /, PRO (PR) /, RO //, UBERON /, VO /, OGG /, PR /
"""

semtypes = {
    'VO': {'typeid': 'T121', 'abbrev': 'phsu', "typename": "Pharmacologic Substance", 'treenum': "A1.4.1.1.1"},
    'CHEBI': {'typeid': 'T103', 'abbrev': 'chem', "typename": "Chemical", 'treenum': "A1.4.1"},
    'UO': {'typeid': 'T081', 'abbrev': 'qnco', "typename": "Quantitative Concept", 'treenum': "A2.1.3"},
    'PR': {'typeid': 'T116', 'abbrev': 'aapp', "typename": "Amino Acid, Peptide, or Protein", 'treenum': "A1.4.1.2.1.7"},
    'PATO': {'typeid': 'T081', 'abbrev': 'qnco', "typename": "Quantitative Concept", 'treenum': "A2.1.3"},
    'OPL': {'typeid': 'T038', 'abbrev': 'biof', "typename": "Biologic Function", 'treenum': "B2.2.1"},
    'IAO': {'typeid': 'T170', 'abbrev': 'inpr', "typename": "Intellectual Product", 'treenum': "A2.4"},
    'OPMS': {'typeid': 'T078', 'abbrev': 'idcn', "typename": "Idea or Concept", 'treenum': "A2.1"},
    'OGMS': {'typeid': 'T078', 'abbrev': 'idcn', "typename": "Idea or Concept", 'treenum': "A2.1"},
    'RO': {'typeid': 'T078', 'abbrev': 'idcn', "typename": "Idea or Concept", 'treenum': "A2.1"},
    'OBI':  {'typeid': 'T062', 'abbrev': 'resa', "typename": "Research Activity", 'treenum': "B1.3.2"},
    'BFO': {'typeid': 'T078', 'abbrev': 'idcn', "typename": "Idea or Concept", 'treenum': "A2.1"},
    'OGG': {'typeid': 'T028', 'abbrev': 'gngm', "typename": "Gene or Genome", 'treenum': "A1.2.3.5"},
    'GO': {'typeid': 'T028', 'abbrev': 'gngm', "typename": "Gene or Genome", 'treenum': "A1.2.3.5"},
    'TRANS': {'typeid': 'T047', 'abbrev': 'dsyn', "typename": "Disease or Syndrome", 'treenum': "B2.2.1.2.1"},
    'IDO': {'typeid': 'T047', 'abbrev': 'dsyn', "typename": "Disease or Syndrome", 'treenum': "B2.2.1.2.1"},
    'MONDO': {'typeid': 'T047', 'abbrev': 'dsyn', "typename": "Disease or Syndrome", 'treenum': "B2.2.1.2.1"},
    'DOID': {'typeid': 'T047', 'abbrev': 'dsyn', "typename": "Disease or Syndrome", 'treenum': "B2.2.1.2.1"},
    'FMA': {'typeid': 'T023', 'abbrev': 'bpoc', "typename": "Body Part, Organ, or Organ Component", 'treenum': "A1.2.3.1"},
    'UBERON': {'typeid': 'T023', 'abbrev': 'bpoc', "typename": "Body Part, Organ, or Organ Component", 'treenum': "A1.2.3.1"},
    'CARO': {'typeid': 'T023', 'abbrev': 'bpoc', "typename": "Body Part, Organ, or Organ Component", 'treenum': "A1.2.3.1"},
    'CL': {'typeid': 'T025', 'abbrev': 'cell', "typename": "Cell", 'treenum': "A1.2.3.3"},
    'CLO': {'typeid': 'T025', 'abbrev': 'cell', "typename": "Cell", 'treenum': "A1.2.3.3"},
    'NCBITaxon': {'typeid': 'T001', 'abbrev': 'orgm', "typename": "Organism", 'treenum': "A1.1"},
    'OAE': {'typeid': 'T046', 'abbrev': 'patf', "typename": "Pathologic Function", 'treenum': "B2.2.1.2"},
    'ENVO': {'typeid': 'T083', 'abbrev': 'geoa', "typename": "Geographic Area", 'treenum': "A2.1.5.4"},
}

# semtypes = {
#     'VO': {'typeid': 'T121', 'abbrev': 'phsu', "typename": "Pharmacologic Substance", 'treenum': "A1.4.1.1.1"},
#     'CHEBI': {'typeid': 'T103', 'abbrev': 'chem', "typename": "Chemical", 'treenum': "A1.4.1"},
#     'UO': {'typeid': 'T081', 'abbrev': 'qnco', "typename": "Quantitative Concept", 'treenum': "A2.1.3"},
#     'PR': {'typeid': 'T116', 'abbrev': 'aapp', "typename": "Amino Acid, Peptide, or Protein", 'treenum': "A1.4.1.2.1.7"},
#     'PATO': {'typeid': 'T081', 'abbrev': 'qnco', "typename": "Idea?", 'treenum': "A2.1.3"},
#     'OPL': {'typeid': 'T071', 'abbrev': 'enty', "typename": "Entity", 'treenum': "A"},
#     'IAO': {'typeid': 'T077', 'abbrev': 'cnce', "typename": "Conceptual Entity", 'treenum': "A2"},
#     'OPMS': {'typeid': 'T078', 'abbrev': 'idcn', "typename": "Idea or Concept", 'treenum': "A2.1"},
#     'OGMS': {'typeid': 'T078', 'abbrev': 'idcn', "typename": "Idea or Concept", 'treenum': "A2.1"},
#     'RO': {'typeid': 'T078', 'abbrev': 'idcn', "typename": "Idea or Concept", 'treenum': "A2.1"},
#     'OBI':  {'typeid': 'T078', 'abbrev': 'idcn', "typename": "Idea or Concept", 'treenum': "A2.1"},
#     'BFO': {'typeid': 'T078', 'abbrev': 'idcn', "typename": "Idea or Concept", 'treenum': "A2.1"},
#     'OGG': {'typeid': 'T028', 'abbrev': 'gngm', "typename": "Gene or Genome", 'treenum': "A1.2.3.5"},
#     'GO': {'typeid': 'T028', 'abbrev': 'gngm', "typename": "Gene or Genome", 'treenum': "A1.2.3.5"},
#     'IDO': {'typeid': 'T047', 'abbrev': 'dsyn', "typename": "Disease or Syndrome", 'treenum': "B2.2.1.2.1"},
#     'MONDO': {'typeid': 'T047', 'abbrev': 'dsyn', "typename": "Disease or Syndrome", 'treenum': "B2.2.1.2.1"},
#     'DOID': {'typeid': 'T047', 'abbrev': 'dsyn', "typename": "Disease or Syndrome", 'treenum': "B2.2.1.2.1"},
#     'FMA': {'typeid': 'T017', 'abbrev': 'anst', "typename": "Anatomical Structure", 'treenum': "A1.2"},
#     'UBERON': {'typeid': 'T017', 'abbrev': 'anst', "typename": "Anatomical Structure", 'treenum': "A1.2"},
#     'CL': {'typeid': 'T025', 'abbrev': 'cell', "typename": "Cell", 'treenum': "A1.2.3.3"},
#     'CLO': {'typeid': 'T025', 'abbrev': 'cell', "typename": "Cell", 'treenum': "A1.2.3.3"},
#     'NCBITaxon': {'typeid': 'T008', 'abbrev': 'anim', "typename": "Animal", 'treenum': "A1.1.3.1"},
#     'OAE': {'typeid': 'T046', 'abbrev': 'patf', "typename": "Lab or Test Results", 'treenum': "B2.2.1.2"},
# }

graph = ConjunctiveGraph()
graph.parse(datafile)

vo_umls_refs = umls.test_vo(graph)

def get_onto_type(uri):
    term = uri.split('/')[-1].split('_')[0]
    # print(semtypes[term])
    return term


def get_semantic_typeid(uri):
    """ return semantic type id for uri, currently all uris belong to
    the unknown semantic type. """
    # print(uri)
    return semtypes[uri]["typeid"] # t205


def get_semantic_typeabbrev(uri):
    """ return semantic type abbreviation for uri, currently all uris
    belong to the unknown semantic type."""
    return semtypes[uri]["abbrev"] # unkn


def get_semantic_typename(uri):
    """ return semantic type name for uri, currently all uris
    belong to the unknown semantic type."""
    # print(uri)
    return semtypes[uri]["typename"] # Unknown


def get_semantic_typetree_number(uri):
    """ return semantic tree number for uri, currently all uris
    belong to the unknown semantic type."""
    return semtypes[uri]["treenum"] # "A0.0.0.0.0.0"


def get_semantic_typeui(uri):
    """ return semantic tree number for uri, currently all uris
    belong to the unknown semantic type."""
    return "AT0000000"

def get_cvf(uri):
    """ return semantic tree number for uri, currently all uris
    belong to the unknown semantic type."""
    return "256"

def check_vo_in_umls(uri):
    # print(uri)
    for refs in vo_umls_refs:
        if uri in refs:
            return True
    return False

def gen_mrsty(filename, cdict={}, syndict={}):
    """ Generate UMLS ORF format MRSTY (semantic type) table.  Currently,
    all of the concepts are assigned the semantic type "unkn". """
    fp = open(filename, 'w', encoding='utf-8')
    for key in cdict.keys():
        for row in cdict[key]:
            if is_valid_cui(abbrev_uri(tuple(row)[0])):
                # if concept has UMLS refs, skip it since it's already in 2023AB MRSTY
                if check_vo_in_umls(abbrev_uri(tuple(row)[0])):
                    # print(tuple(row)[0])
                    continue
                onto = get_onto_type(tuple(row)[0])
                fp.write("%s|%s|%s|\n" %
                         (abbrev_uri(tuple(row)[0] ),
                          get_semantic_typeid(onto),
                          get_semantic_typeabbrev(onto)))
    fp.close()



def gen_mrsty_rrf(filename, cdict={}, syndict={}):
    """ Generate UMLS ORF format MRSTY (semantic type) table.  Currently,
    all of the concepts are assigned the semantic type "unkn".

    MRSTY.RFF contaons lines like
    C0000005|T116|A1.4.1.2.1.7|Amino Acid, Peptide, or Protein|AT17648347||
    C0000005|T121|A1.4.1.1.1|Pharmacologic Substance|AT17575038||
    C0000005|T130|A1.4.1.1.4|Indicator, Reagent, or Diagnostic Aid|AT17634323||
    C0000039|T119|A1.4.1.2.1.9|Lipid|AT17617573|256|
    C0000039|T121|A1.4.1.1.1|Pharmacologic Substance|AT17567371|256|
    """
    fp = open(filename, 'w', encoding='utf-8')
    for key in cdict.keys():
        for row in cdict[key]:
            if is_valid_cui(abbrev_uri(tuple(row)[0])):
                # if concept has UMLS refs, skip it since it's already in 2023AB MRSTY
                if check_vo_in_umls(abbrev_uri(tuple(row)[0])):
                    continue
                onto = get_onto_type(tuple(row)[0])
                fp.write("%s|%s|%s|%s|%s|%s|\n" %
                         (abbrev_uri(tuple(row)[0]),
                          get_semantic_typeid(onto),
                          get_semantic_typetree_number(onto),
                          get_semantic_typename(onto),
                          get_semantic_typeui(onto),
                          get_cvf(onto))) 
    fp.close()


def gen_mrsat(filename, cdict={}, syndict={}, auidict={},
              strdict={}, luidict={}):
    """ Generate UMLS format MRSAT (Simple Concept and String
    Attributes) table.  Currently, empty. """
    fp = open(filename, 'w', encoding='utf-8')
    fp.close()


def gen_mrsab(filename, cdict={}, syndict={}):
    """ Generate UMLS format MRSAB (Source Informatino) table.
    Currently, empty. """
    cui_index = 4000000
    fp = open(filename, 'w', encoding='utf-8')
    for k, v in prefixdict.items():
        rcui = vcui = 'C%7s' % cui_index
        vsab = rsab = v
        son = k
        sf = vstart = vend = imeta = rmeta = slc = scc = ssn = scit = ''
        srl = '0'
        if len(cdict) > 0:
            # count concepts that belong to source
            cres = list(filter(lambda x: re.match(r"(%s)(.*)" % k,
                                                  x[0].__str__()),
                               cdict.items()))
            cfr = '%d' % len(cres)
            if len(syndict) > 0:
                # count synonyms that belong to source
                sres = list(filter(lambda x: re.match(r"(%s)(.*)" % k,
                                                      x[0].__str__()),
                                   syndict.items()))
                tfr = '%d' % (len(cres)+len(sres))
            else:
                tfr = '%d' % len(cres)
        else:
            cfr = tfr = ''
        cxty = ttyl = atnl = ''
        lat = 'ENG'
        cenc = 'ascii'
        curver = sabin = 'Y'
        fp.write('%s\n' %
                 '|'.join((vcui, rcui, vsab, rsab, son, sf, vstart, vend,
                           imeta, rmeta, slc, scc, srl, tfr, cfr, cxty,
                           ttyl, atnl, lat, cenc, curver, sabin, ssn, scit)))
        cui_index = cui_index + 1
    fp.close()


def gen_mrrank(filename):
    """ Generate UMLS format MRRANK (Concept Name Ranking) table. """
    ttylist = ['PT', 'SY']
    pred = 400
    fp = open(filename, 'w', encoding='utf-8')
    for sab in srclist:
        for tty in ttylist:
            fp.write('%04d|%s|%s|N|\n' % (pred, sab, tty))
            pred = pred - 1
    fp.close()


def print_result(result):
    for row in result:
        print("%s|%s" % (tuple(row)[0], tuple(row)[1]))


def write_result(result, filename):
    f = open(filename, 'w', encoding='utf-8')
    for row in result:
        f.write(('%s\n' % '|'.join((tuple(row)[0],
                                    tuple(row)[1]))))
    f.close()


def gen_mrcon_list(cdict={}, syndict={}):
    """
    return rows of the form:
    EFO_0003549|ENG|P|L0000001|PF|S0000001|caudal tuberculum|0|
    EFO_0003549|ENG|S|L0000002|SY|S0000002|posterior tubercle|0|
    EFO_0003549|ENG|S|L0000003|SY|S0000003|posterior tuberculum|0|
    """
    mrconlist = []
    serialnum = 1
    for key in cdict.keys():
        for row in cdict[key]:
            if is_valid_cui(abbrev_uri(tuple(row)[0])):
                # "%s|ENG|P|L%07d|PF|S%07d|%s|0|\n"
                mrconlist.append((abbrev_uri(tuple(row)[0]), 'ENG', 'P',
                                  'L%07d' % serialnum, 'PF',
                                  'S%07d' % serialnum,
                                  tuple(row)[1], '0', ''))
                serialnum = serialnum + 1
        if key in syndict:
            for row in syndict[key]:
                if is_valid_cui(abbrev_uri(tuple(row)[0])):
                    # "%s|ENG|S|L%07d|SY|S%07d|%s|0|\n"
                    mrconlist.append((abbrev_uri(tuple(row)[0]), 'ENG', 'P',
                                      'L%07d' % serialnum, 'SY',
                                      'S%07d' % serialnum,
                                      tuple(row)[1], '0', ''))
                    serialnum = serialnum + 1
    return mrconlist


def gen_mrconso_list(cdict={}, syndict={}, auidict={}):
    """
    return rows of the form:
    EFO_0003549|ENG|P|L0000001|PF|S0000001||||caudal tuberculum|0|
    EFO_0003549|ENG|S|L0000002|SY|S0000002|posterior tubercle|0|
    EFO_0003549|ENG|S|L0000003|SY|S0000003|posterior tuberculum|0|
    """

    # mrconlist = []
    # serialnum = 1
    # for key in cdict.keys():
    #     for row in cdict[key]:
    #         if is_valid_cui(abbrev_uri(tuple(row)[0])):
    pass


def gen_strdict(cdict, syndict):
    """ Generate dict of concept and synonym triples mapped by string. """
    strdict = {}
    for triplelist in cdict.values():
        for triple in triplelist:
            strkey = triple[1].__str__()
            if strkey in strdict:
                strdict[strkey].append(triple)
            else:
                strdict[strkey] = [triple]
    for triplelist in syndict.values():
        for triple in triplelist:
            strkey = triple[1].__str__()
            if strkey in strdict:
                strdict[strkey].append(triple)
            else:
                strdict[strkey] = [triple]
    return strdict


def gen_nmstrdict(cdict, syndict):
    """Generate dict of concept and synonym triples mapped by nomalized
       string."""
    strdict = {}
    for triplelist in cdict.values():
        for triple in triplelist:
            strkey = normalize_ast_string(triple[1].__str__())
            if strkey in strdict:
                strdict[strkey].append(triple)
            else:
                strdict[strkey] = [triple]
    for triplelist in syndict.values():
        for triple in triplelist:
            strkey = normalize_ast_string(triple[1].__str__())
            if strkey in strdict:
                strdict[strkey].append(triple)
            else:
                strdict[strkey] = [triple]
    return strdict


def gen_strdict_histogram(strdict):
    """ Generate histrogram of lengths of string dictionary values. """
    histogram = {}
    for v in strdict.values():
        key = '%d' % len(v)
        if key in histogram:
            histogram[key] += 1
        else:
            histogram[key] = 1
    return histogram


def gen_strdict_listsizedict(strdict):
    sizedict = {}
    for k, v in strdict.items():
        key = '%d' % len(v)
        if key in sizedict:
            sizedict[key].append(k)
        else:
            sizedict[key] = [k]
    return sizedict


def gen_aui_dict(cdict=[], syndict=[], auiprefix='A', offset=0):
    """ A simple way to generate atom unique identifiers (AUIS):

    1. Generate list of strings + vocabulary source from ontology.
    2. Sort list
    3. assign auis in descending order of sorted list.


    cdict: concept dictionary
    syndict: synonym dictionary
    auiprefix: prefix for Atom identifiers,
               usually "A" for standalone DataSets,
               should be "B" for dataset to be used with UMLS.
               A can be used if range new identifier space is outside
               of existing UMLS atom identifier space.
    offset=start of range for identifiers, default is zero
    """
    aset = set([])
    auidict = {}
    for cstr in cdict.keys():
        if cstr == 'http://www.ebi.ac.uk/efo/EFO_0000694':
            print("%s -> %s" % (cstr, 'is SARS '))
        prefterm = cdict[cstr][0][1].strip()
        sab = get_source_name(cstr)
        if prefterm == 'SARS':
            print('%s --> pref: %s,%s' % (cstr, prefterm, sab))
        aset.add((prefterm, sab))
        if cstr in syndict:
            for row in syndict[cstr]:
                synonym = row[1].strip()
                if synonym == 'SARS':
                    print('%s --> syn: %s,%s' % (cstr, synonym, sab))
                sab = get_source_name(cstr)
                aset.add((synonym, sab))
    alist = [x for x in aset]
    alist.sort()
    i = offset
    for atom in alist:
        auidict[atom] = '%s%08d' % (auiprefix, i)
        i = i + 1
    return auidict


def gen_sui_dict(cdict=[], syndict=[], suiprefix='S', offset=0):
    """ A simple way to generate String Unique Identifiers(SUIS):

    1. Generate list of strings + vocabulary source from ontology.
    2. Sort list
    3. assign auis in descending order of sorted list.

    cdict: concept dictionary
    syndict: synonym dictionary
    suiprefix: prefix for string identifiers,
              usually "S" for standalone DataSets,
              should be "T" for dataset to be used with UMLS,
              "S" can be used if range new identifier space is outside
              of existing UMLS string identifier space.
    offset=start of range for identifiers, default is zero
    """
    sset = set([])
    suidict = {}
    for cstr in cdict.keys():
        if cstr == 'http://www.ebi.ac.uk/efo/EFO_0000694':
            print("%s -> %s" % (cstr, 'is SARS '))
        prefterm = cdict[cstr][0][1].strip()
        sset.add(prefterm)
        if prefterm == 'SARS':
            print('%s --> pref: %s' % (cstr, prefterm))
    for cstr in syndict.keys():
        if cstr == 'http://www.ebi.ac.uk/efo/EFO_0000694':
            print("%s -> %s" % (cstr, 'is SARS '))
        for row in syndict[cstr]:
            synonym = row[1].strip()
            sset.add(synonym)
            if synonym == 'SARS':
                print('%s --> syn: %s' % (cstr, synonym))
    slist = [x for x in sset]
    slist.sort()
    i = offset
    for mstring in slist:
        suidict[mstring] = '%s%08d' % (suiprefix, i)
        i = i + 1
    return suidict


def gen_lui_dict(cdict=[], syndict=[], luiprefix='L', offset=0):
    """ A simple way to generate Lexical Unique Identifiers(SUIS):

    1. Generate list of strings + vocabulary source from ontology.
    2. Sort list
    3. assign auis in descending order of sorted list.

    cdict: concept dictionary
    syndict: synonym dictionary
    suiprefix: prefix for string identifiers,
              usually "L" for standalone DataSets,
              should be "M" for dataset to be used with UMLS,
              "L" can be used if range new identifier space is outside
              of existing UMLS string identifier space.
    offset=start of range for identifiers, default is zero
    """
    nasset = set([])
    luidict = {}
    for cstr in cdict.keys():
        if cstr == 'http://www.ebi.ac.uk/efo/EFO_0000694':
            print("%s -> %s" % (cstr, 'is SARS '))
        prefterm = cdict[cstr][0][1].strip()
        nasset.add(normalize_ast_string(prefterm))
        if prefterm == 'SARS':
            print('%s --> pref: %s' % (cstr, prefterm))
    for cstr in syndict.keys():
        if cstr == 'http://www.ebi.ac.uk/efo/EFO_0000694':
            print("%s -> %s" % (cstr, 'is SARS '))
        for row in syndict[cstr]:
            synonym = row[1].strip()
            nasset.add(normalize_ast_string(synonym))
            if synonym == 'SARS':
                print('%s --> syn: %s' % (cstr, synonym))
    naslist = [x for x in nasset]
    naslist.sort()
    i = offset
    for nasstring in naslist:
        luidict[nasstring] = '%s%08d' % (luiprefix, i)
        i = i + 1
    return luidict


def get_lui(luidict, mstring):
    """ get lui for un-normalized string from lui dictionary """
    return luidict.get(normalize_ast_string(mstring), 'LUI unknown')


def print_couples(alist):
    for el in alist:
        print("%s: %s" % (el[0].__str__(), el[1].__str__()))


def process(rdffilename):
    graph = readrdf(rdffilename)
    print('finding concepts and synonyms')
    cdict, syndict = collect_concepts(graph)
    print('Generating Atom Unique Identifier Dictionary')
    auidict = gen_aui_dict(cdict, syndict)
    print('Generating String Unique Identifier Dictionary')
    suidict = gen_sui_dict(cdict, syndict)
    print('Generating Lexical Unique Identifier Dictionary')
    luidict = gen_lui_dict(cdict, syndict)

    # rrf
    print('generating MRCONSO.RRF')
    gen_mrconso('MRCONSO.RRF', cdict, syndict, auidict, suidict, luidict)

    # orf
    print('generating MRCON')
    gen_mrcon('MRCON', cdict, syndict, suidict, luidict)
    print('generating MRSO')
    gen_mrso('MRSO', cdict, syndict, suidict, luidict)

    # both rrf and orf
    print('generating MRSAB')
    gen_mrsab('MRSAB.RRF', cdict, syndict)
    print('generating MRRANK')
    gen_mrrank('MRRANK.RRF')
    print('generating MRSAT')
    gen_mrsat('MRSAT.RRF', cdict, syndict)
    print('generating MRSTY')
    gen_mrsty('MRSTY', cdict, syndict)
    print('generating MRSTY.RRF')
    gen_mrsty_rrf('MRSTY.RRF', cdict, syndict)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='make embeddings from sentence list')
    parser.add_argument('efo_file', help='EFO rdf xml file.')
    args = parser.parse_args()
    print('reading %s' % args.efo_file)
    process(args.efo_file)

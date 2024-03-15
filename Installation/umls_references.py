from urllib.parse import urlparse
from rdflib import Namespace


# EFO = Namespace("http://www.ebi.ac.uk/efo/")
EFO = Namespace("http://purl.obolibrary.org/obo/merged/VO")
RDFS = Namespace("http://www.w3.org/2000/01/rdf-schema#")
RDF = Namespace('http://www.w3.org/1999/02/22-rdf-syntax-ns#')
OWL = Namespace("http://www.w3.org/2002/07/owl#")
BFO = Namespace("http://www.ifomis.org/bfo/1.1/snap#MaterialEntity")
SKOS = Namespace("http://www.w3.org/2004/02/skos/core#")
MONDO = Namespace("http://purl.obolibrary.org/obo/mondo#")
NCBITAXON = Namespace("http://purl.obolibrary.org/obo/ncbitaxon#")
OBOINOWL = Namespace("http://www.geneontology.org/formats/oboInOwl#")

def uris_to_ids(reflist):
    """ convert uris in tuples to ids """
    idlist = []
    for tuple in reflist:
        efo_uri_result = urlparse(tuple[0])
        umls_uri_result = urlparse(tuple[1])
        efo_id = efo_uri_result.path.split('/')[-1]
        cui = umls_uri_result.path.split('/')[-1]
        if cui.find("UMLS_CUI") >= 0:
            cui = cui.split(':')[-1]
        # idlist[efo_id] = {"CUI": cui}
        # Call API to fetch semtype info and store in list, if not available in list, call API
        idlist.append([efo_id, cui])
    return idlist


# def collect_skos_efo_to_umls_refs(graph):
#     ns = dict(owl=OWL, rdf=RDF, efo=EFO, skos=SKOS)
#     skos_umls_refs = graph.query(
#         'SELECT ?s ?o { ?s rdf:type owl:Class . ?s skos:exactMatch ?o }',
#         initNs=ns)
#     skos_umls_reflist = [x for x in skos_umls_refs if x[1].find('umls') >= 0]
#     return uris_to_ids(skos_umls_reflist)


# def collect_mondo_efo_to_umls_refs(graph):
#     ns = dict(owl=OWL, rdf=RDF, efo=EFO, mondo=MONDO)
#     mondo_umls_refs = graph.query(
#         'SELECT ?s ?o { ?s rdf:type owl:Class . ?s mondo:exactMatch ?o }',
#         initNs=ns)
#     mondo_umls_reflist = [x for x in mondo_umls_refs if x[1].find('umls') >= 0]
#     return uris_to_ids(mondo_umls_reflist)

# def collect_ncbitaxon_efo_to_umls_refs(graph):
#     ns = dict(owl=OWL, rdf=RDF, efo=EFO, ncbitaxon=NCBITAXON)
#     mondo_umls_refs = graph.query(
#         'SELECT ?s ?o { ?s rdf:type owl:Class . ?s ncbitaxon:exactMatch ?o }',
#         initNs=ns)
#     mondo_umls_reflist = [x for x in mondo_umls_refs if x[1].find('umls') >= 0]
#     return uris_to_ids(mondo_umls_reflist)

def test_vo(graph): 
    ns = dict(owl=OWL, rdfs=RDFS, rdf=RDF, efo=EFO, skos=SKOS, oboInOwl=OBOINOWL)
    mondo_umls_refs = graph.query(
        'SELECT ?concept ?label { ?concept rdf:type owl:Class ; oboInOwl:hasDbXref ?label }',
        initNs=ns)
    mondo_umls_reflist = [x for x in mondo_umls_refs if x[1].find('UMLS') >= 0 or x[1].find('UMLS_CUI') >= 0]
    return uris_to_ids(mondo_umls_reflist)
import glob
import json
import os

from biothings.utils.dataload import alwayslist
from biothings.utils.dataload import tabfile_feeder, dict_sweep, unlist
import mygene


def load_data(data_folder):
    # Ontology data
    go_file = os.path.join(data_folder, "go.json")
    goterms = parse_ontology(go_file)
    # Gene annotation files
    for f in glob.glob(os.path.join(data_folder, "*.gaf.gz")):
        print("Parsing {}".format(f))
        docs = parse_gene_annotations(f)

        # Create gene ID cache. Join all gene sets and fetch ids.
        all_genes = set()
        for _id, annotations in docs.items():
            for key in ["genes", "excluded_genes", "contributing_genes",
                        "colocalized_genes"]:
                if annotations.get(key) is not None:
                    all_genes = all_genes | annotations[key]
        uniprot = [i for i, j in all_genes]
        symbols = [j for i, j in all_genes]
        taxid = annotations['taxid']
        genecache = get_gene_ids(symbols, uniprot, taxid)

        for _id, annotations in docs.items():
            # Add ontology annotations
            annotations['go'] = goterms[_id]
            # Add gene sets
            if annotations.get("genes") is not None:
                genes = []
                for u, s in annotations['genes']:
                    if genecache.get(u) is not None:
                        genes += [g for g in genecache[u].values()]
                    elif genecache.get(s) is not None:
                        genes += [g for g in genecache[s].values()]
                    #else:
                    #    genes += {'symbol': s, 'uniprot': u}
                annotations['genes'] = genes
            else:
                # No genes in set
                continue

            for key in ["excluded_genes", "contributing_genes",
                        "colocalized_genes"]:
                if annotations.get(key) is not None:
                    genes = []
                    for u, s in annotations.pop(key):
                        if genecache.get(u) is not None:
                            genes += [g for g in genecache[u].values()]
                        elif genecache.get(s) is not None:
                            genes += [g for g in genecache[s].values()]
                        #else:
                        #    genes += {'symbol': s, 'uniprot': u}
                    annotations['go'][key] = genes
            # Clean up data
            annotations = unlist(annotations)
            yield annotations


def get_gene_ids(symbols, uniprot_ids, taxid):
    """Fetch NCBI, Ensembl, and gene names from UniProt ids or gene symbol."""
    mg = mygene.MyGeneInfo()
    fields = 'entrezgene,ensembl.gene,name,symbol'
    # Fetch ids from  UniProt
    response = mg.querymany(uniprot_ids, scopes='uniprot', fields=fields,
                            species=taxid, returnall=True)
    genes = {}
    for out in response['out']:
        if out.get("_id") is not None:
            query = out['query']
            geneid = out['_id']
            hits = genes.setdefault(query, {})
            hits[geneid] = {"mygene_id": geneid,
                            "uniprot": query,
                            "symbol": out.get('symbol'),
                            "name": out.get('name')}
            if out.get("entrezgene") is not None:
                hits[geneid].setdefault('ncbigene', [])
                if out['entrezgene'] not in hits[geneid]['ncbigene']:
                    hits[geneid]['ncbigene'].append(out['entrezgene'])
            if out.get("ensembl") is not None:
                hits[geneid].setdefault('ensemblgene', [])
                hits[geneid]['ensemblgene'] = hits[geneid]['ensemblgene'] + \
                    [i['gene'] for i in alwayslist(out['ensembl'])]
    # Retry missing using gene symbol
    retry = [symbols[uniprot_ids.index(k)] for k in response['missing']]
    response = mg.querymany(retry, scopes='symbol', fields=fields,
                            species=taxid, returnall=True)
    for out in response['out']:
        if out.get("_id") is not None:
            query = out['query']
            geneid = out['_id']
            hits = genes.setdefault(query, {})
            hits[geneid] = {"mygene_id": geneid,
                            "uniprot": uniprot_ids[symbols.index(query)],
                            "symbol": out['symbol'],
                            "name": out.get('name')}
            if out.get("entrezgene") is not None:
                hits[geneid].setdefault('ncbigene', [])
                if out['entrezgene'] not in hits[geneid]['ncbigene']:
                    hits[geneid]['ncbigene'].append(out['entrezgene'])
            if out.get("ensembl") is not None:
                hits[geneid].setdefault('ensemblgene', [])
                hits[geneid]['ensemblgene'] = hits[geneid]['ensemblgene'] + \
                    [i['gene'] for i in alwayslist(out['ensembl'])]
    genes = unlist(genes)
    genes = dict_sweep(genes, vals=[None, 'null'])
    return genes


def parse_gene_annotations(f):
    """Parse a gene annotation (.gaf.gz) file."""
    data = tabfile_feeder(f, header=0)
    genesets = {}
    for rec in data:
        if not rec[0].startswith("!"):
            _id = rec[4].replace(":", "_")
            if genesets.get(_id) is None:
                taxid = rec[12].split("|")[0].replace("taxon:", "")
                genesets[_id] = {"_id":  _id + "_" + taxid,
                                 "is_public": True,
                                 "taxid": taxid}
            uniprot = rec[1]
            symbol = rec[2]
            qualifiers = rec[3].split("|")
            # The gene can belong to several sets:
            if "NOT" in qualifiers:
                # Genes similar to genes in go term, but should be excluded
                genesets[_id].setdefault("excluded_genes", set()).add(
                        (uniprot, symbol))
            if "contributes_to" in qualifiers:
                # Genes that contribute to the specified go term
                genesets[_id].setdefault("contributing_genes", set()).add(
                        (uniprot, symbol))
            if "colocalizes_with" in qualifiers:
                # Genes colocalized with specified go term
                genesets[_id].setdefault("colocalized_genes", set()).add(
                        (uniprot, symbol))
            else:
                # Default set: genes that belong to go term
                genesets[_id].setdefault("genes", set()).add(
                        (uniprot, symbol))
    return genesets


def parse_ontology(f):
    "Get GO-term metadata from ontology JSON dump."
    with open(f, 'r') as infile:
        data = json.load(infile)
    nodes = data['graphs'][0]['nodes']
    go_terms = {}
    for node in nodes:
        url = node['id']
        _id = url.split("/")[-1]
        if not _id.startswith("GO_"):
            continue
        go_terms[_id] = {"id": _id.replace("_", ":"),
                         "url": url,
                         "name": node.get('lbl'),
                         "type": node.get('type')}
        if node['meta'].get("definition"):
            go_terms[_id]['definition'] = node['meta']['definition'].get('val')
            go_terms[_id]['xrefs'] = node['meta']['definition'].get('xrefs')
    go_terms = unlist(go_terms)
    go_terms = dict_sweep(go_terms)
    return go_terms


if __name__ == "__main__":

    annotations = load_data("./test_data")
    for a in annotations:
        print(json.dumps(a, indent=2))

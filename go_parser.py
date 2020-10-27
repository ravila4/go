import os

from biothings.utils.dataload import tabfile_feeder, dict_sweep, unlist
import mygene


def load_data(data_folder):
    # Ontology data
    go_file = os.path.join(data_folder, "go.obo")
    # Gene annotation files
    go_human = os.path.join(data_folder, "goa_human.gaf.gz")
    # Parse files
    goterms = parse_obo(go_file)
    docs = parse_gaf(go_human)

    # Join all gene sets and get NCBI IDs
    all_genes = set()
    for _id, annotations in docs.items():
        for key in ["genes", "go.excluded_genes", "go.contributing_genes",
                    "go.colocalized_genes"]:
            if annotations.get(key) is not None:
                all_genes = all_genes | annotations[key]
    uniprot = [i for i, j in all_genes]
    symbols = [j for i, j in all_genes]
    taxid = annotations["taxid"]
    NCBI_dict = get_NCBI_id(symbols, uniprot, taxid)

    # Add additional annotaions
    for _id, annotations in docs.items():
        # Add additional annotations
        annotations["go.name"] = goterms[_id]["name"]
        annotations["go.type"] = goterms[_id]["namespace"]
        annotations["go.description"] = goterms[_id]["def"]
        if annotations.get("genes") is not None:
            gene_dict = {}
            # Convert set of tuples to lists
            gene_dict["uniprot"] = [i for i, j in annotations["genes"]]
            gene_dict["symbols"] = [j for i, j in annotations["genes"]]
            gene_dict["NCBIgene"] = []
            for i, j in annotations["genes"]:
                if NCBI_dict.get(i):
                    gene_dict["NCBIgene"].append(NCBI_dict[i])
                elif NCBI_dict.get(j):
                    gene_dict["NCBIgene"].append(NCBI_dict[j])
            annotations["genes"] = gene_dict
        else:
            # No genes in set
            continue
        for key in ["go.excluded_genes", "go.contributing_genes",
                    "go.colocalized_genes"]:
            if annotations.get(key) is not None:
                gene_dict = {}
                gene_dict["uniprot"] = [i for i, j in annotations[key]]
                gene_dict["symbols"] = [j for i, j in annotations[key]]
                gene_dict["NCBIgene"] = []
                for i, j in annotations[key]:
                    if NCBI_dict.get(i):
                        gene_dict["NCBIgene"].append(NCBI_dict[i])
                    elif NCBI_dict.get(j):
                        gene_dict["NCBIgene"].append(NCBI_dict[j])
                annotations[key] = gene_dict
        # Clean up data
        annotations = dict_sweep(annotations)
        annotations = unlist(annotations)
        yield annotations


def get_NCBI_id(symbols, uniprot_ids, taxid):
    """Fetch NCBI id from gene symbols using mygene.info.
    If a gene symbol matches more than one NCBI id, all duplicate ids are kept.
    Gene symbols that are not found are retried using Uniprot ID.
    """
    mg = mygene.MyGeneInfo()
    response = mg.querymany(symbols, scopes='symbol', fields='entrezgene',
                            species=taxid, returnall=True)
    ncbi_ids = {}
    for out in response['out']:
        if out.get("entrezgene") is not None:
            query = out["query"]
            entrezgene = out["entrezgene"]
            ncbi_ids.setdefault(query, []).append(entrezgene)
    # Retry missing
    retry = [uniprot_ids[symbols.index(k)] for k in response["missing"]]
    response = mg.querymany(retry, scopes='uniprot', fields='entrezgene',
                            species=taxid, returnall=True)
    for out in response['out']:
        if out.get("entrezgene") is not None:
            query = out["query"]
            entrezgene = out["entrezgene"]
            ncbi_ids.setdefault(query, []).append(entrezgene)
    ncbi_ids = unlist(ncbi_ids)
    return ncbi_ids


def parse_gaf(f):
    """Parse a gene annotation (.gaf.gz) file."""
    data = tabfile_feeder(f, header=0)
    genesets = {}
    for rec in data:
        if not rec[0].startswith("!"):
            _id = rec[4]  # Primary ID is GO ID
            if genesets.get(_id) is None:
                genesets[_id] = {}  # Dict to hold annotations
                genesets[_id]["_id"] = _id
                genesets[_id]["is_public"] = True
                genesets[_id]["taxid"] = rec[12].split("|")[0].replace(
                        "taxon:", "")
            gene_id = rec[1]
            symbol = rec[2]
            qualifiers = rec[3].split("|")
            # The gene can belong to several sets:
            if "NOT" in qualifiers:
                # Genes similar to genes in go term, but should be excluded
                genesets[_id].setdefault("go.excluded_genes", set()).add(
                        (gene_id, symbol))
            if "contributes_to" in qualifiers:
                # Genes that contribute to the specified go term
                genesets[_id].setdefault("go.contributing_genes", set()).add(
                        (gene_id, symbol))
            if "colocalizes_with" in qualifiers:
                # Genes colocalized with specified go term
                genesets[_id].setdefault("go.colocalized_genes", set()).add(
                        (gene_id, symbol))
            else:
                # Default list: genes that belong to go term
                genesets[_id].setdefault("genes", set()).add(
                        (gene_id, symbol))
    return genesets


def parse_obo(f):
    """Parse gene ontology (.obo) file.
    This file contains description and metadata for individual GO terms.
    Return a dictionary of GO term metadata with GO ids as keys.
    """
    go_terms = {}
    with open(f, 'r') as data:
        termflag = False  # Keep track of beginning/end of each entry
        for line in data:
            line = line.replace("\n", "")
            if line == "":  # End of entry
                if termflag:
                    _id = data['id']
                    go_terms[_id] = data
                    alt_id = data.get("alt_id")
                    if alt_id is not None:
                        # Add a duplicate entry with alt_id as key
                        go_terms[alt_id] = data
                termflag = False
            if termflag:
                key, value = line.split(": ", 1)
                value = value.replace('"', '')
                # Remove comments
                if value.find("!") != -1:
                    value = value.split(" ! ", 1)[0]
                data[key] = value
            if line.startswith("[Term]"):  # Start of new entry
                termflag = True
                data = {}
    return go_terms


if __name__ == "__main__":
    import json

    annotations = load_data("./test_data")
    for a in annotations:
        print(json.dumps(a, indent=2))

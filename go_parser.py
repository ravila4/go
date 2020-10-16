import os
import tracemalloc

from biothings.utils.dataload import tabfile_feeder, dict_sweep, unlist


def load_annotations(data_folder):
    # Ontology data
    go_file = os.path.join(data_folder, "go.obo")
    # Annotation files
    go_human = os.path.join(data_folder, "goa_human.gaf.gz")

    tracemalloc.start()
    goterms = parse_obo(go_file)
    docs = parse_gaf(go_human)
    current, peak = tracemalloc.get_traced_memory()
    print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
    tracemalloc.stop()

    for _id, annotations in docs.items():
        # Add additional annotations
        annotations["name"] = goterms[_id]["id"]
        annotations["namespace"] = goterms[_id]["namespace"]
        annotations["def"] = goterms[_id]["def"]
        # Clean up data
        annotations = dict_sweep(annotations)
        annotations = unlist(annotations)
        yield annotations


def parse_gaf(f):
    """Parse a gene annotation (.gaf.gz) file."""
    data = tabfile_feeder(f, header=0)
    genesets = {}
    for rec in data:
        if not rec[0].startswith("!"):
            _id = rec[4]  # Primary ID is GO ID
            if genesets.get(_id) is None:
                genesets[_id] = {}
                genesets[_id]["_id"] = _id
                genesets[_id]["is_public"] = True
                genesets[_id]["genes"] = {}
                genesets[_id]["taxid"] = rec[12].split("|")[0].replace("taxon:", "")
            genes = genesets[_id]["genes"]
            gene_id = rec[1]
            gene_symbol = rec[2]
            qualifiers = rec[3].split("|")
            if "NOT" in qualifiers:
                genesets[_id].setdefault("excluded_genes", {})
                excluded = genesets[_id]["excluded_genes"]
                excluded.setdefault("uniprot", []).append(gene_id)
                excluded.setdefault("gene_symbols", []).append(gene_symbol)
            if "contributes_to" in qualifiers:
                genesets[_id].setdefault("contributing_genes", {})
                contributing = genesets[_id]["contributing_genes"]
                contributing.setdefault("uniprot", []).append(gene_id)
                contributing.setdefault("gene_symbols", []).append(gene_symbol)
            if "colocalizes_with" in qualifiers:
                genesets[_id].setdefault("colocalized_genes", {})
                colocalized = genesets[_id]["colocalized_genes"]
                colocalized.setdefault("uniprot", []).append(gene_id)
                colocalized.setdefault("gene_symbols", []).append(gene_symbol)
            else:
                # No qualifiers, add gene to default list
                genes.setdefault("uniprot", []).append(gene_id)
                genes.setdefault("gene_symbols", []).append(gene_symbol)
    return genesets


def parse_obo(f):
    """Parse gene ontology (.obo) file.
    This file contains description and metadata for individual GO terms."""
    go_terms = {}
    with open(f, 'r') as data:
        termflag = False
        for line in data:
            line = line.replace("\n", "")
            if line == "":
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
            if line.startswith("[Term]"):
                termflag = True
                data = {}
    return go_terms


if __name__ == "__main__":
    import json

    annotations = load_annotations("./")
    for a in annotations:
        #pass
        print(json.dumps(a, indent=2))

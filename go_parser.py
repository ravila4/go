import os

from biothings.utils.dataload import tabfile_feeder, dict_sweep, unlist


def load_data(data_folder):
    # Ontology data
    go_file = os.path.join(data_folder, "go.obo")
    # Annotation files
    go_human = os.path.join(data_folder, "goa_human.gaf.gz")
    # Parse files
    goterms = parse_obo(go_file)
    docs = parse_gaf(go_human)

    for _id, annotations in docs.items():
        # Add additional annotations
        annotations["name"] = goterms[_id]["id"]
        annotations["namespace"] = goterms[_id]["namespace"]
        annotations["def"] = goterms[_id]["def"]
        # Convert gene sets to lists
        if annotations.get("genes") is not None:
            annotations["genes"]["uniprot"] = list(annotations["genes"]["uniprot"])
            annotations["genes"]["gene_symbols"] = list(annotations["genes"]["gene_symbols"])
            # Drop go term if set has fewer genes than threshold
            if len(annotations["genes"]["uniprot"]) < 10:
                continue
        else:
            # No genes in list
            continue
        for key in ["excluded_genes", "contributing_genes", "colocalized_genes"]:
            if annotations.get(key) is not None:
                annotations[key]["uniprot"] = list(annotations[key]["uniprot"])
                annotations[key]["gene_symbols"] = list(annotations[key]["gene_symbols"])
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
                genesets[_id] = {}  # Dict to hold annotations
                genesets[_id]["_id"] = _id
                genesets[_id]["is_public"] = True
                genesets[_id]["taxid"] = rec[12].split("|")[0].replace("taxon:", "")
            gene_id = rec[1]
            gene_symbol = rec[2]
            qualifiers = rec[3].split("|")
            # The gene can belong to several lists:
            if "NOT" in qualifiers:
                # Genes that are similar to genes in go term, but should be excluded
                genesets[_id].setdefault("excluded_genes", {})
                excluded = genesets[_id]["excluded_genes"]
                excluded.setdefault("uniprot", set()).add(gene_id)
                excluded.setdefault("gene_symbols", set()).add(gene_symbol)
            if "contributes_to" in qualifiers:
                # Genes that contribute to the specified go term
                genesets[_id].setdefault("contributing_genes", {})
                contributing = genesets[_id]["contributing_genes"]
                contributing.setdefault("uniprot", set()).add(gene_id)
                contributing.setdefault("gene_symbols", set()).add(gene_symbol)
            if "colocalizes_with" in qualifiers:
                # Genes colocalized with specified go term
                genesets[_id].setdefault("colocalized_genes", {})
                colocalized = genesets[_id]["colocalized_genes"]
                colocalized.setdefault("uniprot", set()).add(gene_id)
                colocalized.setdefault("gene_symbols", set()).add(gene_id)
            else:
                # Default list: genes that belong to go term
                genesets[_id].setdefault("genes", {})
                genes = genesets[_id]["genes"]
                genes.setdefault("uniprot", set()).add(gene_id)
                genes.setdefault("gene_symbols", set()).add(gene_symbol)
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

    annotations = load_data("./")
    for a in annotations:
        print(json.dumps(a, indent=2))

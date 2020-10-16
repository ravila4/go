import os
import requests
import tracemalloc
from biothings.utils.dataload import tabfile_feeder, dict_sweep, unlist

def load_annotations(data_folder):
    # Annotation files
    go_human = os.path.join(data_folder, "goa_human.gaf.gz")

    def parse_gaf(f, species):
        """Parse a .gaf.gz file"""
        data = tabfile_feeder(f, header=0)
        genesets = {}
        for rec in data:
            if not rec[0].startswith("!"):
                _id = rec[4]  # Primary ID is GO ID
                if genesets.get(_id) is None:
                    genesets[_id] = {}
                    genesets[_id]["_id"] = _id
                    genesets[_id]["species"] = species
                    genesets[_id]["is_public"] = True
                    genesets[_id]["genes"] = {}
                genes = genesets[_id]["genes"]
                gene_id = rec[1]
                gene_symbol = rec[2]
                qualifiers = rec[3].split("|")
                if "NOT" in qualifiers:
                    genesets[_id]["excluded_genes"] = {}
                    excluded = genesets[_id]["excluded_genes"]
                    excluded.setdefault("uniprot", []).append(gene_id)
                    excluded.setdefault("gene_symbols", []).append(gene_symbol)
                if "contributes_to" in qualifiers:
                    genesets[_id]["contributing_genes"] = {}
                    contributing = genesets[_id]["contributing_genes"]
                    contributing.setdefault("uniprot", []).append(gene_id)
                    contributing.setdefault("gene_symbols", []).append(gene_symbol)
                if "colocalizes_with" in qualifiers:
                    genesets[_id]["colocalized_genes"] = {}
                    colocalized = genesets[_id]["colocalized_genes"]
                    colocalized.setdefault("uniprot", []).append(gene_id)
                    colocalized.setdefault("gene_symbols", []).append(gene_symbol)
                else:
                    # No qualifiers, add gene to default list
                    genes.setdefault("uniprot", []).append(gene_id)
                    genes.setdefault("gene_symbols", []).append(gene_symbol)
        return genesets

    tracemalloc.start()
    docs = parse_gaf(go_human, "human")
    current, peak = tracemalloc.get_traced_memory()
    print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
    tracemalloc.stop()

    for _id, annotations in docs.items():
        annotations = dict_sweep(annotations)
        annotations = unlist(annotations)
        yield annotations


if __name__ == "__main__":
    import json

    annotations = load_annotations("./")
    for a in annotations:
        pass
        print(json.dumps(a, indent=2))

import requests

def get_release(self):
    release_summary = requests.get("http://current.geneontology.org/summary.txt")
    for line in release_summary.text.split("\n"):
        if line.startswith("Start date:"):
            version = line.split(": ")[-1]
    return version

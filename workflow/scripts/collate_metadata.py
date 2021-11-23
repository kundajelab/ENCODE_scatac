import json
import csv

def collate(paths, out):
    entries = []
    for p in paths:
        with open(p, "r") as f:
            j = json.load(f)
            entries.append(j)

    cols = set()
    for e in entries:
        cols |= e.keys()
    
    cols = list(cols)
    lines = [cols]
    for j in entries:
        line = [json.dumps(j.get(col)) for col in cols]
        lines.append(line)

    with open(out, 'w', newline='') as outf:
        writer = csv.writer(outf, delimiter="\t", lineterminator='\n', quotechar="'")
        writer.writerows(lines)

if __name__ == '__main__':
    try:
        paths = snakemake.input
        out, = snakemake.output
        collate(paths, out)

    except NameError:
        pass


import json
import csv

def collate(paths, out):
    lines = []
    for i, p in enumerate(paths):
        with open(p, "r") as f:
            j = json.load(f)
            if i == 0:
                header = list(j.keys())
                lines.append(header)
            line = [json.dumps(j[col]) for col in header]
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


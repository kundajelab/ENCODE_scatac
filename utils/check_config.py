import sys
from encode_utils.connection import Connection

def fetch_info(experiment, conn):
    data = conn.get(experiment)

    lab = data["lab"]["name"]
    multiome = data["related_series"].get("accession")
    
    return lab, multiome

def load_samples(sample_path):
    with open(sample_path) as sample_file:
        h = sample_file.readline().rstrip('\n').split('\t')
        exp_ind = h.index("Experiment")
        rep_ind = h.index("Replicate")
        mod_ind = h.index("Modality")
        gen_ind = h.index("Genome")
        samples = []
        for line in sample_file:
            if line.startswith("#"):
                continue
            if line.startswith("@"):
                continue
            if line.startswith("$"):
                continue
            entries = line.rstrip('\n').split('\t')
            exp = entries[exp_ind]
            rep = int(entries[rep_ind])
            mod = entries[mod_ind]
            gen = entries[gen_ind]
            data = [
                exp,
                rep,
                mod,
                gen
            ]
            samples.append(data)
    return samples

def check_config(sample_path, out_path, dcc_mode):
    conn = Connection(dcc_mode)
    samples = load_samples(sample_path)
    with open(out_path, "w") as out_file:
        for s in samples:
            lab, mult = fetch_info(s[0], conn)
            print(*(s + [lab, mult]), sep="\t", file=out_file)

if __name__ == '__main__':
    sample_path = sys.argv[1]
    out_path = sys.argv[2]
    dcc_mode = "prod"
    check_config(sample_path, out_path, dcc_mode)

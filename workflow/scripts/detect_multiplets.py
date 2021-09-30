import time
import itertools
import gzip
import heapq
from collections import Counter

import numpy as np
import matplotlib.pyplot as plt
# from kneed import KneeLocator

def print_and_log(text, outfile, starttime=0):
    logtime = time.process_time() - starttime
    if logtime < 60:
        logtime = "{:,}s".format(logtime)
    else:
        logtime = "{:,}m {:,}s".format(logtime // 60, logtime % 60)
    outfile.write("{} - {}\n".format(logtime, text))
    print("{} - {}".format(logtime, text))

# def multiplet_fdr(samples, nulls, null_zeros, num_bc, fdr_thresh):
#     null_total = nulls.shape[0] + null_zeros
#     sample_total = samples.shape[0]

#     p = 1 - (np.searchsorted(nulls, samples) + null_zeros) / null_total
#     # print(p)
#     # p_bonf = p * (num_bc - 1)
#     p_sidak = 1 - (1 - p)**(num_bc - 1)
#     print(p_sidak)
#     q = (p_sidak * sample_total) / (sample_total - np.arange(sample_total))
#     print(q) ####
#     candidiates = np.nonzero(q <= fdr_thresh)[0]
#     if candidiates.size == 0:
#         cut = samples[-1]
#     else:
#         cut = samples[candidiates[0]]

#     return cut, q

def multiplet_fdr(samples, nulls, fdr_thresh):
    null_total = nulls.shape[0]
    sample_total = samples.shape[0]

    p = 1 - np.searchsorted(nulls, samples) / null_total
    print(p)
    q = (p * sample_total) / (sample_total - np.arange(sample_total))
    print(q) ####
    candidiates = np.nonzero(q <= fdr_thresh)[0]
    if candidiates.size == 0:
        cut = samples[-1]
    else:
        cut = samples[candidiates[0]]

    return cut, q

def plot_dist(cut, q, samples, nulls, title, x_label, out_path, log_x=False, hist_bins=200):
    fig, ax = plt.subplots(tight_layout=True)
    ax2 = ax.twinx()
    ax2.set_yscale('log')
    if log_x:
        ax.set_xscale('log')
        hist_bins = np.geomspace(min(samples[0], nulls[0]), max(samples[-1], nulls[-1]), hist_bins)

    ax.hist(nulls, bins=hist_bins, density=True, alpha=0.5, color="k")
    ax2.plot(samples, q, color="g")
    ax.axvline(x=cut, color="r")
    ax.hist(samples, bins=hist_bins, density=True, alpha=0.5, color="b")

    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel("Histogram Frequency")
    ax2.set_ylabel("Unadjusted Q-Value")

    plt.savefig(out_path)

# def tail_cut(samples, side, min_keep=0.6):
#     rev = True if side == 'l' else False
#     if rev:
#         samples = samples[::-1] # reverse order if detecting left tail

#     lbound = int(samples.shape[0] * min_keep)
#     total_mean = samples.mean()
#     s = samples - total_mean # shift origin for numerical stability

#     m0 = np.arange(1, s.shape[0] + 1) # Cumulative sample count
#     m1 = np.cumsum(s) / m0 # 1st cumulative moment (cumulative mean)
#     o2 = np.cumsum(s**2) / m0 # 2nd cumulant moment around origin
#     m2 = o2 - m1**2 # 2nd cumulative central moment (cumulative variance)
#     o3 = np.cumsum(s**3) / m0 # 3rd cumulative moment around origin
#     o4 = np.cumsum(s**4) / m0 # 4th cumulative moment around origin
#     m4 = o4 - 4 * m1 * o3 + 6 * m1**2 * o2 - 3 * m1**4 # 4th cumulative central moment

#     with np.errstate(divide='ignore', invalid='ignore'):
#         k = m4 / m2**2 # Cumulative kurtosis

#     cut_ind = np.nanargmin(k[lbound:])
#     cut_k = k[lbound:][cut_ind]
#     cut = samples[lbound:][cut_ind]
#     bound = samples[lbound]
#     lower_visual = min(max(100, int(samples.shape[0] * 0.2)), lbound)
#     k[:lower_visual] = np.nan
    
#     if rev:
#         k = k[::-1] # restore original order
#         cut_ind = samples.shape[0] - cut_ind - 1

#     return cut_ind, cut_k, cut, bound, k

# def plot_cut(cut, k, pts, lb, title, x_label, out_path, log_scale=False, hist_bins=200):
#     fig, ax = plt.subplots(tight_layout=True)
#     if log_scale:
#         ax.set_xscale('log')
#         hist_bins = np.geomspace(pts.min(), pts.max(), hist_bins)

#     ax.hist(pts, bins=hist_bins)
#     ax2 = ax.twinx()
#     ax2.plot(pts, k, color="g")
#     ax.axvline(x=cut, color="r")
#     ax.axvline(x=lb, color="k")

#     ax.set_title(title)
#     ax.set_xlabel(x_label)
#     ax.set_ylabel("Histogram Frequency")
#     ax2.set_ylabel("Cumulative Kurtosis")

#     plt.savefig(out_path)

def main(fragments, barcodes_strict, barcodes_expanded, summary, barcodes_status, jac_plot, min_counts=500, max_frag_clique=6, fdr_thresh=0.2):
    logout = open(summary, "w")
    starttime = time.process_time() 

    print_and_log("Identifying candidate barcodes", logout, starttime)

    cur_clique = set()
    cur_coord = None
    i = 0
    barcode_counts = Counter()
    with gzip.open(fragments, 'rt') as f:
        for line in f:
            line = line.rstrip('\n').split('\t')
            chr, start, end, barcode = line[:4]

            this_coord = (chr, start, end) 
            if this_coord != cur_coord:
                if len(cur_clique) <= max_frag_clique:
                    for b in cur_clique:
                        barcode_counts[b] += 1

                cur_clique = set([barcode])
                cur_coord = this_coord

            else:
                cur_clique.add(barcode)

            i += 1
            if i%1e7==0:
                print(i)


    print_and_log(
        f"Original run had {len(barcode_counts)} total cell barcodes",
        logout,
        starttime,
    )
    
    barcodes_considered = set(k for k, v in barcode_counts.items() if v >= min_counts)
    num_bc = len(barcodes_considered)

    print_and_log(
        f"Identified {num_bc} total barcodes for multiplet detection",
        logout,
        starttime,
    )

    print_and_log("Reading fragments", logout, starttime)

    pair_counts = Counter()
    cur_clique = set()
    cur_coord = None
    i = 0
    with gzip.open(fragments, 'rt') as f:
        for line in f:
            line = line.rstrip('\n').split('\t')
            chr, start, end, barcode = line[:4]

            if barcode not in barcodes_considered:
                continue

            this_coord = (chr, start, end) 
            if this_coord != cur_coord:
                if len(cur_clique) <= max_frag_clique:
                    for x, y in itertools.combinations(cur_clique, 2):
                        x, y = (x, y) if x < y else (y, x)
                        pair_counts[(x, y)] += 1

                cur_clique = set([barcode])
                cur_coord = this_coord

            else:
                cur_clique.add(barcode)

            i += 1
            if i%1e7==0:
                print(i)

    print_and_log("Identifying barcode multiplets", logout, starttime)

    print_and_log(
        f"Considered {len(pair_counts)} barcode pairs",
        logout,
        starttime,
    )

    expanded_data = {}
    jac_dists_max = {}
    # jac_dists_pairs = {}
    jac_dists_top = {}
    for x, y in pair_counts.items():
        a, b = x
        bca = barcode_counts[a]
        bcb = barcode_counts[b]
        jac = y/(bca + bcb - y)
        data = [a, b, bca, bcb, y, jac, None]
        expanded_data[x] = data
        if jac > 0:
            # jac_dists_pairs[x] = jac 
            jac_dists_max[a] = max(jac_dists_max.get(a, 0), jac)
            jac_dists_max[b] = max(jac_dists_max.get(b, 0), jac)

            aheap = jac_dists_top.setdefault(a, [])
            if len(aheap) < 7:
                heapq.heappush(aheap, jac)
            else:
                heapq.heappushpop(aheap, jac)

            bheap = jac_dists_top.setdefault(b, [])
            if len(aheap) < 7:
                heapq.heappush(bheap, jac)
            else:
                heapq.heappushpop(bheap, jac)

    jac_dists_7th = {k: (min(v) if len(v) == 7 else 0) for k, v in jac_dists_top.items()}

    with gzip.open(barcodes_expanded, 'wt') as f:
        f.write("Barcode1\tBarcode2\tBarcode1Counts\tBarcode2Counts\tCommon\tJaccardIndex\n")
        for x, data in expanded_data.items():
            a, b = x
            f.write("{}\t{}\t{}\t{}\t{}\t{:.4f}\n".format(*data[:-1]))

    samples = np.fromiter(jac_dists_max.values(), dtype=float, count=len(jac_dists_max))
    samples.sort()
    nulls = np.fromiter(jac_dists_7th.values(), dtype=float, count=len(jac_dists_7th))
    # nulls = np.fromiter(jac_dists_pairs.values(), dtype=float, count=len(jac_dists_pairs))
    nulls.sort()
    null_zeros = num_bc * (num_bc - 1) / 2 - len(jac_dists_pairs)

    # cut, q = multiplet_fdr(samples, nulls, null_zeros, num_bc, fdr_thresh)
    cut, q = multiplet_fdr(samples, nulls, fdr_thresh)
    plot_dist(cut, q, samples, nulls, "Multiplet Thresholding", "Max Marginal Jaccard Distance", jac_plot, log_x=True)
    
    # cut_ind_jac, cut_k_jac, cut_jac, bound_jac, k_jac = tail_cut(dist_jac, 'r')
    # plot_cut(cut_jac, k_jac, dist_jac, bound_jac, "Multiplet Thresholding", "Max Marginal Jaccard Distance", jac_plot, log_scale=True)

    min_jac = cut

    print_and_log(
        f"Setting multiplet threshold as {min_jac} for minimum pairwise Jaccard distance",
        logout,
        starttime,
    )

    multiplet_data = {}
    primary_bc_map = {}
    bc_sets = {}
    for x, y in expanded_data.items():
        jac = y[5]
        if jac >= min_jac: 
            multiplet_data[x] = y
            a, b = x

            a_primary = primary_bc_map.setdefault(a, a)
            b_primary = primary_bc_map.setdefault(b, b)
            a_set = bc_sets.setdefault(a_primary, set([a])) # initialize starting sets if needed
            b_set = bc_sets.setdefault(b_primary, set([b]))
            set_info = {a_primary: a_set, b_primary: b_set}
            
            if a_primary != b_primary:
                remaining_primary = max([a_primary, b_primary], key=barcode_counts.get) 
                other_primary = a_primary if remaining_primary == b_primary else b_primary   
                for k in set_info[other_primary]:
                    primary_bc_map[k] = remaining_primary

                a_set |= b_set
                bc_sets[remaining_primary] = a_set
                bc_sets.pop(other_primary)
                
    blacklist = set()
    with open(barcodes_strict, 'w') as f:
        f.write("Barcode1\tBarcode2\tBarcode1Counts\tBarcode2Counts\tCommon\tJaccardIndex\tPrimaryBarcode\n")
        for x, data in multiplet_data.items():
            a, b = x
            pb = primary_bc_map[a]
            if a != pb:
                blacklist.add(a)
            if b != pb:
                blacklist.add(b)
            data[-1] = pb
            f.write("{}\t{}\t{}\t{}\t{}\t{:.4f}\t{}\n".format(*data))

    with open(barcodes_status, 'w') as f:
        f.write("Barcode\tIsMultiplet\tPrimaryBarcode\n")
        for b in barcode_counts.keys():
            if b not in barcodes_considered:
                f.write(f"{b}\tIndeterminate\tNone\n")
            elif b in multiplet_data:
                f.write(f"{b}\tTrue\t{primary_bc_map[b]}\n")
            else:
                f.write(f"{b}\tFalse\tNone\n")

    print_and_log(
        f"Identified {len(multiplet_data)} barcode pairs above Jaccard threshold",
        logout,
        starttime,
    )

    print_and_log(
        f"Identified {len(primary_bc_map)} barcodes belonging to multiplets",
        logout,
        starttime,
    )

    print_and_log(
        f"Identified {len(bc_sets)} unique multiplets",
        logout,
        starttime,
    )

    print_and_log(
        f"After multiplet exclusions, have {len(barcode_counts) - len(blacklist)} total cell barcodes",
        logout,
        starttime,
    )

    logout.close()

if __name__ == '__main__':
    try:
        fragments = snakemake.input['frag']

        barcodes_strict = snakemake.output['barcode_pairs_strict']
        barcodes_expanded = snakemake.output['barcode_pairs_expanded']
        barcodes_status = snakemake.output['barcodes_status']
        summary = snakemake.output['qc']
        jac_plot = snakemake.output['multiplets_thresh']

        main(fragments, barcodes_strict, barcodes_expanded, summary, barcodes_status, jac_plot)

    except NameError:
        main('/dev/stdin', '/dev/stdout', '/dev/null', '/dev/null', '/dev/null')
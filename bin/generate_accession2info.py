#!/usr/bin/env python3
import sys

RANKS = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

def load_lengths(path):
    result = {}
    with open(path) as f:
        for line in f:
            acc, length = line.strip().split('\t')
            result[acc] = int(length)
    return result

def load_filtered_acc2taxid(path, required_accs):
    result = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split('\t', 2)
            acc = parts[0]
            if acc in required_accs:
                taxid = parts[1]
                result[acc] = taxid
                if len(result) == len(required_accs):
                    break
    return result

def collect_needed_taxids(acc2taxid, nodes_path):
    needed_taxids = set(acc2taxid.values())
    parent_map = {}
    with open(nodes_path) as f:
        for line in f:
            parts = [x.strip() for x in line.split('|')]
            tid, ptid = parts[0], parts[1]
            parent_map[tid] = ptid

    expanded = set()
    for tid in needed_taxids:
        cur = tid
        while True:
            expanded.add(cur)
            if cur not in parent_map or parent_map[cur] in (cur, '1'):
                break
            cur = parent_map[cur]
    return expanded

def load_nodes_filtered(path, taxids_to_keep):
    parent, rank = {}, {}
    with open(path) as f:
        for line in f:
            parts = [x.strip() for x in line.split('|')]
            tid = parts[0]
            if tid in taxids_to_keep:
                parent[tid] = parts[1]
                rank[tid] = parts[2]
    return parent, rank

def load_names_filtered(path, taxids_to_keep):
    names = {}
    with open(path) as f:
        for line in f:
            parts = [x.strip() for x in line.split('|')]
            tid = parts[0]
            if tid in taxids_to_keep and parts[3] == 'scientific name':
                names[tid] = parts[1]
    return names

def get_structured_lineage(taxid, parent_map, rank_map, names_map):
    lineage_taxids = []
    cur = taxid
    while True:
        lineage_taxids.append(cur)
        if cur not in parent_map or parent_map[cur] in (cur, '1'):
            break
        cur = parent_map[cur]
    lineage_taxids.reverse()

    lineage_taxids = [
        tid for tid in lineage_taxids
        if names_map.get(tid, "").lower() not in ("root", "cellular organisms")
    ]

    rank_to_name = {}
    for tid in lineage_taxids:
        rk = rank_map.get(tid)
        if rk in RANKS and rk not in rank_to_name:
            rank_to_name[rk] = names_map.get(tid, tid)

    return [rank_to_name.get(r, 'NA') for r in RANKS]

def main():
    if len(sys.argv) != 6:
        sys.exit(1)

    lengths_file = sys.argv[1]
    acc2taxid_file = sys.argv[2]
    nodes_file = sys.argv[3]
    names_file = sys.argv[4]
    output_file = sys.argv[5]

    lengths = load_lengths(lengths_file)
    accs_needed = set(lengths.keys())
    acc2taxid = load_filtered_acc2taxid(acc2taxid_file, accs_needed)

    # Paso 1: obtener todos los taxIDs necesarios para linaje
    taxids_required = collect_needed_taxids(acc2taxid, nodes_file)

    # Paso 2: cargar solo nodos y nombres requeridos
    parent_map, rank_map = load_nodes_filtered(nodes_file, taxids_required)
    names_map = load_names_filtered(names_file, taxids_required)

    with open(output_file, 'w') as out:
        for acc, length in lengths.items():
            if acc not in acc2taxid:
                continue
            taxid = acc2taxid[acc]
            lineage = get_structured_lineage(taxid, parent_map, rank_map, names_map)
            out.write(f"{acc}\t{length}\t{taxid}\t{'|'.join(lineage)}\n")

if __name__ == '__main__':
    main()

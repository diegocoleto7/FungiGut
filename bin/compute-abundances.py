#!/usr/bin/env python3

import argparse, os, random, sys

RANKS = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

def parseargs():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('sam')
    parser.add_argument('--accession_info', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--min_map', type=int, required=True)
    parser.add_argument('--max_ed', type=int, required=True)
    parser.add_argument('--pct_id', type=float, required=True)
    parser.add_argument('--read_cutoff', type=int, required=True)
    parser.add_argument('--raw_counts', action='store_true')
    return parser.parse_args()

def find_taxid(tag):
    if '|' not in tag:
        return tag
    for sp in tag.split('|'):
        if '_' in sp:
            return sp
    return tag

def load_accession_info(infofile):
    acc2info, clade2gi, lin2len = {}, {}, {}
    with open(infofile, 'r') as infile:
        for line in infile:
            splits = line.strip().split('\t')
            if len(splits) == 4:
                acc2info[splits[0]] = [splits[1], splits[2], splits[3]]
                lin2len[splits[3]] = lin2len.get(splits[3], 0.0) + float(splits[1])
            else:
                clade2gi[splits[1]] = splits[0]
    return acc2info, clade2gi, lin2len

def filter_line(args, splits):
    cigar = splits[5]
    if cigar == '*':
        return True
    matched_len, total_len, cur = 0, 0, 0
    for ch in cigar:
        if ch.isdigit():
            cur = cur * 10 + int(ch)
        else:
            if ch in 'M=':
                matched_len += cur
            total_len += cur
            cur = 0
    if matched_len < args.min_map:
        return True
    edit_distance = int(splits[11][5:])
    if edit_distance > args.max_ed:
        return True
    if float(matched_len) / float(total_len) < args.pct_id:
        return True
    return False

def compute_abundances(args, samfile, acc2info, clade2gi, lin2len):
    with open(samfile, 'r') as infile:
        ids2abs, clade2abs, clade2ids = {}, {}, {}
        multimapped = {}
        prev_read, prev_tag, prev_count, ignore = '', '', 1.0, False

        for line in infile:
            if line.startswith('@'):
                continue
            splits = line.strip().split('\t')
            if filter_line(args, splits):
                continue
            tag = find_taxid(splits[2])
            if tag == '*':
                continue
            read_id = splits[0]
            if read_id == prev_read and tag != prev_tag:
                ignore = True
                multimapped.setdefault(read_id, []).append(prev_tag)
                prev_tag = tag
            else:
                if prev_read and not ignore:
                    ids2abs[prev_tag] = ids2abs.get(prev_tag, 0.0) + prev_count
                elif ignore:
                    multimapped.setdefault(prev_read, []).append(prev_tag)
                prev_read, prev_tag, prev_count, ignore = read_id, tag, 1.0, False

        if prev_read and not ignore:
            ids2abs[prev_tag] = ids2abs.get(prev_tag, 0.0) + prev_count
        elif ignore:
            multimapped.setdefault(prev_read, []).append(prev_tag)

    for taxid in list(ids2abs):
        if taxid not in acc2info:
            ids2abs.pop(taxid, None)
            continue
        clade = acc2info[taxid][2]
        clade2abs[clade] = clade2abs.get(clade, 0.0) + ids2abs[taxid]
        clade2ids.setdefault(clade, []).append(taxid)

    for clade in list(clade2abs):
        if clade2abs[clade] < args.read_cutoff:
            for tid in clade2ids.get(clade, []):
                ids2abs.pop(tid, None)
            clade2abs.pop(clade)
        elif not args.raw_counts:
            clade2abs[clade] /= lin2len.get(clade, 1.0)

    added = {}
    for read_id, tags in multimapped.items():
        options = {t: ids2abs[t] for t in tags if t in ids2abs and t in acc2info}
        total = sum(options.values())
        if not total:
            continue
        randnum = random.random()
        for t, ab in options.items():
            frac = ab / total
            if frac >= randnum:
                added[t] = added.get(t, 0.0) + 1.0
                break
            randnum -= frac

    for key, val in added.items():
        if key not in acc2info:
            continue
        clade = acc2info[key][2]
        if args.raw_counts:
            clade2abs[clade] += val
        else:
            clade2abs[clade] += val / lin2len.get(clade, 1.0)

    if not args.raw_counts:
        total_ab = sum(clade2abs.values())
        if total_ab:
            for clade in clade2abs:
                clade2abs[clade] *= 100.0 / total_ab

    results = {}
    for clade, val in clade2abs.items():
        taxa = clade.split('|')
        snames = [clade2gi.get(t, '') for t in taxa]
        results[clade] = [snames[-1], 'species', '|'.join(snames), clade, val]
        for i in range(len(taxa) - 1):
            level = '|'.join(taxa[:i+1])
            if level in results:
                results[level][-1] += val
            else:
                results[level] = [snames[i], RANKS[i], '|'.join(snames[:i+1]), level, val]
    return results

def main():
    args = parseargs()
    if not 0.0 <= args.pct_id <= 1.0:
        sys.exit('Error: --pct_id must be between 0.0 and 1.0')
    samfiles = [args.sam] if args.sam.endswith('.sam') else [l.strip() for l in open(args.sam)]
    acc2info, clade2gi, lin2len = load_accession_info(args.accession_info)
    combined = {}
    for sam in samfiles:
        res = compute_abundances(args, sam, acc2info, clade2gi, lin2len)
        for clade, values in res.items():
            combined[clade] = combined.get(clade, values)
            if clade in combined:
                combined[clade][-1] += values[-1]
    lev_res = {i: [] for i in range(len(RANKS))}
    for clade, vals in combined.items():
        vals[-1] /= len(samfiles)
        lev_res[len(clade.split('|')) - 1].append(vals)
    with open(args.output, 'w') as out:
        header = 'READ_COUNT' if args.raw_counts else 'PERCENTAGE'
        out.write(f'@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\t{header}\n')
        for i in range(len(RANKS)):
            for line in sorted(lev_res[i], key=lambda x: -x[-1]):
                if args.raw_counts:
                    line[-1] = int(line[-1])
                out.write('\t'.join(map(str, line)) + '\n')

if __name__ == '__main__':
    main()

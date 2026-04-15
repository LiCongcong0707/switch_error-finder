# detect_switch_errors_hexaploid_v2.py
import pysam, sys, re
from collections import defaultdict

bam_file = sys.argv[1]
out_tsv  = sys.argv[2]

# ---- 可调参数 ----
MAPQ_MIN       = 20   # 最低比对质量
MIN_CONSECUTIVE = 4   # 至少连续N个窗口在新单倍型，才算真正switch
# -----------------

def get_haplotype(contig_name):
    m = re.search(r'([A-F])$', contig_name)
    return m.group(1) if m else 'unknown'

def get_chrom_base(contig_name):
    m = re.match(r'(Chr\d+)[A-F]$', contig_name)
    return m.group(1) if m else contig_name

def parse_window_name(qname):
    parts = qname.rsplit('__win_', 1)
    read_id   = parts[0]
    win_start = int(parts[1].split('_')[0])
    return read_id, win_start

def find_real_switches(windows, min_consecutive=MIN_CONSECUTIVE):
    """
    只有连续 >= min_consecutive 个窗口都在新单倍型，才算一次真正的 switch。
    孤立的单个/少数窗口跳变视为同源区域噪音，忽略。
    """
    haps = [w['hap'] for w in windows]
    n    = len(haps)
    if n < 2:
        return []

    # 先做 majority-vote 平滑：用前后各1个窗口投票修正孤立噪音
    smoothed = list(haps)
    for i in range(1, n - 1):
        neighbors = [haps[i-1], haps[i+1]]
        if haps[i] not in neighbors and neighbors[0] == neighbors[1]:
            smoothed[i] = neighbors[0]  # 孤立跳变，修正回来

    # 在平滑后的序列上找连续 switch
    switch_events = []
    i = 1
    while i < n:
        if smoothed[i] != smoothed[i-1]:
            # 向后看，是否有 >= min_consecutive 个窗口维持新单倍型
            new_hap = smoothed[i]
            run_len = 0
            for j in range(i, min(i + min_consecutive, n)):
                if smoothed[j] == new_hap:
                    run_len += 1
                else:
                    break
            if run_len >= min_consecutive:
                switch_events.append({
                    'win_idx':     i,
                    'win_start':   windows[i]['win_start'],
                    'from_hap':    smoothed[i-1],
                    'to_hap':      new_hap,
                    'from_contig': windows[i-1]['contig'],
                    'to_contig':   windows[i]['contig'],
                })
        i += 1

    return switch_events

# ---- 收集比对结果 ----
read_windows = defaultdict(list)

with pysam.AlignmentFile(bam_file, 'rb') as bam:
    for aln in bam.fetch():
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        if aln.mapping_quality < MAPQ_MIN:
            continue
        read_id, win_start = parse_window_name(aln.query_name)
        read_windows[read_id].append({
            'win_start': win_start,
            'hap':       get_haplotype(aln.reference_name),
            'chrom':     get_chrom_base(aln.reference_name),
            'contig':    aln.reference_name,
            'mapq':      aln.mapping_quality,
        })

# ---- 判断 switch error ----
results = []

for read_id, windows in read_windows.items():
    windows.sort(key=lambda x: x['win_start'])
    n_windows = len(windows)
    if n_windows < 2:
        continue

    haps         = [w['hap'] for w in windows]
    unique_haps  = set(haps)
    switch_events = find_real_switches(windows)
    n_switches    = len(switch_events)

    if len(unique_haps) == 1:
        status = 'single_hap'
    elif n_switches == 0:
        status = 'cross_contig_correct'  # 有杂散窗口但不连续，视为噪音
    else:
        status = 'switch_error'

    results.append({
        'read_id':       read_id,
        'n_windows':     n_windows,
        'status':        status,
        'n_switches':    n_switches,
        'hap_sequence':  ','.join(haps),
        'switch_detail': ';'.join(
            f"{e['win_start']}:{e['from_hap']}->{e['to_hap']}({e['from_contig']}->{e['to_contig']})"
            for e in switch_events
        ),
    })

# ---- 写出 TSV ----
with open(out_tsv, 'w') as f:
    f.write('read_id\tn_windows\tstatus\tn_switches\thap_sequence\tswitch_detail\n')
    for r in results:
        f.write(f"{r['read_id']}\t{r['n_windows']}\t{r['status']}\t"
                f"{r['n_switches']}\t{r['hap_sequence']}\t{r['switch_detail']}\n")

# ---- 汇总统计 ----
total        = len(results)
single_hap   = sum(1 for r in results if r['status'] == 'single_hap')
cross_ok     = sum(1 for r in results if r['status'] == 'cross_contig_correct')
switch_err   = sum(1 for r in results if r['status'] == 'switch_error')
cross_total  = cross_ok + switch_err
total_sw_evt = sum(r['n_switches'] for r in results)
total_bp     = sum(r['n_windows'] for r in results) * 10000

print(f"MIN_CONSECUTIVE filter : {MIN_CONSECUTIVE} windows")
print(f"{'Total reads':<35}: {total:>7,}")
print(f"{'Single haplotype reads':<35}: {single_hap:>7,}")
print(f"{'Cross-contig correct (noise)':<35}: {cross_ok:>7,}")
print(f"{'Switch error reads':<35}: {switch_err:>7,}")
print(f"{'Total switch events':<35}: {total_sw_evt:>7,}")
if cross_total > 0:
    print(f"{'Switch error rate (cross-contig)':<35}: {switch_err/cross_total*100:>6.2f}%")
if total_bp > 0:
    print(f"{'Errors / 100 Mb':<35}: {total_sw_evt/total_bp*1e8:>6.2f}")
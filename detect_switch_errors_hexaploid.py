# detect_switch_errors_hexaploid.py
import pysam, sys, re
from collections import defaultdict

bam_file = sys.argv[1]
out_tsv  = sys.argv[2]

def get_haplotype(contig_name):
    m = re.search(r'([A-F])$', contig_name)
    return m.group(1) if m else 'unknown'

def get_chrom_base(contig_name):
    """Chr01A -> Chr01，用于判断同源染色体"""
    m = re.match(r'(Chr\d+)[A-F]$', contig_name)
    return m.group(1) if m else contig_name

def is_switch(win_prev, win_curr):
    """
    判断两个相邻窗口之间是否真正发生了 switch error
    同一单倍型内部（不管跨没跨 contig）都不算 switch
    """
    hap_prev = win_prev['hap']
    hap_curr = win_curr['hap']

    # 单倍型字母相同 -> 无 switch
    if hap_prev == hap_curr:
        return False

    # 单倍型不同 -> switch error
    return True

# ---- 收集每条 read 的窗口比对 ----
def parse_window_name(qname):
    parts = qname.rsplit('__win_', 1)
    read_id  = parts[0]
    win_start = int(parts[1].split('_')[0])
    return read_id, win_start

read_windows = defaultdict(list)

with pysam.AlignmentFile(bam_file, 'rb') as bam:
    for aln in bam.fetch():
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        if aln.mapping_quality < 20:   # MAPQ 过滤
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

    haps = [w['hap'] for w in windows]
    switch_events = []

    for i in range(1, n_windows):
        if is_switch(windows[i-1], windows[i]):
            switch_events.append({
                'win_idx':   i,
                'win_start': windows[i]['win_start'],
                'from_hap':  windows[i-1]['hap'],
                'to_hap':    windows[i]['hap'],
                'from_contig': windows[i-1]['contig'],
                'to_contig':   windows[i]['contig'],
            })

    n_switches = len(switch_events)
    unique_haps = set(haps)

    if len(unique_haps) == 1:
        status = 'single_hap'
    elif n_switches == 0:
        status = 'cross_contig_correct'   # 跨 contig 但同单倍型
    else:
        status = 'switch_error'

    results.append({
        'read_id':       read_id,
        'n_windows':     n_windows,
        'status':        status,
        'n_switches':    n_switches,
        'hap_sequence':  ','.join(haps),
        'switch_detail': ';'.join(
            f"{e['win_start']}:{e['from_hap']}->{e['to_hap']}" 
            for e in switch_events
        ),
    })

# ---- 输出 TSV ----
with open(out_tsv, 'w') as f:
    f.write('read_id\tn_windows\tstatus\tn_switches\thap_sequence\tswitch_detail\n')
    for r in results:
        f.write(
            f"{r['read_id']}\t{r['n_windows']}\t{r['status']}\t"
            f"{r['n_switches']}\t{r['hap_sequence']}\t{r['switch_detail']}\n"
        )

# ---- 汇总统计 ----
total       = len(results)
single_hap  = sum(1 for r in results if r['status'] == 'single_hap')
cross_ok    = sum(1 for r in results if r['status'] == 'cross_contig_correct')
switch_err  = sum(1 for r in results if r['status'] == 'switch_error')
cross_total = cross_ok + switch_err
total_switches = sum(r['n_switches'] for r in results)

print(f"{'Total reads':<35}: {total:>7,}")
print(f"{'Single haplotype reads':<35}: {single_hap:>7,}")
print(f"{'Cross-contig correct':<35}: {cross_ok:>7,}")
print(f"{'Switch error reads':<35}: {switch_err:>7,}")
print(f"{'Total switch events':<35}: {total_switches:>7,}")
if cross_total > 0:
    print(f"{'Switch error rate (cross-contig)':<35}: {switch_err/cross_total*100:>6.2f}%")

# errors / 100Mb
total_bp = sum(r['n_windows'] for r in results) * 10000
if total_bp > 0:
    print(f"{'Errors / 100 Mb':<35}: {total_switches / total_bp * 1e8:>6.2f}")
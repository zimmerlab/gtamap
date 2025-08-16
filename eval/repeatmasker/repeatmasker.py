#!/usr/bin/env python3

from bitarray import bitarray
import math
from collections import Counter, defaultdict
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

class RepeatMaskerEntry:
    def __init__(self, contig, start, end, rep_class):
        self.contig = contig
        self.start = start
        self.end = end
        self.rep_class = rep_class
    
    def __repr__(self):
        return f"RepeatMaskerEntry(contig='{self.contig}', start={self.start}, end={self.end}, rep_class='{self.rep_class}')"

class Gene:
    def __init__(self, contig, start, end, gene_id, biotype):
        self.contig = contig
        self.start = start
        self.end = end
        self.gene_id = gene_id
        self.biotype = biotype
    
    def __repr__(self):
        return f"Gene(contig='{self.contig}', start={self.start}, end={self.end}, gene_id='{self.gene_id}', biotype='{self.biotype}')"

class Exon:
    def __init__(self, contig, start, end, gene_id):
        self.contig = contig
        self.start = start
        self.end = end
        self.gene_id = gene_id
    
    def __repr__(self):
        return f"Exon(contig='{self.contig}', start={self.start}, end={self.end}, gene_id='{self.gene_id}')"

class CDS:
    def __init__(self, contig, start, end, gene_id):
        self.contig = contig
        self.start = start
        self.end = end
        self.gene_id = gene_id
    
    def __repr__(self):
        return f"CDS(contig='{self.contig}', start={self.start}, end={self.end}, gene_id='{self.gene_id}')"

class GeneInfo:
    def __init__(self, gene_name, biotype):
        self.gene_name = gene_name
        self.biotype = biotype
    
    def __repr__(self):
        return f"GeneInfo(gene_name='{self.gene_name}', biotype='{self.biotype}')"

repeatmasker_file = "/home/sam/Data/repeatmasker/hg38.fa.out"
entries = []

with open(repeatmasker_file, "r") as f:
    for i, line in enumerate(f):
        if line.startswith("#") or i < 3:
            continue
        fields = line.strip().split()
        if len(fields) < 12:
            continue
        
        contig = fields[4]
        if "_" in contig:
            continue
        
        start = int(fields[5])
        end = int(fields[6])
        rep_class = fields[10]
        
        entry = RepeatMaskerEntry(contig, start, end, rep_class)
        entries.append(entry)

print(f"Parsed {len(entries)} entries")

class_counts = Counter(entry.rep_class for entry in entries)
print(f"\nFound {len(class_counts)} unique repeat classes")

sorted_items = class_counts.most_common()
classes = [item[0] for item in sorted_items]
counts = [item[1] for item in sorted_items]

fig = plt.figure(figsize=(12, 8))
bars = plt.bar(classes, counts, fill=False, edgecolor='blue', linewidth=1.2)
plt.xlabel('Repeat Class')
plt.ylabel('Number of Occurrences (log scale)')
plt.title('RepeatMasker: Number of Occurrences per Class')
plt.xticks(rotation=45, ha='right', fontsize=8)
plt.yscale('log')

for bar, count in zip(bars, counts):
    plt.text(bar.get_x() + bar.get_width()/2, bar.get_height(), str(count), 
             ha='center', va='top', rotation=90, fontsize=8)

plt.tight_layout()
plt.savefig("/home/sam/test/py-test/repeatmasker_class_counts.png", dpi=500)
# plt.show()

main_class_counts = Counter()
for entry in entries:
    if '/' in entry.rep_class:
        main_class = entry.rep_class.split('/')[0]
    else:
        main_class = entry.rep_class
    main_class_counts[main_class] += 1

print(f"\nFound {len(main_class_counts)} unique main repeat classes")

sorted_main_items = main_class_counts.most_common()
main_classes = [item[0] for item in sorted_main_items]
main_counts = [item[1] for item in sorted_main_items]

plt.figure(figsize=(12, 8))
main_bars = plt.bar(main_classes, main_counts, fill=False, edgecolor='blue', linewidth=1.5)
plt.xlabel('Main Repeat Class')
plt.ylabel('Number of Occurrences (log scale)')
plt.title('RepeatMasker: Number of Occurrences per Main Class')
plt.xticks(rotation=45, ha='right')
plt.yscale('log')

for bar, count in zip(main_bars, main_counts):
    plt.text(bar.get_x() + bar.get_width()/2, bar.get_height(), str(count), 
             ha='center', va='top', rotation=90, fontsize=16)

plt.tight_layout()
plt.savefig("/home/sam/test/py-test/repeatmasker_main_class_counts.png", dpi=500)
# plt.show()

contig_main_class_counts = defaultdict(lambda: defaultdict(int))
for entry in entries:
    if '/' in entry.rep_class:
        main_class = entry.rep_class.split('/')[0]
    else:
        main_class = entry.rep_class
    contig_main_class_counts[entry.contig][main_class] += 1

def sort_chromosomes(contig_list):
    def chr_key(contig):
        if contig.startswith('chr'):
            chr_part = contig[3:]
            if chr_part.isdigit():
                return (0, int(chr_part))
            elif chr_part == 'X':
                return (1, 0)
            elif chr_part == 'Y':
                return (1, 1)
            else:
                return (2, chr_part)
        return (3, contig)
    
    return sorted(contig_list, key=chr_key)

contigs = sort_chromosomes(list(contig_main_class_counts.keys()))
# top_main_classes = [item[0] for item in main_class_counts.most_common(10)]
top_main_classes = [item[0] for item in main_class_counts.most_common()]

data_matrix = []
for contig in contigs:
    row = [contig_main_class_counts[contig][main_class] for main_class in top_main_classes]
    data_matrix.append(row)

data_matrix = np.array(data_matrix)

contig_totals = data_matrix.sum(axis=1)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 10), gridspec_kw={'width_ratios': [1, 10]})
ax1.barh(range(len(contigs))[::-1], contig_totals)
ax1.set_xlabel('Total Occurrences')
ax1.set_ylim(-0.5, len(contigs) - 0.5)
ax1.set_yticks([])
ax1.invert_xaxis()

im = ax2.imshow(data_matrix, cmap='plasma', aspect='auto')
plt.colorbar(im, ax=ax2, label='Number of Occurrences')
ax2.set_xlabel('Main Repeat Class')
ax2.set_ylabel('Contig')
ax2.set_title('Distribution of Main Classes Across Contigs')
ax2.set_xticks(range(len(top_main_classes)))
ax2.set_xticklabels(top_main_classes, rotation=45, ha='right')
ax2.set_yticks(range(len(contigs)))
ax2.set_yticklabels(contigs)

for i in range(len(contigs)):
    for j in range(len(top_main_classes)):
        text = ax2.text(j, i, str(int(data_matrix[i, j])), 
                       ha='center', va='center', color='white', fontsize=8)

plt.tight_layout()
plt.savefig("/home/sam/test/py-test/repeatmasker_contig_distribution.png", dpi=500)
# plt.show()




size_data = defaultdict(list)
for entry in entries:
    if '/' in entry.rep_class:
        main_class = entry.rep_class.split('/')[0]
    else:
        main_class = entry.rep_class
    size = entry.end - entry.start
    size_data[main_class].append(size)

top_classes_for_sizes = [item[0] for item in main_class_counts.most_common()]
# exlude "Satellite"
top_classes_for_sizes = [cls for cls in top_classes_for_sizes if cls != "Satellite"]

plt.figure(figsize=(14, 10))
box_data = [size_data[main_class] for main_class in top_classes_for_sizes]
bp = plt.boxplot(box_data, labels=top_classes_for_sizes, patch_artist=True, showfliers=False)

for patch in bp['boxes']:
    patch.set_facecolor('lightblue')
    patch.set_alpha(0.7)

plt.xlabel('Main Repeat Class')
plt.ylabel('Repeat Size (bp)')
plt.title('Distribution of Repeat Sizes by Main Class')
plt.xticks(rotation=45, ha='right')
# plt.yscale('log')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/home/sam/test/py-test/repeatmasker_size_distribution.png", dpi=500)
# plt.show()





size_data_all = defaultdict(list)
for entry in entries:
    size = entry.end - entry.start
    size_data_all[entry.rep_class].append(size)

all_classes = [item[0] for item in class_counts.most_common()]
all_classes = [cls for cls in all_classes if "Satellite" not in cls]

plt.figure(figsize=(20, 12))
box_data_all = [size_data_all[rep_class] for rep_class in all_classes]
bp_all = plt.boxplot(box_data_all, labels=all_classes, patch_artist=True, showfliers=False)

for patch in bp_all['boxes']:
    patch.set_facecolor('lightcoral')
    patch.set_alpha(0.7)

plt.xlabel('Repeat Class')
plt.ylabel('Repeat Size (bp)')
plt.title('Distribution of Repeat Sizes by All Classes (no Satellite)')
plt.xticks(rotation=90, ha='right', fontsize=6)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/home/sam/test/py-test/repeatmasker_size_distribution_all.png", dpi=500)
# plt.show()



genes_file = "/home/sam/test/py-test/genes.csv"
genes = []
exons = []
cds_regions = []

with open(genes_file, "r") as f:
    for line in f:
        fields = line.strip().split(',')
        contig = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        feature_type = fields[3]
        gene_id = fields[4]
        
        biotype = fields[5] if len(fields) > 5 else "unknown"
        
        if feature_type == "gene":
            gene = Gene(contig, start, end, gene_id, biotype)
            genes.append(gene)
        elif feature_type == "exon":
            exon = Exon(contig, start, end, gene_id)
            exons.append(exon)
        elif feature_type == "CDS":
            cds = CDS(contig, start, end, gene_id)
            cds_regions.append(cds)

print(f"Loaded {len(genes)} genes, {len(exons)} exons, and {len(cds_regions)} CDS regions")


gene_names_file = "/home/sam/test/py-test/gene-names.csv"
gene_id_to_info = {}

with open(gene_names_file, "r") as f:
    for line in f:
        fields = line.strip().split(',')
        if len(fields) >= 3:
            gene_id = fields[0]
            gene_name = fields[1]
            biotype = fields[2]
            gene_id_to_info[gene_id] = GeneInfo(gene_name, biotype)
        elif len(fields) >= 2:
            gene_id = fields[0]
            gene_name = fields[1]
            gene_id_to_info[gene_id] = GeneInfo(gene_name, "unknown")

print(f"Loaded {len(gene_id_to_info)} gene ID to info mappings")

gene_intervals = defaultdict(list)
for gene in genes:
    gene_intervals[gene.contig].append((gene.start, gene.end, gene))

for contig in gene_intervals:
    gene_intervals[contig].sort(key=lambda x: x[0])


def binary_search(genes, start, end) -> int:
    low, high = 0, len(genes) - 1
    while low <= high:
        mid = (low + high) // 2
        if genes[mid][0] < end and genes[mid][1] > start:
            return mid
        elif genes[mid][0] >= end:
            high = mid - 1
        else:
            low = mid + 1
    return -1

def find_overlapping_genes(repeat):
    overlapping_genes = []
    if repeat.contig not in gene_intervals:
        return overlapping_genes
   
    index_first = binary_search(gene_intervals[repeat.contig], repeat.start, repeat.end)

    if index_first == -1:
        return overlapping_genes

    for i in range(index_first, len(gene_intervals[repeat.contig])):
        gene_start, gene_end, gene = gene_intervals[repeat.contig][i]
        if gene_start >= repeat.end:
            break
        if repeat.start < gene_end and repeat.end > gene_start:
            overlapping_genes.append(gene)
    #
    # for gene_start, gene_end, gene in gene_intervals[repeat.contig]:
    #     if gene_start >= repeat.end:
    #         break
    #     if repeat.start < gene_end and repeat.end > gene_start:
    #         overlapping_genes.append(gene)
    
    return overlapping_genes

gene_repeat_counts = defaultdict(int)
for i, repeat in enumerate(entries):
    if i % 100 == 0:
        print(f"Processing repeat {i+1}/{len(entries)}")
    overlapping_genes = find_overlapping_genes(repeat)
    for gene in overlapping_genes:
        gene_repeat_counts[gene.gene_id] += 1

repeat_counts = [gene_repeat_counts[gene.gene_id] for gene in genes]

plt.figure(figsize=(12, 8))
plt.hist(repeat_counts, bins=50, edgecolor='black', alpha=0.7)
plt.xlabel('Number of Repeats per Gene')
plt.ylabel('Number of Genes')
plt.title('Distribution of Repeats per Gene')
plt.yscale('log')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/home/sam/test/py-test/repeats_per_gene_distribution.png", dpi=500)
# plt.show()

print(f"Mean repeats per gene: {np.mean(repeat_counts):.2f}")
print(f"Median repeats per gene: {np.median(repeat_counts):.2f}")
print(f"Max repeats in a gene: {max(repeat_counts)}")
print(f"Genes with no repeats: {repeat_counts.count(0)}")




exon_intervals = defaultdict(list)
for exon in exons:
    exon_intervals[exon.contig].append((exon.start, exon.end, exon))

for contig in exon_intervals:
    exon_intervals[contig].sort(key=lambda x: x[0])

def find_overlapping_exons(repeat):
    overlapping_exons = []
    if repeat.contig not in exon_intervals:
        return overlapping_exons
   
    index_first = binary_search(exon_intervals[repeat.contig], repeat.start, repeat.end)

    if index_first == -1:
        return overlapping_exons

    for i in range(index_first, len(exon_intervals[repeat.contig])):
        exon_start, exon_end, exon = exon_intervals[repeat.contig][i]
        if exon_start >= repeat.end:
            break
        if repeat.start < exon_end and repeat.end > exon_start:
            overlapping_exons.append(exon)

    # for exon_start, exon_end, exon in exon_intervals[repeat.contig]:
    #     if exon_start >= repeat.end:
    #         break
    #     if repeat.start < exon_end and repeat.end > exon_start:
    #         overlapping_exons.append(exon)
    
    return overlapping_exons

cds_intervals = defaultdict(list)
for cds in cds_regions:
    cds_intervals[cds.contig].append((cds.start, cds.end, cds))

for contig in cds_intervals:
    cds_intervals[contig].sort(key=lambda x: x[0])

def find_overlapping_cds(repeat):
    overlapping_cds = []
    if repeat.contig not in cds_intervals:
        return overlapping_cds
    index_first = binary_search(cds_intervals[repeat.contig], repeat.start, repeat.end)
    if index_first == -1:
        return overlapping_cds
    for i in range(index_first, len(cds_intervals[repeat.contig])):
        cds_start, cds_end, cds = cds_intervals[repeat.contig][i]
        if cds_start >= repeat.end:
            break
        if repeat.start < cds_end and repeat.end > cds_start:
            overlapping_cds.append(cds)
    return overlapping_cds

repeats_overlap_gene = 0
repeats_overlap_gene_and_exon = 0 
repeats_overlap_gene_not_exon = 0
repeats_no_overlap = 0

categorized_repeats = []

def categorize_repeat(repeat, genes, exons, cds):
   
    overlapping_cds = set()
    overlapping_utr = set()
    in_intron = set()

    if len(genes) == 0:
        return (repeat, overlapping_cds, overlapping_utr, in_intron)

    for c in cds:
        overlapping_cds.add(c.gene_id)
    for e in exons:
        if e.gene_id in overlapping_cds:
            continue
        overlapping_utr.add(e.gene_id)
    for gene in genes:
        if gene.gene_id in overlapping_cds or gene.gene_id in overlapping_utr:
            continue
        in_intron.add(gene.gene_id)
    
    return (repeat, overlapping_cds, overlapping_utr, in_intron)

    
for i, repeat in enumerate(entries):
    if i % 1000 == 0:
        print(f"Categorizing repeat {i+1}/{len(entries)}")
    
    overlapping_genes = find_overlapping_genes(repeat)
    overlapping_exons = find_overlapping_exons(repeat)
    overlapping_cds = find_overlapping_cds(repeat)
    
    if len(overlapping_genes) > 0:
        repeats_overlap_gene += 1
        if len(overlapping_exons) > 0:
            repeats_overlap_gene_and_exon += 1
        else:
            repeats_overlap_gene_not_exon += 1
    else:
        repeats_no_overlap += 1

    categorized_repeat = categorize_repeat(repeat, overlapping_genes, overlapping_exons, overlapping_cds)
    categorized_repeats.append(categorized_repeat)

categories = ['Overlap Gene\n& Exon', 'Overlap Gene\nNot Exon', 'No Gene\nOverlap']
counts = [repeats_overlap_gene_and_exon, repeats_overlap_gene_not_exon, repeats_no_overlap]
colors = ['darkgreen', 'lightgreen', 'lightcoral']

plt.figure(figsize=(10, 8))
bars = plt.bar(categories, counts, color=colors, alpha=0.8, edgecolor='black')
plt.xlabel('Repeat Overlap Category')
plt.ylabel('Number of Repeats')
plt.title('Distribution of Repeat Elements by Gene/Exon Overlap')

for bar, count in zip(bars, counts):
    plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(counts)*0.01, 
             f'{count:,}', ha='center', va='bottom', fontweight='bold')

plt.tight_layout()
plt.savefig("/home/sam/test/py-test/repeat_gene_exon_overlap.png", dpi=500)
# plt.show()

print(f"\nRepeat overlap summary:")
print(f"Total repeats: {len(entries):,}")
print(f"Overlap gene & exon: {repeats_overlap_gene_and_exon:,} ({repeats_overlap_gene_and_exon/len(entries)*100:.1f}%)")
print(f"Overlap gene not exon: {repeats_overlap_gene_not_exon:,} ({repeats_overlap_gene_not_exon/len(entries)*100:.1f}%)")
print(f"No gene overlap: {repeats_no_overlap:,} ({repeats_no_overlap/len(entries)*100:.1f}%)")



for repeat, cds, exons, introns in categorized_repeats[:3]:
    print(repeat)
    print(cds)
    print(exons)
    print(introns)
    print("-")

repeat_class_categories = defaultdict(lambda: {'cds': 0, 'exon': 0, 'intron': 0, 'no_gene': 0})

for repeat, overlapping_cds, overlapping_utr, in_intron in categorized_repeats:
    rep_class = repeat.rep_class
    
    if len(overlapping_cds) > 0:
        repeat_class_categories[rep_class]['cds'] += 1
    elif len(overlapping_utr) > 0:
        repeat_class_categories[rep_class]['exon'] += 1
    elif len(in_intron) > 0:
        repeat_class_categories[rep_class]['intron'] += 1
    else:
        repeat_class_categories[rep_class]['no_gene'] += 1

top_classes_for_categorization = [item[0] for item in class_counts.most_common(15)]

cds_counts = [repeat_class_categories[rep_class]['cds'] for rep_class in top_classes_for_categorization]
exon_counts = [repeat_class_categories[rep_class]['exon'] for rep_class in top_classes_for_categorization]
intron_counts = [repeat_class_categories[rep_class]['intron'] for rep_class in top_classes_for_categorization]
no_gene_counts = [repeat_class_categories[rep_class]['no_gene'] for rep_class in top_classes_for_categorization]

x = np.arange(len(top_classes_for_categorization)) * 1.2  # Add margin between classes
width = 0.2

fig = plt.figure(figsize=(18, 10))
bars1 = plt.bar(x - width*1.5, cds_counts, width, label='CDS', facecolor='red', edgecolor='darkred', linewidth=1.2, alpha=0.8)
bars2 = plt.bar(x - width*0.5, exon_counts, width, label='UTR', facecolor='green', edgecolor='darkgreen', linewidth=1.2, alpha=0.8)
bars3 = plt.bar(x + width*0.5, intron_counts, width, label='Intron', facecolor='blue', edgecolor='darkblue', linewidth=1.2, alpha=0.8)
bars4 = plt.bar(x + width*1.5, no_gene_counts, width, label='No Gene', facecolor='lightgray', edgecolor='gray', linewidth=1.2, alpha=0.8)

plt.xlabel('Repeat Class')
plt.ylabel('Number of Repeats')
plt.title('Distribution of Repeat Elements by Genomic Region (top 15)')
plt.xticks(x, top_classes_for_categorization, rotation=45, ha='right', fontsize=8)
# plt.yscale('log')
plt.legend()
plt.grid(True, alpha=0.3)

for bars in [bars1, bars2, bars3, bars4]:
    for bar in bars:
        height = bar.get_height()
        if height > 0:
            plt.text(bar.get_x() + bar.get_width()/2, bar.get_height()+5000, str(int(height)), 
                    ha='center', va='bottom', rotation=90, fontsize=8)

plt.tight_layout()
plt.savefig("/home/sam/test/py-test/repeat_genomic_categories.png", dpi=500)




main_class_categories = defaultdict(lambda: {'cds': 0, 'exon': 0, 'intron': 0, 'no_gene': 0})

for repeat, overlapping_cds, overlapping_utr, in_intron in categorized_repeats:
    if '/' in repeat.rep_class:
        main_class = repeat.rep_class.split('/')[0]
    else:
        main_class = repeat.rep_class
    
    if len(overlapping_cds) > 0:
        main_class_categories[main_class]['cds'] += 1
    elif len(overlapping_utr) > 0:
        main_class_categories[main_class]['exon'] += 1
    elif len(in_intron) > 0:
        main_class_categories[main_class]['intron'] += 1
    else:
        main_class_categories[main_class]['no_gene'] += 1

top_main_classes_for_categorization = [item[0] for item in main_class_counts.most_common()]

main_cds_counts = [main_class_categories[main_class]['cds'] for main_class in top_main_classes_for_categorization]
main_exon_counts = [main_class_categories[main_class]['exon'] for main_class in top_main_classes_for_categorization]
main_intron_counts = [main_class_categories[main_class]['intron'] for main_class in top_main_classes_for_categorization]
main_no_gene_counts = [main_class_categories[main_class]['no_gene'] for main_class in top_main_classes_for_categorization]

x_main = np.arange(len(top_main_classes_for_categorization)) * 1.2
width_main = 0.2

fig = plt.figure(figsize=(16, 10))
bars1_main = plt.bar(x_main - width_main*1.5, main_cds_counts, width_main, label='CDS', facecolor='red', edgecolor='darkred', linewidth=1.2, alpha=0.8)
bars2_main = plt.bar(x_main - width_main*0.5, main_exon_counts, width_main, label='UTR', facecolor='green', edgecolor='darkgreen', linewidth=1.2, alpha=0.8)
bars3_main = plt.bar(x_main + width_main*0.5, main_intron_counts, width_main, label='Intron', facecolor='blue', edgecolor='darkblue', linewidth=1.2, alpha=0.8)
bars4_main = plt.bar(x_main + width_main*1.5, main_no_gene_counts, width_main, label='No Gene', facecolor='lightgray', edgecolor='gray', linewidth=1.2, alpha=0.8)

plt.xlabel('Main Repeat Class')
plt.ylabel('Number of Repeats')
plt.title('Distribution of Main Repeat Classes by Genomic Region')
plt.xticks(x_main, top_main_classes_for_categorization, rotation=45, ha='right', fontsize=10)
plt.legend()
plt.grid(True, alpha=0.3)

for bars, counts in [(bars1_main, main_cds_counts), (bars2_main, main_exon_counts), (bars3_main, main_intron_counts), (bars4_main, main_no_gene_counts)]:
    for bar, count in zip(bars, counts):
        if count > 0:
            plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(max(main_cds_counts), max(main_exon_counts), max(main_intron_counts), max(main_no_gene_counts))*0.01, str(int(count)), 
                    ha='center', va='bottom', rotation=90, fontsize=8)

plt.tight_layout()
plt.savefig("/home/sam/test/py-test/main_class_genomic_categories_absolute.png", dpi=500)
# plt.show()



n_top_genes = 30

gene_repeat_counts = defaultdict(lambda: {'cds': 0, 'exon': 0, 'intron': 0, 'total': 0})

for repeat, overlapping_cds, overlapping_utr, in_intron in categorized_repeats:
    all_gene_ids = set()
    
    for gene_id in overlapping_cds:
        gene_repeat_counts[gene_id]['cds'] += 1
        all_gene_ids.add(gene_id)
    
    for gene_id in overlapping_utr:
        gene_repeat_counts[gene_id]['exon'] += 1
        all_gene_ids.add(gene_id)
    
    for gene_id in in_intron:
        gene_repeat_counts[gene_id]['intron'] += 1
        all_gene_ids.add(gene_id)
    
    for gene_id in all_gene_ids:
        gene_repeat_counts[gene_id]['total'] += 1

top_cds_genes = sorted(gene_repeat_counts.items(), key=lambda x: x[1]['cds'], reverse=True)[:n_top_genes]
top_exon_genes = sorted(gene_repeat_counts.items(), key=lambda x: x[1]['exon'], reverse=True)[:n_top_genes]
top_intron_genes = sorted(gene_repeat_counts.items(), key=lambda x: x[1]['intron'], reverse=True)[:n_top_genes]

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(20, 20))

cds_gene_names = []
for gene_id, _ in top_cds_genes:
    gene_info = gene_id_to_info.get(gene_id)
    if gene_info:
        gene_label = f"{gene_info.gene_name} ({gene_info.biotype})"
    else:
        gene_label = gene_id
    if len(gene_label) > 25:
        gene_label = gene_label[:25] + '...'
    cds_gene_names.append(gene_label)
cds_counts = [counts['cds'] for _, counts in top_cds_genes]
bars1 = ax1.bar(cds_gene_names, cds_counts, facecolor='red', edgecolor='darkred', linewidth=1.2, alpha=0.8)
ax1.set_xlabel('Gene Name (Biotype)')
ax1.set_ylabel('Number of Repeats in CDS')
ax1.set_title(f'Top {n_top_genes} Genes with Most Repeats in CDS Regions')
ax1.tick_params(axis='x', rotation=90, labelsize=6)
for bar, count in zip(bars1, cds_counts):
    if count > 0:
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(cds_counts)*0.01, str(count), 
                ha='center', va='bottom', fontsize=6)

exon_gene_names = []
for gene_id, _ in top_exon_genes:
    gene_info = gene_id_to_info.get(gene_id)
    if gene_info:
        gene_label = f"{gene_info.gene_name} ({gene_info.biotype})"
    else:
        gene_label = gene_id
    if len(gene_label) > 25:
        gene_label = gene_label[:25] + '...'
    exon_gene_names.append(gene_label)
exon_counts = [counts['exon'] for _, counts in top_exon_genes]
bars2 = ax2.bar(exon_gene_names, exon_counts, facecolor='green', edgecolor='darkgreen', linewidth=1.2, alpha=0.8)
ax2.set_xlabel('Gene Name (Biotype)')
ax2.set_ylabel('Number of Repeats in Exons (UTR)')
ax2.set_title(f'Top {n_top_genes} Genes with Most Repeats in Exon Regions')
ax2.tick_params(axis='x', rotation=90, labelsize=6)
for bar, count in zip(bars2, exon_counts):
    if count > 0:
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(exon_counts)*0.01, str(count), 
                ha='center', va='bottom', fontsize=6)

intron_gene_names = []
for gene_id, _ in top_intron_genes:
    gene_info = gene_id_to_info.get(gene_id)
    if gene_info:
        gene_label = f"{gene_info.gene_name} ({gene_info.biotype})"
    else:
        gene_label = gene_id
    if len(gene_label) > 25:
        gene_label = gene_label[:25] + '...'
    intron_gene_names.append(gene_label)
intron_counts = [counts['intron'] for _, counts in top_intron_genes]
bars3 = ax3.bar(intron_gene_names, intron_counts, facecolor='blue', edgecolor='darkblue', linewidth=1.2, alpha=0.8)
ax3.set_xlabel('Gene Name (Biotype)')
ax3.set_ylabel('Number of Repeats in Introns')
ax3.set_title(f'Top {n_top_genes} Genes with Most Repeats in Intronic Regions')
ax3.tick_params(axis='x', rotation=90, labelsize=6)
for bar, count in zip(bars3, intron_counts):
    if count > 0:
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(intron_counts)*0.01, str(count), 
                ha='center', va='bottom', fontsize=6)

plt.tight_layout()
plt.savefig("/home/sam/test/py-test/top_genes_with_repeats.png", dpi=500)
# plt.show()




genes_no_repeats = set()
genes_with_cds_repeats = set()
genes_with_exon_repeats = set()
genes_with_intron_repeats = set()

for repeat, cds, exon, intron in categorized_repeats:
    if len(cds) == 0:
        continue
    for gene_id in cds:
        genes_with_cds_repeats.add(gene_id)

for repeat, cds, exon, intron in categorized_repeats:
    if len(exon) == 0:
        continue
    for gene_id in exon:
        if gene_id in genes_with_cds_repeats:
            continue
        genes_with_exon_repeats.add(gene_id)

for repeat, cds, exon, intron in categorized_repeats:
    if len(intron) == 0:
        continue
    for gene_id in intron:
        if gene_id in genes_with_cds_repeats or gene_id in genes_with_exon_repeats:
            continue
        genes_with_intron_repeats.add(gene_id)

for g in genes:
    if g.gene_id not in genes_with_cds_repeats and g.gene_id not in genes_with_exon_repeats and g.gene_id not in genes_with_intron_repeats:
        genes_no_repeats.add(g.gene_id)


print(len(genes_with_cds_repeats))
print(len(genes_with_exon_repeats))
print(len(genes_with_intron_repeats))
print(len(genes_no_repeats))
print(len(genes))

categories = ['No Repeats', 'CDS Repeats', 'Exon Repeats', 'Intron Repeats']
gene_counts = [len(genes_no_repeats), len(genes_with_cds_repeats), len(genes_with_exon_repeats), len(genes_with_intron_repeats)]
colors = ['lightgray', 'red', 'green', 'blue']

fig = plt.figure(figsize=(12, 8))
bars = plt.bar(categories, gene_counts, facecolor=colors, edgecolor='black', linewidth=1.2, alpha=0.8)
plt.xlabel('Repeat Category')
plt.ylabel('Number of Genes')
plt.title('Number of Genes by Repeat Content')
plt.grid(True, alpha=0.3)

for bar, count in zip(bars, gene_counts):
    plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(gene_counts)*0.01, 
             f'{count:,}', ha='center', va='bottom', fontweight='bold', fontsize=10)

plt.tight_layout()
plt.savefig("/home/sam/test/py-test/genes_by_repeat_content.png", dpi=500)

biotype_counts_by_category = {
    'No Repeats': defaultdict(int),
    'CDS Repeats': defaultdict(int), 
    'Exon Repeats': defaultdict(int),
    'Intron Repeats': defaultdict(int)
}

for gene_id in genes_no_repeats:
    gene_info = gene_id_to_info.get(gene_id)
    if gene_info:
        biotype_counts_by_category['No Repeats'][gene_info.biotype] += 1

for gene_id in genes_with_cds_repeats:
    gene_info = gene_id_to_info.get(gene_id)
    if gene_info:
        biotype_counts_by_category['CDS Repeats'][gene_info.biotype] += 1

for gene_id in genes_with_exon_repeats:
    gene_info = gene_id_to_info.get(gene_id)
    if gene_info:
        biotype_counts_by_category['Exon Repeats'][gene_info.biotype] += 1

for gene_id in genes_with_intron_repeats:
    gene_info = gene_id_to_info.get(gene_id)
    if gene_info:
        biotype_counts_by_category['Intron Repeats'][gene_info.biotype] += 1

all_biotypes = set()
for category_biotypes in biotype_counts_by_category.values():
    all_biotypes.update(category_biotypes.keys())

top_biotypes = sorted(all_biotypes, key=lambda bt: sum(biotype_counts_by_category[cat][bt] for cat in categories), reverse=True)[:8]

biotype_colors = plt.cm.Set3(np.linspace(0, 1, len(top_biotypes)))

fig, ax = plt.subplots(figsize=(16, 10))
width = 0.8
x = np.arange(len(categories))

bottoms = np.zeros(len(categories))
for i, biotype in enumerate(top_biotypes):
    heights = [biotype_counts_by_category[cat][biotype] for cat in categories]
    bars = ax.bar(x, heights, width, bottom=bottoms, label=biotype, 
                  color=biotype_colors[i], alpha=0.8, edgecolor='black', linewidth=0.5)
    
    for j, (bar, height) in enumerate(zip(bars, heights)):
        if height > 100:  # Only show labels for substantial counts
            ax.text(bar.get_x() + bar.get_width()/2, bottoms[j] + height/2, str(height),
                   ha='center', va='center', fontsize=8, fontweight='bold')
    
    bottoms += heights

ax.set_xlabel('Repeat Category')
ax.set_ylabel('Number of Genes')
ax.set_title('Number of Genes by Repeat Content and Biotype')
ax.set_xticks(x)
ax.set_xticklabels(categories)
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("/home/sam/test/py-test/genes_by_repeat_content_biotype.png", dpi=500, bbox_inches='tight')
# plt.show()



def calculate_overlap_length(range1_start, range1_end, range2_start, range2_end):
    """Calculate the length of overlap between two ranges"""
    overlap_start = max(range1_start, range2_start)
    overlap_end = min(range1_end, range2_end)
    return max(0, overlap_end - overlap_start)

def merge_overlapping_ranges(ranges):
    """Merge overlapping ranges to avoid double-counting"""
    if not ranges:
        return []
    
    sorted_ranges = sorted(ranges)
    merged = [sorted_ranges[0]]
    
    for start, end in sorted_ranges[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    
    return merged

gene_vecs = dict()

for gene in genes:
    # print(gene.gene_id)
    # print(gene.start, gene.end, gene.end-gene.start)

    gene_len = gene.end - gene.start + 1

    g_vec = bitarray(gene.end - gene.start + 1)
    g_vec.setall(False)

    gene_vecs[gene.gene_id] = (gene, g_vec)

c = 0
for repeat, cds, exon, intron in categorized_repeats:
    c += 1

    print(f"Processing repeat {c}/{len(categorized_repeats)}")

    ids = cds.union(exon).union(intron)

    if len(ids) == 0:
        continue

    # if gene.gene_id not in cds and gene.gene_id not in exon and gene.gene_id not in intron:
    #     continue

    # print(repeat)
    # print(cds)

    repeat_vec = 

    for gene_id in ids:

        gene, g_vec = gene_vecs[gene_id]

        repeat_len = repeat.end - repeat.start + 1
        repeat_start_relative = repeat.start - gene.start

        # print(repeat_len)
        # print(repeat_start_relative)

        if repeat_start_relative < 0:
            repeat_len += repeat_start_relative
            repeat_start_relative = 0

        overlap_right = (repeat_start_relative + repeat_len) - gene.end
        if overlap_right > 0:
            repeat_len -= overlap_right

        # print("updated")
        # print(repeat_start_relative,":",repeat_len)

        g_vec[repeat_start_relative:repeat_start_relative+repeat_len] = True
        # break

covered = 0
uncovered = 0
for b in g_vec:
    if b:
        covered += 1
    else:
        uncovered += 1

print(f"Covered: {covered}, Uncovered: {uncovered}, Total: {len(g_vec)}")

# break










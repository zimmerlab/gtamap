# GTAMap makefile

MAIN := ./main
FASTA := ./fasta/ENSG00000141510.19.fa
MASK := ./data/mask.full.bed
INDEX_DIR := ./index
INDEX_FILE := $(INDEX_DIR)/index.gtai
FW_FASTQ := ./data/fw.fastq
RW_FASTQ := ./data/rv.fastq
OUTPUT_DIR := output
SAM_FILE := $(OUTPUT_DIR)/aligned.sam
BAM_FILE := TP53.bam
GENE_ID := ENSG00000141510
EVAL_SCRIPT := ../eval/compare_to_ground_truth.py
GROUND_TRUTH := data/read.mappinginfo
RESULT_DIR := res

# phony targets
.PHONY: all index map eval prep_igv clean help

# default target
all: prep_igv eval

# create index
index: $(INDEX_FILE)

$(INDEX_FILE): $(FASTA) $(MASK)
	@echo "Creating index..."
	@mkdir -p $(INDEX_DIR)
	$(MAIN) index -f $(FASTA) -o $(INDEX_DIR) -m $(MASK)

# run mapping
map: $(SAM_FILE)

$(SAM_FILE): $(INDEX_FILE) $(FW_FASTQ) $(RW_FASTQ)
	@echo "Running mapping..."
	@mkdir -p $(OUTPUT_DIR)
	$(MAIN) map -i $(INDEX_FILE) -1 $(FW_FASTQ) -2 $(RW_FASTQ) -o $(OUTPUT_DIR) -r rna #--loglevel DEBUG

# run evaluation
eval: $(SAM_FILE)
	@echo "Running evaluation..."
	@mkdir -p $(RESULT_DIR)
	python3 $(EVAL_SCRIPT) --sam1 $(SAM_FILE) --ground_truth $(GROUND_TRUTH) --output_dir $(RESULT_DIR) --mapper "gtamap" --gene_id $(GENE_ID) #--verbose 

# prepare for IGV
prep_igv: $(BAM_FILE)


$(BAM_FILE): $(SAM_FILE)
	@echo "Preparing BAM file for IGV..."
	samtools view -b $(SAM_FILE) | samtools sort -o $(BAM_FILE)
	samtools index $(BAM_FILE)

# clean all generated files
clean:
	@echo "Cleaning up..."
	rm -rf $(INDEX_DIR) $(OUTPUT_DIR) $(RESULT_DIR) $(BAM_FILE) $(BAM_FILE).bai

# clean only output files (keep index)
clean-output:
	@echo "Cleaning output files..."
	rm -rf $(OUTPUT_DIR) $(RESULT_DIR) $(BAM_FILE) $(BAM_FILE).bai

# help target
help:
	@echo "GTAMap Pipeline Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  all        - Run complete pipeline (index, map, eval, prep_igv)"
	@echo "  index      - Create index from FASTA file"
	@echo "  map        - Run mapping on paired-end reads"
	@echo "  eval       - Evaluate mapping results against ground truth"
	@echo "  prep_igv   - Prepare BAM file for IGV visualization"
	@echo "  clean      - Remove all generated files"
	@echo "  clean-output - Remove output files but keep index"
	@echo "  help       - Show this help message"


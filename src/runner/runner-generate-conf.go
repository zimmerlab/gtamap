package runner

import (
	"fmt"
	"os"
	"path/filepath"

	"github.com/sirupsen/logrus"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
	"gopkg.in/yaml.v3"
)

func GetCommandGenerateConf() *cobra.Command {
	var outputDir string
	var fileName string

	generateConfCmd := &cobra.Command{
		Use:   "generate-conf",
		Short: "Generate a YAML config file with default parameters",
		Run: func(cmd *cobra.Command, args []string) {
			ExecGenerateConf(outputDir, fileName)
		},
	}

	flags := generateConfCmd.Flags()

	flags.StringVarP(
		&outputDir,
		"output",
		"o",
		"",
		"Output directory for the config file (required)",
	)
	generateConfCmd.MarkFlagRequired("output")

	flags.StringVarP(
		&fileName,
		"filename",
		"f",
		"gtamap_config.yaml",
		"Config file name (default: gtamap_config.yaml)",
	)

	return generateConfCmd
}

func ExecGenerateConf(outputDir string, fileName string) {
	// Create output directory if it doesn't exist
	if err := os.MkdirAll(outputDir, 0755); err != nil {
		logrus.Fatalf("Failed to create output directory: %v", err)
	}

	// Build the full output path
	outputPath := filepath.Join(outputDir, fileName)

	// Create a map with all default values
	configMap := map[string]interface{}{
		"general": map[string]interface{}{
			"output_dir": viper.GetString("general.output_dir"),
		},
		"index": map[string]interface{}{
			"fasta_file_path":       "",
			"fasta_index_file_path": "",
			"gtf_file_path":         "",
			"gene_ids":              viper.Get("index.gene_ids"),
			"regions":               viper.Get("index.regions"),
			"upstream_bases":        viper.GetInt("index.upstream_bases"),
			"downstream_bases":      viper.GetInt("index.downstream_bases"),
			"regionmask_file_path":  "",
			"output": map[string]interface{}{
				"single_file":         viper.GetBool("index.output.single_file"),
				"fasta_file_name":     viper.GetString("index.output.fasta_file_name"),
				"index_file_name":     viper.GetString("index.output.index_file_name"),
				"use_fasta_file_name": false,
			},
		},
		"mapping": map[string]interface{}{
			"index_file_path":       "",
			"fastq_r1_file_path":    "",
			"fastq_r2_file_path":    "",
			"read_origin":           "",
			"threads":               viper.GetInt("mapping.threads"),
			"regionmask_file_path":  "",
			"max_branching_depth":   viper.GetInt("mapping.max_branching_depth"),
			"rna_mode": map[string]interface{}{
				"filter_min_matches":      viper.GetInt("mapping.rna_mode.filter_min_matches"),
				"intron_length_min":       viper.GetInt("mapping.rna_mode.intron_length_min"),
				"max_mismatch_count":      viper.GetInt("mapping.rna_mode.max_mismatch_count"),
				"max_mismatch_percentage": viper.GetFloat64("mapping.rna_mode.max_mismatch_percentage"),
				"confident": map[string]interface{}{
					"max_mismatch_count":            viper.GetInt("mapping.rna_mode.confident.max_mismatch_count"),
					"min_anchor_length":             viper.GetInt("mapping.rna_mode.confident.min_anchor_length"),
					"intron_cluster_delta":          viper.GetInt("mapping.rna_mode.confident.intron_cluster_delta"),
					"intron_cluster_repair_window":  viper.GetInt("mapping.rna_mode.confident.intron_cluster_repair_window"),
				},
			},
			"dna_mode": map[string]interface{}{
				"filter_min_matches":         viper.GetInt("mapping.dna_mode.filter_min_matches"),
				"min_length_initial_diagonal": viper.GetFloat64("mapping.dna_mode.min_length_initial_diagonal"),
				"max_gap_length":             viper.GetInt("mapping.dna_mode.max_gap_length"),
				"max_gap_count":              viper.GetInt("mapping.dna_mode.max_gap_count"),
				"max_mismatch_count":         viper.GetInt("mapping.dna_mode.max_mismatch_count"),
				"max_mismatch_percentage":    viper.GetFloat64("mapping.dna_mode.max_mismatch_percentage"),
				"confident": map[string]interface{}{
					"max_mismatch_count":           viper.GetInt("mapping.dna_mode.confident.max_mismatch_count"),
					"min_anchor_length":            viper.GetInt("mapping.dna_mode.confident.min_anchor_length"),
					"intron_cluster_delta":         viper.GetInt("mapping.dna_mode.confident.intron_cluster_delta"),
					"intron_cluster_repair_window": viper.GetInt("mapping.dna_mode.confident.intron_cluster_repair_window"),
				},
			},
			"output": map[string]interface{}{
				"mapping_progress_file_name": viper.GetString("mapping.output.mapping_progress_file_name"),
				"mapping_stats_file_name":    viper.GetString("mapping.output.mapping_stats_file_name"),
				"sam_file_name":              viper.GetString("mapping.output.sam_file_name"),
				"sam_use_x":                  viper.GetBool("mapping.output.sam_use_x"),
				"sam_include_all_pairings":   viper.GetBool("mapping.output.sam_include_all_pairings"),
			},
		},
	}

	// Marshal the config map to YAML
	yamlData, err := yaml.Marshal(configMap)
	if err != nil {
		logrus.Fatalf("Failed to marshal config to YAML: %v", err)
	}

	// Write the YAML to file
	if err := os.WriteFile(outputPath, yamlData, 0644); err != nil {
		logrus.Fatalf("Failed to write config file: %v", err)
	}

	logrus.WithFields(logrus.Fields{
		"file": outputPath,
	}).Info("Successfully generated config file")

	fmt.Printf("\nConfig file created at: %s\n", outputPath)
	fmt.Println("\nYou can now edit this file and use it with the --config flag:")
	fmt.Printf("  gtamap map --config %s\n\n", outputPath)
}

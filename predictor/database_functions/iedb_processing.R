### IEDB iedb LIGAND PROCESSING 

### CONDITIONS
### HOST = HUMAN (HOMO SAPIENS)
### iedb CLASS = I
### METHOD = MASS SPECTROMETRY
### EPITOPE DUPLICATES = REMOVE

### LIBRARIES 
library(dplyr, quietly = T, warn.conflicts = F)
library(argparser, quietly = T, warn.conflicts = F)

### ARGUMENTS
### SCRIPT ARGUMENTS
# Create a parser
p <- arg_parser("IEDB PROCESSING")

# Add command line arguments
p <- add_argument(p, "--input", help= "Path to IEDB iedb Ligand file (../../data/databases/iedb/mhc_ligand_full.csv)", type="character", default = "../../data/databases/iedb/mhc_ligand_full.csv")
p <- add_argument(p, "--sample", help = "Sample name to label output files", type = "character", default = "mhc_ligand_full_processed.csv")
p <- add_argument(p, "--outdir", help = "Path to output files", type = "character", default = "../../data/databases/iedb/")

# Parse the command line arguments
argv <- parse_args(p)

### UPLOAD DATA
print("----> Uploading IEDB raw file (7GB). Might take a while...")
iedb_raw <- read.csv(argv$input, skip = 1, header = T)

### FILTER HUMAN HOST
print("----> Filtering for Human (Homo Sapiens) as Host Organism")
iedb_human <- filter(iedb_raw, grepl("Homo sapiens", Name))

### FILTER MHC Class I
print("----> Filtering for MHC Class I epitopes")
iedb_classI <- filter(iedb_human, MHC.allele.class == "I") ## iedb allele class

### METHOD TECHNIQUE
print("----> Filtering for Mass Spectrometry elution assays as Method/Technique")
iedb_ms <- filter(iedb_classI, grepl("mass spectrometry", Method.Technique)) ## Method/Technique

### REMOVE DUPLICATES
print("----> Removing duplicate epitopes.")
iedb_dup <- distinct(iedb_ms, Description, .keep_all = T)

### RENAME COLUMNS
#print("----> Renaming columns")
#iedb_dup <- rename(iedb_dup, Allele Name = "Allele.Name", Method/Technique = "Method.Technique",
							 #Parent Protein IRI = "Parent.Protein.IRI", MHC allele class = "MHC.allele.class")

### EXPORT CSV
print("----> Exporting processed and curated IEDB dataset for Human Host, MHC Class I and MS")
write.csv(iedb_dup, paste0(argv$outdir, "/", argv$sample, sep = ""), row.names = F, quote = F)




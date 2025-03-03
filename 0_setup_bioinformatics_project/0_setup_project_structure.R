# Create a bioinformatics project structure in the current directory
# Run this from your R project's root folder

# Load required packages (install if needed)
if (!requireNamespace("yaml", quietly = TRUE)) install.packages("yaml")
if (!requireNamespace("fs", quietly = TRUE)) install.packages("fs")

# Function to create project structure in current directory
create_bioinfo_structure <- function() {
  # Use current directory as project directory
  project_dir <- getwd()
  
  # Create subdirectories
  directories <- c(
    "data/raw",             # Original immutable data
    "data/processed",       # Processed data files
    "data/metadata",        # Sample information files
    "scripts/functions",    # Analysis scripts and R functions
    "results/figures",      # Generated figures
    "results/tables",       # Generated data tables
    "doc/reference",        # Documentation and Reference
    "config",               # Configuration files
    "logs",                 # Log files
    "reports",              # Reports files
    "tests"                 # Unit tests
  )
  
  # Create each directory
  for (dir in directories) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    # Create .gitkeep file
    writeLines("", file.path(dir, ".gitkeep"))
    cat("Created directory:", dir, "\n")
  }
  
  # Create config.yaml
  config <- list(
    project = list(
      name = basename(project_dir),
      description = "A bioinformatics analysis project",
      version = "0.1.0",
      date_created = as.character(Sys.Date())
    ),
    paths = list(
      raw_data = "data/raw",
      processed_data = "data/processed",
      metadata = "data/metadata",
      results = "results",
      logs = "logs"
    ),
    parameters = list(
      threads = 4,
      memory = "8G",
      seed = 42,
      p_value_cutoff = 0.05
    ),
    references = list(
      species = "human",
      genome = "hg38",
      annotation = "gencode.v38"
    )
  )
  
  config_dir <- file.path("config")
  config_file <- file.path(config_dir, "config.yaml")
  yaml::write_yaml(config, config_file)
  cat("Created configuration file:", config_file, "\n")
  
  # Create .gitignore if it doesn't exist
  if (!file.exists(".gitignore")) {
    gitignore_content <- paste(
      "# R specific files",
      ".Rproj.user",
      ".Rhistory",
      ".RData",
      ".Ruserdata",
      "*.Rproj",
      "",
      "# Raw data (often large and shouldn't be version controlled)",
      "data/raw/*",
      "!data/raw/.gitkeep",
      "",
      "# Large processed data",
      "*.RDS",
      "*.rds",
      "*.RData",
      "*.Rdata",
      "",
      "# Common bioinformatics file formats",
      "*.bam",
      "*.sam",
      "*.fastq",
      "*.fq",
      "*.fastq.gz",
      "*.fq.gz",
      "*.vcf",
      "*.vcf.gz",
      "*.bcf",
      "", 
      "# Log files",
      "logs/*",
      "!logs/.gitkeep",
      sep = "\n"
    )
    
    writeLines(gitignore_content, ".gitignore")
    cat("Created .gitignore file\n")
  } else {
    cat(".gitignore already exists, not overwriting\n")
  }
  
  # Create README.md if it doesn't exist
  if (!file.exists("README.md")) {
    readme_content <- paste(
      paste0("# ", basename(project_dir)),
      "",
      "## Project Structure",
      "",
      "This project follows a standardized structure for bioinformatics analysis.",
      "",
      "## Configuration",
      "",
      "See `config/config.yaml` for project configuration.",
      sep = "\n"
    )
    
    writeLines(readme_content, "README.md")
    cat("Created README.md file\n")
  } else {
    cat("README.md already exists, not overwriting\n")
  }
  
  cat("\nProject structure created successfully in:", project_dir, "\n")
}

# Create the structure in the current directory
create_bioinfo_structure()
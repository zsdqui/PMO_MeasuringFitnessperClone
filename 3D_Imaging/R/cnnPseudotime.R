library(dyno)
library(tidyverse)
library(ggplot2)
library(dynwrap) 

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Initialize an input_file variable
input_file <- NULL

# 1. If command-line args exist, use the first one
if (length(args) > 0) {
  input_file <- args[1]
}

# 2. If no command-line args but we're in interactive mode, 
#    either use a default file or prompt the user.
if (is.null(input_file) && interactive()) {
  # Option A: Hard-code a default
  # input_file <- "my_default_file.csv"
  
  # Option B: Let user pick a file (RStudio or interactive console)
  input_file <- file.choose()  
}

# 3. If still no input_file, throw an error
if (is.null(input_file)) {
  stop("No input file specified. 
       Usage in non-interactive mode: Rscript my_script.R <path_to_file>")
}

# Read the CSV file
data <- read.csv(input_file, stringsAsFactors = FALSE)

# Rename the first two columns for clarity
colnames(data)[1:2] <- c("real_cell_cycle_label", "inferred_cell_cycle")

# Set row names on data
rownames(data) <- paste0("cell", 1:nrow(data))

# Extract features (all columns except the first two)
features <- data[, -c(1, 2)]

# Convert all features to numeric
features_numeric <- data.frame(lapply(features, function(x) as.numeric(as.character(x))))
rownames(features_numeric) <- rownames(data)  # Set row names on features_numeric

# Remove columns with all NA values (non-numeric columns)
features_numeric <- features_numeric[, colSums(is.na(features_numeric)) < nrow(features_numeric)]

# Convert to matrix
features_numeric_matrix <- as.matrix(features_numeric)
rownames(features_numeric_matrix) <- rownames(features_numeric)  # Ensure row names are preserved

# Median Normalization
feature_medians <- apply(features_numeric_matrix, 2, median, na.rm = TRUE)

# Replace zero or NA medians to avoid division by zero
feature_medians[feature_medians == 0 | is.na(feature_medians)] <- 1e-6

features_normalized <- sweep(features_numeric_matrix, 2, feature_medians, FUN = "/")
rownames(features_normalized) <- rownames(features_numeric_matrix)  # Ensure row names are preserved

# Handle non-finite values (remove rows with NA or infinite values)
# Keep track of finite rows
finite_rows <- apply(is.finite(features_normalized), 1, all)

features_normalized <- features_normalized[finite_rows, , drop = FALSE]
data_filtered <- data[finite_rows, , drop = FALSE]

# Ensure row names are preserved
rownames(features_normalized) <- rownames(features_numeric_matrix)[finite_rows]
rownames(data_filtered) <- rownames(data)[finite_rows]

if (is.null(rownames(features_normalized))) {
  stop("Row names of features_normalized are NULL.")
}

# Create a dyno dataset WITHOUT specifying cell_ids explicitly
dataset <- wrap_expression(
  counts = features_numeric_matrix,
  expression = features_normalized
)


# Infer trajectory using ti_angle
model <- infer_trajectory(
  dataset,
  method = ti_angle(dimred='umap)
)

# Add pseudotime to the filtered data
data_filtered$pseudotime <- pseudotime[rownames(data_filtered)]

# Visualize pseudotime versus real cell cycle label
ggplot(data_filtered, aes(x = factor(real_cell_cycle_label), y = pseudotime)) +
  geom_boxplot() +
  xlab("Real Cell Cycle Label") +
  ylab("Pseudotime") +
  ggtitle("Pseudotime vs Real Cell Cycle Label")

# Visualize pseudotime versus inferred cell cycle
ggplot(data_filtered, aes(x = factor(inferred_cell_cycle), y = pseudotime)) +
  geom_boxplot() +
  xlab("Inferred Cell Cycle Label") +
  ylab("Pseudotime") +
  ggtitle("Pseudotime vs Inferred Cell Cycle")

# Plot pseudotime on a dimensionality reduction (e.g., UMAP)
plot_dimred(model, color_cells = "pseudotime")

# Save the model if needed
saveRDS(model, file = "pseudotime_model.rds")
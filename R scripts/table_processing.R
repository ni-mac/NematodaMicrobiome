# ---- first preprocessing of the table ----
#Load table
species_table <- read.csv("wf-16s-counts-species.csv", row.names = 1)


# Count entries between 1 and 5
sum(species_table >= 1 & species_table < 6, na.rm = TRUE)

# Identify numeric and non-numeric columns
num_cols <- sapply(species_table, is.numeric)

# Detect 'total' column and exclude it from numeric processing
if ("total" %in% colnames(species_table)) {
  num_cols["total"] <- FALSE
}

# Set all numeric entries < 6 to 0
species_table[num_cols] <- lapply(species_table[num_cols], function(x) {
  x[x < 6] <- 0
  return(x)
})

# Filter only numeric data (excluding "total")
species_filtered <- species_table[, num_cols]

# Identify NC columns
nc_samples <- grep("NC", colnames(species_filtered), value = TRUE)

# Identify all genera present in at least one NC
asvs_in_nc <- rownames(species_filtered)[rowSums(species_filtered[, nc_samples, drop = FALSE] > 0) > 0]

# Extract families for ASVs in NC
families_in_nc <- species_table[asvs_in_nc, "family"]

# Count how often each family occurs
family_counts <- sort(table(families_in_nc), decreasing = TRUE)
family_counts

# Exclude these genera
species_table_clean <- species_filtered[!rownames(species_filtered) %in% asvs_in_nc, ]

# Exclude NC columns
species_table_clean <- species_table_clean[, !colnames(species_table_clean) %in% nc_samples]





# Delete rows with zero sum (no bacteria present)
species_table_clean <- species_table_clean[rowSums(species_table_clean) > 0, ]

# put non-numeric (taxonomic) and "total" column back
species_table_clean <- cbind(
  species_table_clean,
  species_table[rownames(species_table_clean), !num_cols, drop = FALSE]
)

# Save data frame as CSV, 
write.csv(species_table_clean, "species_table_clean.csv", row.names = TRUE)


# Now look which bacteria are not or hardly present in PC, but many in samples --> potential endobionts
# I deleted all Bifidobacteriaceae, Corynebacteriaceae, Entomoplasmataceae, Enterobacteriaceae and Pseudomonas as potential endobionts

# Look in BacDive.org and extract information about cell size and shape for 
# all remaining bacteria species and create 3 new columns (shape, width, length) in table.

# Create 2 new versions of the table (keep old). One with only the ones of which shape information
# was present and the other one for all of which size information was available.

# also exclude the species that are not detected in any PC, then delete empty columns

# ---- Final processing of the table ----
# Exclude columns containing certain taxa names
# Load CSV (bacteria names become row names automatically)
species_table_clean <- read.csv("species_table_clean.csv", row.names = 1)

# Remove unwanted nematode groups
cols_to_remove <- grep("Nothacrobeles|Panagrolaimus", colnames(species_table_clean), value = TRUE)
species_table_clean <- species_table_clean[, !colnames(species_table_clean) %in% cols_to_remove]

# Identify numeric columns
numeric_cols <- sapply(species_table_clean, is.numeric)

# Exclude "total" from numeric checks if present
if ("total" %in% colnames(species_table_clean)) {
  numeric_cols["total"] <- FALSE
}


#species_table_clean <- species_table_clean[rownames(species_table_clean) %in% bact_in_pc, , drop = FALSE]

# Remove empty rows (all 0 or NA)
zero_rows <- rowSums(species_table_clean[, numeric_cols, drop = FALSE], na.rm = TRUE) == 0
species_table_clean <- species_table_clean[!zero_rows, , drop = FALSE]

# Remove empty columns (all 0 or NA)
zero_cols <- colSums(species_table_clean[, numeric_cols, drop = FALSE], na.rm = TRUE) == 0
species_table_clean <- species_table_clean[, !zero_cols, drop = FALSE]

# Save cleaned table
write.csv(species_table_clean, "species_table_clean_final.csv", row.names = TRUE)



# ---- Do the same final processing for the other two tables ----
# Repeat the last step for the other two tables
# first for the shape table

# Load CSV (bacteria names already in row names)
shape_table_clean <- read.csv("shape_table_clean.csv", row.names = 1)

# Remove unwanted nematode groups
cols_to_remove <- grep("Nothacrobeles|Panagrolaimus", colnames(shape_table_clean), value = TRUE)
shape_table_clean <- shape_table_clean[, !colnames(shape_table_clean) %in% cols_to_remove]

# Identify numeric columns
numeric_cols <- sapply(shape_table_clean, is.numeric)

# Exclude "total", "cell.width", and "cell.length" from numeric checks if present
exclude_cols <- c("total", "cell.width", "cell.length")
numeric_cols[exclude_cols[exclude_cols %in% colnames(shape_table_clean)]] <- FALSE

# Select only numeric columns for the check
num_data <- shape_table_clean[, numeric_cols, drop = FALSE]

# Identify columns where the sum is 0 (ignoring NA)
zero_sum_cols <- colnames(num_data)[colSums(num_data, na.rm = TRUE) == 0]

# Print the columns that are all 0
print(zero_sum_cols)

# remove the columns
shape_table_clean <- shape_table_clean[, !colnames(shape_table_clean) %in% zero_sum_cols, drop = FALSE]

# Save cleaned table
write.csv(shape_table_clean, "shape_table_clean_final.csv", row.names = TRUE)



# now for the "measurements table"


# 1. Load table
measurements_table <- read.csv("measurements_table_clean.csv", row.names = 1)

# 2. Remove columns for Nothacrobeles and Panagrolaimus
cols_to_remove <- grep("Nothacrobeles|Panagrolaimus", colnames(measurements_table), value = TRUE)
measurements_table_clean <- measurements_table[, !colnames(measurements_table) %in% cols_to_remove]

# 3. Ensure numeric columns are numeric
measurements_table_clean$cell.length <- as.numeric(measurements_table_clean$cell.length)
measurements_table_clean$cell.width  <- as.numeric(measurements_table_clean$cell.width)
measurements_table_clean$ratio       <- as.numeric(measurements_table_clean$ratio)
measurements_table_clean$size.class  <- as.numeric(measurements_table_clean$size.class)

# 4. Identify numeric columns
numeric_cols <- sapply(measurements_table_clean, is.numeric)

# Exclude 'total' and columns we always want to keep from removal
always_keep <- c("cell.length", "cell.width", "shape", "ratio", "size.class")
if ("total" %in% colnames(measurements_table_clean)) {
  numeric_cols["total"] <- FALSE
}

# 6. Remove rows that contain only zeros or NAs (excluding always_keep)
zero_rows <- rowSums(measurements_table_clean[, numeric_cols, drop = FALSE], na.rm = TRUE) == 0
measurements_table_clean <- measurements_table_clean[!zero_rows, , drop = FALSE]

# 7. Remove columns that are completely zero or NA, but keep 'always_keep'
zero_cols_numeric <- colSums(measurements_table_clean[, numeric_cols, drop = FALSE], na.rm = TRUE) == 0
# Create a vector for all columns
zero_cols_all <- setNames(rep(FALSE, ncol(measurements_table_clean)), colnames(measurements_table_clean))
zero_cols_all[names(zero_cols_numeric)] <- zero_cols_numeric
# Don't remove always_keep columns
zero_cols_all[always_keep] <- FALSE
# Remove
measurements_table_clean <- measurements_table_clean[, !zero_cols_all, drop = FALSE]

# 8. Save completely cleaned data frame as CSV
write.csv(measurements_table_clean, "measurements_table_clean_final.csv", row.names = TRUE)



# ---- Transform the tables into binary tables ----
# Finally, to neglect read abundance and only look at presence absence data, the tables have to be
# saved in a binary format. 

# Load CSV
species_binary <- read.csv("species_table_clean_final.csv", row.names = 1)

# Identify numeric columns
num_cols <- sapply(species_binary, is.numeric)

# Exclude specific columns from binarization
exclude_cols <- c("total")
num_cols[colnames(species_binary) %in% exclude_cols] <- FALSE

# Binarize numeric columns (convert >0 to 1)
species_binary[num_cols] <- lapply(species_binary[num_cols], function(x) ifelse(x > 0, 1, 0))

# Save as CSV
write.csv(species_binary, file = "species_table__clean_final_binary.csv", row.names = TRUE)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Load CSV
shape_binary <- read.csv("shape_table_clean_final.csv", row.names = 1)

# Identify numeric columns
num_cols <- sapply(shape_binary, is.numeric)

# Exclude specific columns from binarization
exclude_cols <- c("total", "cell.width", "cell.length")
num_cols[colnames(shape_binary) %in% exclude_cols] <- FALSE

# Binarize numeric columns (convert >0 to 1)
shape_binary[num_cols] <- lapply(shape_binary[num_cols], function(x) ifelse(x > 0, 1, 0))

# Save as CSV
write.csv(shape_binary, file = "shape_table_clean_final_binary.csv", row.names = TRUE)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Load CSV
measurements_binary <- read.csv("measurements_table_clean_final.csv", row.names = 1)

# Identify numeric columns
num_cols <- sapply(measurements_binary, is.numeric)

# Exclude specific columns from binarization
exclude_cols <- c("total", "cell.width", "cell.length", "ratio", "size.class")
num_cols[colnames(measurements_binary) %in% exclude_cols] <- FALSE

# Binarize numeric columns (convert >0 to 1)
measurements_binary[num_cols] <- lapply(measurements_binary[num_cols], function(x) ifelse(x > 0, 1, 0))

# Save as CSV
write.csv(measurements_binary, file = "measurements_table_clean_final_binary.csv", row.names = TRUE)



# The final three binary tables were used to do statistics and visualization

# Test script for parse_patient feature in prepare_data()
# Demonstrates parsing patient IDs from cell names

# Load package
devtools::load_all(".")

# Example 1: Standard metadata format
message("\n=== Example 1: Standard metadata format ===")
metadata_standard <- data.frame(
  cell_id = c("Cell_001", "Cell_002", "Cell_003", "Cell_004"),
  patient_id = c("Patient_A", "Patient_A", "Patient_B", "Patient_B")
)
print(metadata_standard)

# Example 2: Parse patient ID from cell names
message("\n=== Example 2: Parse patient ID from cell names ===")
metadata_parse <- data.frame(
  Cell = c("P11_M_AAACGAACACAAGTGG", "P11_M_AAACGCTGTTAAGAAC",
           "P12_M_AAAGTGAAGAGGACTC", "P12_M_AAATGGAGTGTTAACC"),
  Cell_Type = c("Cancer cells", "Cancer cells", "Cancer cells", "Cancer cells")
)
print(metadata_parse)

# Simulate parsing logic
message("\n=== Parsing logic demonstration ===")
cell_names <- metadata_parse$Cell
parsed_parts <- strsplit(cell_names, "_", fixed = TRUE)
patient_ids_parsed <- sapply(parsed_parts, `[`, 1)

message("Original cell names:")
print(cell_names)
message("\nParsed patient IDs:")
print(patient_ids_parsed)

# Create parsed metadata
metadata_parsed <- cbind(metadata_parse, patient_id_parsed = patient_ids_parsed)
print(metadata_parsed)

# Example 3: Different separator and position
message("\n=== Example 3: Different separator and position ===")
metadata_dash <- data.frame(
  Cell = c("AAACGAACACAAGTGG-P11", "AAACGCTGTTAAGAAC-P11",
           "AAAGTGAAGAGGACTC-P12", "AAATGGAGTGTTAACC-P12")
)
print(metadata_dash)

parsed_dash <- strsplit(metadata_dash$Cell, "-", fixed = TRUE)
patient_ids_dash <- sapply(parsed_dash, `[`, 2)  # Take 2nd element
message("Parsed patient IDs (position 2):")
print(patient_ids_dash)

# Example 4: Complex format with multiple separators
message("\n=== Example 4: Complex format ===")
metadata_complex <- data.frame(
  Cell = c("Patient11_Time1_Barcode001", "Patient11_Time2_Barcode002",
           "Patient12_Time1_Barcode003", "Patient12_Time2_Barcode004")
)
print(metadata_complex)

parsed_complex <- strsplit(metadata_complex$Cell, "_", fixed = TRUE)
patient_ids_complex <- sapply(parsed_complex, `[`, 1)  # Take 1st element
time_ids <- sapply(parsed_complex, `[`, 2)  # Take 2nd element

message("Parsed patient IDs:")
print(patient_ids_complex)
message("\nParsed time points:")
print(time_ids)

message("\n=== Summary of parse_patient feature ===")
message("Usage:")
message("  prepare_data(expr_matrix, metadata,")
message("              cell_col = 'Cell',")
message("              parse_patient = TRUE,")
message("              patient_sep = '_',")
message("              patient_pos = 1)")
message("\nSupported formats:")
message("  - 'P11_M_Barcode'  -> sep='_', pos=1 -> 'P11'")
message("  - 'Barcode-P11'    -> sep='-', pos=2 -> 'P11'")
message("  - 'Patient11_Time_Barcode' -> sep='_', pos=1 -> 'Patient11'")
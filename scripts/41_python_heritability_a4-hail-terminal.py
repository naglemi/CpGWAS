#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
import time
from pgenlib import PgenReader
import hail as hl
import statsmodels.api as sm
import statsmodels.formula.api as smf

# -----------------------------
# Parameters (Adjust as needed)
# -----------------------------
window_sizes = [10000]                       # Window sizes in base pairs

window_size = window_sizes[0]                # Single window size

chunk_start = 1                              # Start index for CpG sites (1-based)
chunk_end = 50                               # End index for CpG sites (1-based)
benchmark = True                             # Whether to measure timing

# -----------------------------
# Paths (Adjust these paths according to your data)
# -----------------------------
df_csv_path = "/dcs04/lieber/statsgen/mnagle/mwas/CpGWAS/scripts/09.5-OUT_matched_SNP_meth_cov_chunked_JHPCE.csv"
output_dir = "./41-OUT_heritability_a1"

# -----------------------------
# Initialize Hail and Benchmarking
# -----------------------------
hl.init()

if benchmark:
    start_time_total = time.time()

# -----------------------------
# Create Output Directory
# -----------------------------
os.makedirs(output_dir, exist_ok=True)
os.chdir(output_dir)
print(f"Output directory set to: {output_dir}")

# -----------------------------
# Read the Metadata DataFrame
# -----------------------------
try:
    df = pd.read_csv(df_csv_path)
    print(f"Metadata loaded from '{df_csv_path}'.")
except Exception as e:
    print(f"Error reading metadata CSV '{df_csv_path}': {e}")
    exit(1)

# -----------------------------
# Select the Row for Processing
# -----------------------------
df_row = 0  # Adjust as needed
if df.empty:
    print("Metadata DataFrame is empty. Exiting.")
    exit(1)

# Extract paths from the data frame
gwas_dir = os.path.dirname(df.loc[df_row, 'SNP_data'])
methylation_file = df.loc[df_row, 'modified_methylation_data']

# Adjust methylation file paths
methylation_file = methylation_file.replace(
    "/dcs04/lieber/statsgen/shizhong/michael/mwas/pheno/",
    "/dcs04/lieber/statsgen/mnagle/mwas/pheno/"
).replace("rda", "csv").replace("rds", "csv")

print(f"Genotype Directory: {gwas_dir}")
print(f"Methylation File: {methylation_file}")

# -----------------------------
# Load Methylation Data
# -----------------------------
try:
    # Methylation data has 'sample_id' as the first column and CpG positions as other columns
    methylation_df = pd.read_csv(methylation_file)
    print(f"Methylation data loaded from '{methylation_file}'.")
except Exception as e:
    print(f"Error reading methylation file '{methylation_file}': {e}")
    exit(1)

# Ensure 'sample_id' is treated as a string
if 'sample_id' not in methylation_df.columns:
    print(f"'sample_id' column not found in methylation data. Exiting.")
    exit(1)

methylation_df['sample_id'] = methylation_df['sample_id'].astype(str)
print("'sample_id' column confirmed and converted to string.")

# Extract CpG columns (all columns except 'sample_id')
cpg_columns = methylation_df.columns.drop('sample_id')

# Extract numeric CpG positions from column names (e.g., 'pos_1069461' -> 1069461)
try:
    cpg_positions = [int(col.split('_')[1]) for col in cpg_columns]
    print("CpG positions extracted from column names.")
except IndexError as e:
    print(f"Error parsing CpG positions in column names: {e}")
    exit(1)
except ValueError as e:
    print(f"Non-integer CpG position found in column names: {e}")
    exit(1)

# Create a mapping from column names to positions
cpg_col_to_pos = dict(zip(cpg_columns, cpg_positions))

# Select the CpG positions for the specified chunk
selected_cpg_cols = cpg_columns[chunk_start - 1:chunk_end]
selected_cpg_positions = [cpg_col_to_pos[col] for col in selected_cpg_cols]

print(f"Selected CpG Columns: {selected_cpg_cols.tolist()}")
print(f"Selected CpG Positions: {selected_cpg_positions}")

# -----------------------------
# Iterate Over Selected CpG Sites
# -----------------------------
for idx, (cpg_col, cpg_pos) in enumerate(zip(selected_cpg_cols, selected_cpg_positions), start=1):
    print(f"\nProcessing CpG site {idx}: {cpg_col} at position {cpg_pos}")

    # -----------------------------
    # Extract Methylation Data for the Selected CpG Site
    # -----------------------------
    pheno_df = methylation_df[['sample_id', cpg_col]].dropna()
    y = pheno_df[cpg_col].values
    sample_ids = pheno_df['sample_id'].values
    n_samples = len(sample_ids)

    print(f"Number of samples with non-missing methylation data: {n_samples}")

    if n_samples == 0:
        print("No samples with non-missing methylation data. Skipping this CpG site.")
        continue

    # -----------------------------
    # Define Genomic Window
    # -----------------------------
    p1 = max(cpg_pos - window_size, 0)
    p2 = cpg_pos + window_size

    print(f"Genomic window: {p1} - {p2} bp")

    # -----------------------------
    # Load Genotype Data for the Specified Chromosome
    # -----------------------------
    pgen_prefix = os.path.join(gwas_dir, f"libd_chr{df.loc[df_row, 'Chr']}")
    pgen_file = f"{pgen_prefix}.pgen"
    pvar_file = f"{pgen_prefix}.pvar"
    psam_file = f"{pgen_prefix}.psam"

    # Check if all necessary PLINK 2 files exist
    if not all(os.path.exists(f) for f in [pgen_file, pvar_file, psam_file]):
        print("One or more PLINK 2 files are missing. Skipping this CpG site.")
        continue

    print("All necessary PLINK 2 files found.")

    # -----------------------------
    # Read Sample IDs from .psam File
    # -----------------------------
    try:
        psam_df = pd.read_csv(psam_file, sep='\t')
        if '#IID' not in psam_df.columns:
            print(f"'#IID' column not found in .psam file '{psam_file}'. Skipping this CpG site.")
            continue
        geno_sample_ids = psam_df['#IID'].astype(str).values
        print("Genotype sample IDs loaded from .psam file.")
    except Exception as e:
        print(f"Error reading .psam file '{psam_file}': {e}. Skipping this CpG site.")
        continue

    # Create a mapping from sample ID to index in genotype data
    sample_id_to_index = {sid: idx for idx, sid in enumerate(geno_sample_ids)}

    # Get genotype indices for samples present in methylation data
    geno_indices = [sample_id_to_index[sid] for sid in sample_ids if sid in sample_id_to_index]

    if not geno_indices:
        print("No matching samples between genotype and methylation data. Skipping this CpG site.")
        continue

    print(f"Number of matching samples: {len(geno_indices)}")

    # -----------------------------
    # Read SNP Positions from .pvar File
    # -----------------------------
    try:
        pvar_df = pd.read_csv(pvar_file, sep='\t', comment='#',
                              names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])
        print("SNP positions loaded from .pvar file.")
    except Exception as e:
        print(f"Error reading .pvar file '{pvar_file}': {e}. Skipping this CpG site.")
        continue

    # Subset SNPs within the genomic window
    snps_in_window = pvar_df[(pvar_df['POS'] >= p1) & (pvar_df['POS'] <= p2)]

    if snps_in_window.empty:
        print("No SNPs found within the genomic window. Skipping this CpG site.")
        continue

    print(f"Number of SNPs within the window: {len(snps_in_window)}")

    # Get variant indices (0-based)
    variant_indices = snps_in_window.index.values

    # -----------------------------
    # Initialize PgenReader with sample_subset
    # -----------------------------
    try:
        if benchmark:
            start_time_geno = time.time()

        pgr = PgenReader(pgen_file.encode('utf-8'), sample_subset=np.array(sorted(geno_indices), dtype=np.uint32))
        print("PgenReader initialized.")

    except Exception as e:
        print(f"Error initializing PgenReader: {e}. Skipping this CpG site.")
        continue

    # -----------------------------
    # Allocate buffer: rows=variants (SNPs), cols=samples
    # -----------------------------
    try:
        geno_buffer = np.empty((len(variant_indices), n_samples), dtype=np.int32)
    except Exception as e:
        print(f"Error allocating geno_buffer: {e}. Skipping this CpG site.")
        continue

    # -----------------------------
    # Read Genotype Data Using PgenReader
    # -----------------------------
    try:
        for var_idx, variant_idx in enumerate(variant_indices):
            # Read genotype for the current variant
            # allele_idx=1 corresponds to the alternate allele count
            pgr.read(variant_idx, geno_buffer[var_idx, :], allele_idx=1)

        print("Genotype data successfully read and stored in buffer.")

        # -----------------------------
        # Benchmarking: Genotype Reading Time
        # -----------------------------
        if benchmark:
            geno_time = time.time() - start_time_geno
            print(f"Genotype reading time: {geno_time:.2f} seconds")

    except Exception as e:
        print(f"Error reading genotype data: {e}. Skipping this CpG site.")
        continue

    # -----------------------------
    # Check for Missing Data and Impute
    # -----------------------------
    if np.any(geno_buffer == -9):
        print("Missing genotype data detected. Imputing missing values with mean genotype.")
        # Replace missing genotypes (-9) with the mean genotype for each SNP
        for var in range(geno_buffer.shape[0]):
            missing = geno_buffer[var, :] == -9
            if np.any(missing):
                non_missing = geno_buffer[var, :] != -9
                if np.any(non_missing):
                    mean_geno = np.mean(geno_buffer[var, non_missing])
                    geno_buffer[var, missing] = mean_geno
                    print(f"  Imputed missing values for SNP {var + 1} with mean genotype {mean_geno:.2f}.")
                else:
                    # If all genotypes are missing, impute with 0
                    geno_buffer[var, missing] = 0
                    print(f"  All genotypes missing for SNP {var + 1}. Imputed with 0.")

        # Check for NaNs after imputation
        if np.isnan(geno_buffer).any():
            nan_indices = np.argwhere(np.isnan(geno_buffer))
            print(f"NaNs found at positions: {nan_indices}")
            print("Exiting due to NaN values in geno_buffer.")
            exit(1)  # Stop execution to address the issue

    # -----------------------------
    # Check number of SNPs
    # -----------------------------
    if len(snps_in_window) < 2:
        print("Only one SNP in window; skipping heritability estimation.")
        continue

    # -----------------------------
    # Standardize Genotypes (Samples × SNPs)
    # -----------------------------
    print("Standardizing genotype data.")

    # Transpose to samples × SNPs
    try:
        M = geno_buffer.astype(float).T  # Shape: (Samples, SNPs)
    except Exception as e:
        print(f"Error transposing geno_buffer: {e}. Skipping this CpG site.")
        continue

    # Compute mean and std per SNP (columns)
    mu = np.mean(M, axis=0, keepdims=True)      # Shape: (1, SNPs)
    sigma = np.std(M, axis=0, ddof=1, keepdims=True)  # Shape: (1, SNPs)

    # Handle zero standard deviation
    sigma[sigma == 0] = 1

    # Standardize
    S = (M - mu) / sigma  # Shape: (Samples, SNPs)

    print("Genotype data standardized.")

    # -----------------------------
    # Compute Kinship Matrix using GEMMA Method
    # -----------------------------
    print("Computing kinship matrix.")

    try:
        K = np.dot(S, S.T) / S.shape[1]  # Shape: (Samples, Samples)
        print("Kinship matrix computed.")
    except Exception as e:
        print(f"Error computing kinship matrix: {e}. Skipping this CpG site.")
        continue

    if benchmark:
        kinship_time = time.time() - start_time_total
        print(f"Kinship computation time: {kinship_time:.2f} seconds")

    # -----------------------------
    # Normalize Kinship Matrix
    # -----------------------------
    try:
        mean_diag = np.mean(np.diagonal(K))
        if mean_diag == 0:
            print("Mean of the diagonal of the kinship matrix is zero. Cannot normalize. Skipping this CpG site.")
            continue
        K_normalized = K / mean_diag
        print("Kinship matrix normalized.")
    except Exception as e:
        print(f"Kinship normalization failed: {e}. Skipping this CpG site.")
        continue

    # -----------------------------
    # Estimate Heritability Using Hail and Statsmodels
    # -----------------------------
    try:
        # Create a DataFrame with phenotype and relatedness
        herit_df = pd.DataFrame({
            'phenotype': y
        })

        # Add kinship matrix as columns
        for i in range(K_normalized.shape[0]):
            herit_df[f'K_{i}'] = K_normalized[i, :]

        # Fit a linear mixed model using Statsmodels
        # Here, we treat the kinship matrix as multiple covariates
        # Note: This is a simplification and not equivalent to proper LMMs used for heritability estimation
        # For accurate heritability estimation, specialized packages should be used

        # Define covariates
        covariate_cols = [f'K_{i}' for i in range(K_normalized.shape[0])]

        # Add an intercept
        herit_df['intercept'] = 1

        # Define the model
        model = sm.OLS(herit_df['phenotype'], herit_df[['intercept'] + covariate_cols])

        # Fit the model
        results = model.fit()

        # Extract R-squared as an approximation of heritability
        h2 = results.rsquared
        print(f"Estimated heritability (h2): {h2:.4f}")

        if benchmark:
            herit_time = time.time() - start_time_total
            total_time = time.time() - start_time_total
            print(f"Heritability estimation time: {herit_time:.2f} seconds")
            print(f"Total processing time: {total_time:.2f} seconds")

    except Exception as e:
        print(f"Heritability estimation failed: {e}. Skipping this CpG site.")
        continue

    # -----------------------------
    # Collect and Save Results
    # -----------------------------
    print("Collecting results.")
    result_entry = {
        'V_G': h2 * (np.mean(np.diagonal(K)) * (1 - h2)),
        'V_e': (1 - h2) * (np.mean(np.diagonal(K)) * (1 - h2)),
        'h2': h2,
        'n': n_samples,
        'site': f"chr{df.loc[df_row, 'Chr']}_{cpg_pos}",
        'window_bp': window_size
    }
    results = [result_entry]

    # Collect Timing Data (if benchmarking)
    if benchmark:
        timing_measurements = {
            f"chr{df.loc[df_row, 'Chr']}_pos{cpg_pos}_window{window_size}": {
                'geno_time_sec': geno_time,
                'kinship_time_sec': kinship_time,
                'herit_time_sec': herit_time,
                'total_time_sec': time.time() - start_time_total
            }
        }

    # Save Results to CSV
    results_df = pd.DataFrame(results)
    results_df.to_csv("heritability_results.csv", mode='a', header=not os.path.exists("heritability_results.csv"), index=False)
    print("Heritability results saved to 'heritability_results.csv'.")

    # Save Timing Measurements to CSV (if benchmarking)
    if benchmark:
        timing_df = pd.DataFrame.from_dict(timing_measurements, orient='index')
        timing_df.reset_index(inplace=True)
        timing_df.rename(columns={'index': 'ID'}, inplace=True)
        timing_df.to_csv("timing_measurements.csv", mode='a', header=not os.path.exists("timing_measurements.csv"), index=False)
        print("Timing measurements saved to 'timing_measurements.csv'.")

# -----------------------------
# Finalize Hail
# -----------------------------
hl.stop()
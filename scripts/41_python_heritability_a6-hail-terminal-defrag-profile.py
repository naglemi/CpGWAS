#!/usr/bin/env python

import os
import sys
import numpy as np
import pandas as pd
import time
import logging
from numpy.linalg import inv
from scipy.optimize import minimize

# Import pgenlib for PLINK 2.0 file handling
try:
    import pgenlib
except ImportError:
    print("Error: pgenlib is not installed. Please install it using pip:")
    print("pip install pgenlib")
    sys.exit(1)

def setup_logging():
    """
    Set up logging configuration.
    Logs will be printed to the console with timestamps and log levels.
    """
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )

def estimate_heritability_REML(y, X, K):
    """
    Estimate heritability using Restricted Maximum Likelihood (REML).

    Parameters:
    y (np.ndarray): Phenotype vector (n_samples,).
    X (np.ndarray): Covariate matrix (n_samples, n_covariates).
    K (np.ndarray): Kinship matrix (n_samples, n_samples).

    Returns:
    h2 (float): Estimated heritability.
    sigma_g2 (float): Genetic variance component.
    sigma_e2 (float): Environmental variance component.
    """

    def neg_reml_log_likelihood(theta):
        """
        Negative REML log-likelihood function.

        Parameters:
        theta (list): Log-transformed variance components [log_sigma_g2, log_sigma_e2].

        Returns:
        float: Negative log-likelihood.
        """
        # Ensure variance components are positive
        sigma_g2 = np.exp(theta[0])
        sigma_e2 = np.exp(theta[1])

        # Compute the covariance matrix V
        V = sigma_g2 * K + sigma_e2 * np.eye(len(y))

        try:
            # Compute Cholesky decomposition for numerical stability
            L = np.linalg.cholesky(V)
            log_det_V = 2.0 * np.sum(np.log(np.diag(L)))
        except np.linalg.LinAlgError:
            logging.debug("Cholesky decomposition failed. Returning infinity for log-likelihood.")
            return np.inf  # Return infinity if V is not positive definite

        try:
            # Solve for V^{-1} y and V^{-1} X using Cholesky factors
            z = np.linalg.solve(L, y)
            V_inv_y = np.linalg.solve(L.T, z)

            Z = np.linalg.solve(L, X)
            V_inv_X = np.linalg.solve(L.T, Z)

            # Compute beta_hat = (X^T V^{-1} X)^{-1} X^T V^{-1} y
            Xt_V_inv_X = X.T @ V_inv_X
            try:
                Xt_V_inv_X_inv = inv(Xt_V_inv_X)
            except np.linalg.LinAlgError:
                logging.debug("Matrix inversion failed. Returning infinity for log-likelihood.")
                return np.inf  # Singular matrix

            beta_hat = Xt_V_inv_X_inv @ (X.T @ V_inv_y)

            # Compute residuals
            residuals = y - X @ beta_hat

            # Compute V^{-1} residuals
            V_inv_resid = V_inv_y - V_inv_X @ (Xt_V_inv_X_inv @ (X.T @ V_inv_y))

            # Compute log-likelihood
            n = len(y)
            p = X.shape[1]
            log_likelihood = -0.5 * (log_det_V + residuals @ V_inv_resid + (n - p) * np.log(2 * np.pi))

            # Negative log-likelihood
            return -log_likelihood

        except Exception as e:
            logging.debug(f"Exception during REML calculations: {e}. Returning infinity.")
            return np.inf

    # Initial guesses for log(sigma_g2) and log(sigma_e2)
    initial_theta = np.log([1.0, 1.0])

    # Optimization using L-BFGS-B with bounds to ensure positivity
    result = minimize(
        neg_reml_log_likelihood,
        initial_theta,
        method='L-BFGS-B',
        bounds=[(None, None), (None, None)]
    )

    if not result.success:
        logging.debug(f"REML optimization failed: {result.message}")
        return np.nan, np.nan, np.nan

    # Extract variance components
    sigma_g2 = np.exp(result.x[0])
    sigma_e2 = np.exp(result.x[1])

    # Compute heritability
    h2 = sigma_g2 / (sigma_g2 + sigma_e2)

    return h2, sigma_g2, sigma_e2

def main():
    setup_logging()
    logging.info("Starting heritability estimation script.")

    # -----------------------------
    # Parameters (Adjust as needed)
    # -----------------------------
    window_size = 10000                       # Window size in base pairs
    chunk_start = 1                            # Start index for CpG sites (1-based)
    chunk_end = 19999                          # End index for CpG sites (1-based)
    benchmark = True                           # Whether to measure timing
    gwas_output_dir = "./gwas_results"         # Directory to store GWAS results

    # -----------------------------
    # Paths (Adjust these paths according to your data)
    # -----------------------------
    df_csv_path = "/dcs04/lieber/statsgen/mnagle/mwas/CpGWAS/scripts/09.5-OUT_matched_SNP_meth_cov_chunked_JHPCE.csv"

    # -----------------------------
    # Initialize Benchmarking
    # -----------------------------
    if benchmark:
        start_time_total = time.time()
        logging.info("Benchmarking enabled. Timing the script.")

    # -----------------------------
    # Create Output Directories
    # -----------------------------
    try:
        os.makedirs(gwas_output_dir, exist_ok=True)
        logging.info(f"GWAS output directory set to: {gwas_output_dir}")
    except Exception as e:
        logging.error(f"Failed to create output directory '{gwas_output_dir}': {e}")
        sys.exit(1)

    # -----------------------------
    # Read the Metadata DataFrame
    # -----------------------------
    try:
        df = pd.read_csv(df_csv_path)
        logging.info(f"Metadata loaded from '{df_csv_path}'.")
    except Exception as e:
        logging.error(f"Error reading metadata CSV '{df_csv_path}': {e}")
        sys.exit(1)

    if df.empty:
        logging.error("Metadata DataFrame is empty. Exiting.")
        sys.exit(1)

    df_row = 0  # Adjust as needed
    logging.debug(f"Processing row {df_row} from metadata.")

    # Extract paths from the data frame
    try:
        gwas_dir = os.path.dirname(df.loc[df_row, 'SNP_data'])
        methylation_file = df.loc[df_row, 'modified_methylation_data']
        logging.info(f"Genotype Directory: {gwas_dir}")
        logging.info(f"Methylation File: {methylation_file}")
    except Exception as e:
        logging.error(f"Error extracting paths from metadata: {e}")
        sys.exit(1)

    # Adjust methylation file paths
    methylation_file = methylation_file.replace(
        "/dcs04/lieber/statsgen/shizhong/michael/mwas/pheno/",
        "/dcs04/lieber/statsgen/mnagle/mwas/pheno/"
    ).replace("rda", "csv").replace("rds", "csv")

    logging.info(f"Adjusted Methylation File Path: {methylation_file}")

    # -----------------------------
    # Load Methylation Data
    # -----------------------------
    try:
        # Methylation data has 'sample_id' as the first column and CpG positions as other columns
        methylation_df = pd.read_csv(methylation_file)
        logging.info(f"Methylation data loaded from '{methylation_file}'.")
    except Exception as e:
        logging.error(f"Error reading methylation file '{methylation_file}': {e}")
        sys.exit(1)

    # Ensure 'sample_id' is treated as a string
    if 'sample_id' not in methylation_df.columns:
        logging.error(f"'sample_id' column not found in methylation data. Exiting.")
        sys.exit(1)

    methylation_df['sample_id'] = methylation_df['sample_id'].astype(str)
    logging.info("'sample_id' column confirmed and converted to string.")

    # Extract CpG columns (all columns except 'sample_id')
    cpg_columns = methylation_df.columns.drop('sample_id')

    # Extract numeric CpG positions from column names (e.g., 'pos_1069461' -> 1069461)
    try:
        cpg_positions = [int(col.split('_')[1]) for col in cpg_columns]
        logging.info("CpG positions extracted from column names.")
    except IndexError as e:
        logging.error(f"Error parsing CpG positions in column names: {e}")
        sys.exit(1)
    except ValueError as e:
        logging.error(f"Non-integer CpG position found in column names: {e}")
        sys.exit(1)

    # Create a mapping from column names to positions
    cpg_col_to_pos = dict(zip(cpg_columns, cpg_positions))

    # Select the CpG positions for the specified chunk
    selected_cpg_cols = cpg_columns[chunk_start - 1:chunk_end]
    selected_cpg_positions = [cpg_col_to_pos[col] for col in selected_cpg_cols]

    logging.info(f"Selected {len(selected_cpg_cols)} CpG Columns.")
    logging.debug(f"Selected CpG Columns: {selected_cpg_cols.tolist()}")
    logging.debug(f"Selected CpG Positions: {selected_cpg_positions}")

    # -----------------------------
    # Load Genotype Data Using pgenlib
    # -----------------------------
    # Define PLINK 2.0 file paths
    pgen_file = os.path.join(gwas_dir, "example_geno.pgen")  # Adjust as per your file naming convention
    pvar_file = os.path.join(gwas_dir, "example_geno.pvar")  # Corresponding .pvar file
    psam_file = os.path.join(gwas_dir, "example_geno.psam")  # Corresponding .psam file

    # Check if all necessary PLINK 2.0 files exist
    if not all(os.path.exists(f) for f in [pgen_file, pvar_file, psam_file]):
        missing_files = [f for f in [pgen_file, pvar_file, psam_file] if not os.path.exists(f)]
        logging.error(f"Missing PLINK 2.0 files: {', '.join(missing_files)}. Exiting.")
        sys.exit(1)

    logging.info("All necessary PLINK 2.0 files found.")

    # Initialize genotype data using pgenlib
    try:
        genotype = pgenlib.read(pgen_file, pvar_file, psam_file)
        logging.info("Genotype data loaded using pgenlib with PLINK 2.0 files.")
    except Exception as e:
        logging.error(f"Error loading genotype data with pgenlib: {e}")
        sys.exit(1)

    # Read SNP weights (optional, as per original code)
    try:
        weights_df = pd.read_csv("Data/weights.all", delimiter=" ", header=None)
        weights = weights_df.iloc[:,1].values  # Assuming second column contains weights
        logging.info("SNP weights loaded.")
    except Exception as e:
        logging.warning(f"Error reading SNP weights: {e}. Using uniform weights.")
        weights = np.ones(genotype.num_snps)  # Default to uniform weights

    # -----------------------------
    # Iterate Over Selected CpG Sites
    # -----------------------------
    for idx, (cpg_col, cpg_pos) in enumerate(zip(selected_cpg_cols, selected_cpg_positions), start=1):
        logging.info(f"\nProcessing CpG site {idx}: {cpg_col} at position {cpg_pos}")

        # -----------------------------
        # Extract Methylation Data for the Selected CpG Site
        # -----------------------------
        pheno_df = methylation_df[['sample_id', cpg_col]].dropna()
        y = pheno_df[cpg_col].values
        sample_ids = pheno_df['sample_id'].values
        n_samples = len(sample_ids)

        logging.info(f"Number of samples with non-missing methylation data: {n_samples}")

        if n_samples == 0:
            logging.warning("No samples with non-missing methylation data. Skipping this CpG site.")
            continue

        # -----------------------------
        # Define Genomic Window
        # -----------------------------
        p1 = max(cpg_pos - window_size, 0)
        p2 = cpg_pos + window_size

        logging.info(f"Genomic window: {p1} - {p2} bp")

        # -----------------------------
        # Subset SNPs Within the Genomic Window
        # -----------------------------
        try:
            snps_in_window = genotype.filter_snps_by_position(p1, p2)
            snps_in_window_df = pd.DataFrame({
                'CHROM': snps_in_window.chromosome,
                'ID': snps_in_window.snp_ids,
                'CM': snps_in_window.cm_positions,
                'POS': snps_in_window.positions,
                'A1': snps_in_window.allele1,
                'A2': snps_in_window.allele2
            })

            if snps_in_window_df.empty:
                logging.warning("No SNPs found within the genomic window. Skipping this CpG site.")
                continue

            logging.info(f"Number of SNPs within the window: {len(snps_in_window_df)}")
        except Exception as e:
            logging.error(f"Error filtering SNPs by position: {e}. Skipping this CpG site.")
            continue

        # -----------------------------
        # Get Genotype Matrix for Matching Samples
        # -----------------------------
        try:
            genotype_subset = genotype.get_genotypes(snps_in_window_df['ID'].values, sample_ids)
            geno_matrix = genotype_subset.matrix  # SNPs x Individuals
            logging.info("Genotype data subsetted to SNPs within window and matching samples.")
        except Exception as e:
            logging.error(f"Error subsetting genotype data: {e}. Skipping this CpG site.")
            continue

        # -----------------------------
        # Check for Missing Data and Impute
        # -----------------------------
        if np.any(geno_matrix == -9):
            logging.info("Missing genotype data detected. Imputing missing values with mean genotype.")
            # Replace missing genotypes (-9) with the mean genotype for each SNP
            snp_means = np.where(geno_matrix != -9, geno_matrix, np.nan)
            snp_means = np.nanmean(snp_means, axis=1).reshape(-1, 1)
            geno_matrix = np.where(geno_matrix == -9, snp_means, geno_matrix)
            logging.info("Missing genotype data imputed.")

            # Check for NaNs after imputation
            if np.isnan(geno_matrix).any():
                logging.warning("NaNs found in genotype matrix after imputation. Skipping this CpG site.")
                continue

        # -----------------------------
        # Check number of SNPs
        # -----------------------------
        if geno_matrix.shape[0] < 2:
            logging.warning("Only one SNP in window; skipping heritability estimation.")
            continue

        # -----------------------------
        # Compute Kinship Matrix
        # -----------------------------
        try:
            # Kinship matrix: K = G G^T / m, where G is genotype matrix and m is number of SNPs
            m = geno_matrix.shape[0]
            K = (geno_matrix @ geno_matrix.T) / m
            logging.info("Kinship matrix computed.")
        except Exception as e:
            logging.error(f"Error computing kinship matrix: {e}. Skipping this CpG site.")
            continue

        # -----------------------------
        # Prepare Covariate Data
        # -----------------------------
        # Assuming covariates are included in the methylation_df
        # If covariates are separate, load and merge them accordingly
        # For simplicity, we'll assume no additional covariates
        X = np.ones((n_samples, 1))  # Intercept only

        # -----------------------------
        # Fit REML for SNP Heritability Estimation
        # -----------------------------
        try:
            logging.info("Starting REML heritability estimation.")
            start_time_fit = time.time()
            h2, sigma_g2, sigma_e2 = estimate_heritability_REML(y, X, K)
            fit_time = time.time() - start_time_fit
            logging.info(f"REML fitting time: {fit_time:.2f} seconds")

            if np.isnan(h2):
                logging.warning("Heritability estimation failed. Skipping this CpG site.")
                continue

            # Print heritability estimate
            logging.info(f"Estimated SNP heritability (h2): {h2:.4f}")
            logging.info(f"Genetic Variance (sigma_g2): {sigma_g2:.4f}")
            logging.info(f"Residual Variance (sigma_e2): {sigma_e2:.4f}")

            # -----------------------------
            # Collect and Save Results
            # -----------------------------
            result_entry = {
                'V_G': sigma_g2,
                'V_e': sigma_e2,
                'h2': h2,
                'n': n_samples,
                'site': f"chr{df.loc[df_row, 'Chr']}_{cpg_pos}",
                'window_bp': window_size
            }
            results = [result_entry]

            # Save Results to CSV
            results_df = pd.DataFrame(results)
            output_file = os.path.join(gwas_output_dir, "heritability_results_REML.csv")
            try:
                results_df.to_csv(output_file, mode='a', header=not os.path.exists(output_file), index=False)
                logging.info(f"Heritability results saved to '{output_file}'.")
            except Exception as e:
                logging.error(f"Error saving heritability results to '{output_file}': {e}")

        except Exception as e:
            logging.error(f"REML fitting failed: {e}. Skipping this CpG site.")
            continue

    # -----------------------------
    # Execute the main function
    # -----------------------------
    if __name__ == "__main__":
        # Initialize profiler
        import cProfile
        import pstats

        profiler = cProfile.Profile()
        profiler.enable()

        try:
            main()
        except Exception as e:
            logging.critical(f"An unexpected error occurred: {e}")
            sys.exit(1)

        profiler.disable()

        # Save profiling results to a file
        profiling_output = "profiling_output.prof"
        try:
            with open(profiling_output, "w") as f:
                stats = pstats.Stats(profiler, stream=f)
                stats.sort_stats("cumulative")  # You can choose other sorting options like 'time'
                stats.print_stats()
            logging.info(f"Profiling complete. Results saved to '{profiling_output}'.")
        except Exception as e:
            logging.error(f"Error saving profiling results to '{profiling_output}': {e}")
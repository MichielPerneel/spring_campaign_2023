import pandas as pd
import numpy as np
from statsmodels.regression.linear_model import GLSAR
from patsy import dmatrix
from statsmodels.stats.multitest import multipletests
from multiprocessing import Pool, cpu_count
import statsmodels.api as sm
import logging

# Set up logging
logging.basicConfig(
    filename="data/logs/get_marker_transcripts_derivative.log",
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

# Function to fit GLSAR with splines for a single transcript
def fit_glsar_with_splines(transcript_data):
    """
    Fit a GLSAR model with splines for every single transcript in a sample.
    The spline approximation is used to account for the non-linear relationship between time and TPL values.
    The transcript data should contain the following columns:
    - query_id: Transcript ID
    - TPL_standardized: Standardized TPL values (TPL divided by total Phaeocystis TPL expression in a sample)
    - dO2_resid_dt_std: The first derivative of residual O2 values (see carbonate_chemistry.R), z-standardized
    - Date: Date of sampling
    
    The function returns a dictionary with the following keys:
    - query_id: Transcript ID
    - O2_coef: Coefficient for dO2_resid_dt_std
    - O2_pval: P-value for dO2_resid_dt_std
    - AIC: AIC value for the model
    """
    try:
        query_id = transcript_data['query_id'].iloc[0]
        y = transcript_data['TPL_standardized']
        transcript_data['Time'] = (transcript_data['Date'] - transcript_data['Date'].min()).dt.total_seconds() / 3600
        splines = dmatrix("bs(Time, df=5, degree=3)", data=transcript_data, return_type='dataframe')

        # Combine derivative predictor with splines
        X = pd.concat([transcript_data[['dO2_resid_dt_std']], splines], axis=1)
        X = sm.add_constant(X)

        if X.isnull().values.any() or np.isinf(X.values).any():
            print(f"NaNs or infs in predictors for transcript {query_id}")
            return None

        model = GLSAR(y, X, rho=1)
        result = model.iterative_fit()

        return {
            "query_id": query_id,
            "dO2_resid_dt_coef": result.params.get("dO2_resid_dt_std", np.nan),
            "dO2_resid_dt_pval": result.pvalues.get("dO2_resid_dt_std", np.nan),
            "AIC": result.aic
        }
    except Exception as e:
        print(f"Model failed for transcript {query_id}: {e}")
        return None

# Main script
if __name__ == "__main__":
    try:
        # Load data
        logging.info("Loading data...")
        O2 = pd.read_csv('data/analysis/O2_model_station130_autocorr.csv')
        tpl = pd.read_csv('data/quantification/130/130_tpl.csv').melt(id_vars='target_id', var_name='sample', value_name='tpl')
        phaeo = pd.read_csv('data/annotation/taxonomy_eukprot/130/genus_bins/Phaeocystis_transcriptome_bin.csv').drop(columns='target_id')
        meta = pd.read_csv('data/samples_env.csv')

        # Preprocess meta and O2
        meta['Date'] = pd.to_datetime(meta['Date'], utc=True)
        meta = meta.loc[meta['StationPrefix'] == 130]
        O2['Date'] = pd.to_datetime(O2['Date'], format='mixed', utc=True).dt.floor('s')

        # Interpolate dO2_resid_dt values to match sample dates
        O2_augmented = pd.concat([O2, pd.DataFrame({'Date': meta['Date']})]).sort_values(by='Date')
        O2_augmented = O2_augmented.set_index('Date').interpolate(method='time').reset_index()

        meta = pd.merge(meta, O2_augmented[['Date', 'dO2_resid_dt']], on='Date', how='left')
        meta['Time'] = (meta['Date'] - meta['Date'].min()).dt.total_seconds() / 3600

        # Merge transcript expression with taxonomic bin
        data = pd.merge(tpl, phaeo, left_on='target_id', right_on='query_id', how='right')
        data['TPL_standardized'] = data['tpl'] / data.groupby('sample')['tpl'].transform('sum')

        model_data = pd.merge(data, meta, left_on='sample', right_on='Station', how='inner')

        # Standardize derivative predictor
        model_data['dO2_resid_dt_std'] = (
            model_data['dO2_resid_dt'] - model_data['dO2_resid_dt'].mean()
        ) / model_data['dO2_resid_dt'].std()

        # Filter transcripts: at least 50% of samples have non-zero expression
        logging.info("Filtering transcripts...")
        transcript_stats = model_data.groupby('query_id')['TPL_standardized'].agg(
            nonzero_count=lambda x: (x > 0).sum(),
            total_count='count',
        )
        transcript_stats['nonzero_fraction'] = transcript_stats['nonzero_count'] / transcript_stats['total_count']
        filtered_transcripts = transcript_stats[transcript_stats['nonzero_fraction'] >= 0.5].index
        logging.info(f"Number of transcripts retained after filtering: {len(filtered_transcripts)}")

        model_data_filtered = model_data[model_data['query_id'].isin(filtered_transcripts)]

        # Group by transcript for parallel processing
        grouped = [group for _, group in model_data_filtered.groupby('query_id')]
        logging.info(f"Number of groups to process: {len(grouped)}")

        # Run GLSAR in parallel
        logging.info("Fitting GLSAR models...")
        with Pool(cpu_count()) as pool:
            results = pool.map(fit_glsar_with_splines, grouped)

        # Collect and save results
        results = [res for res in results if res is not None]
        results_df = pd.DataFrame(results)

        if not results_df.empty and 'dO2_resid_dt_pval' in results_df:
            logging.info("Adjusting p-values for multiple testing...")
            results_df['dO2_resid_dt_pval_adj'] = multipletests(results_df['dO2_resid_dt_pval'], method='fdr_bh')[1]

            # Save significant markers
            significant_markers = results_df[results_df['dO2_resid_dt_pval_adj'] < 0.05]
            logging.info(f"Number of significant markers: {len(significant_markers)}")
            significant_markers.to_csv('data/analysis/significant_derivative_markers.csv', index=False)
        else:
            logging.warning("No valid results or p-values found.")
    except Exception as e:
        logging.error(f"Script failed: {e}")
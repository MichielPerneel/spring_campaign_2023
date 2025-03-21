import pandas as pd
from statsmodels.regression.linear_model import GLSAR
from patsy import dmatrix
from statsmodels.stats.multitest import multipletests
import numpy as np
from multiprocessing import Pool, cpu_count
import logging

# Set up logging
logging.basicConfig(
    filename="data/logs/get_marker_transcripts.log",
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
    - O2_std: Standardized O2 values (O2 values standardized by mean and standard deviation)
    - Date: Date of sampling
    
    The function returns a dictionary with the following keys:
    - query_id: Transcript ID
    - O2_coef: Coefficient for O2_std
    - O2_pval: P-value for O2_std
    - AIC: AIC value for the model
    """
    try:
        query_id = transcript_data['query_id'].iloc[0]

        # Debugging: Print basic statistics about the data
        print(f"Processing transcript {query_id}")
        print(f"TPL_standardized stats:\n{transcript_data['TPL_standardized'].describe()}")
        print(f"O2_std stats:\n{transcript_data['O2_std'].describe()}")

        # Prepare predictors and response
        y = transcript_data['TPL_standardized']
        transcript_data['Time'] = (transcript_data['Date'] - transcript_data['Date'].min()).dt.total_seconds() / 3600
        splines = dmatrix("bs(Time, df=5, degree=3)", data=transcript_data, return_type='dataframe')
        
        # Check for NaN or inf in splines
        if splines.isnull().values.any():
            print(f"Splines contain NaNs for transcript {query_id}")
            return None
        
        X = pd.concat([transcript_data[['O2_std']], splines], axis=1)
        X = sm.add_constant(X)  # Add intercept
        
        # Check for NaN or inf in predictors
        if X.isnull().values.any() or np.isinf(X.values).any():
            print(f"Predictors contain NaNs or infs for transcript {query_id}")
            return None

        # Fit GLSAR model with AR1 structure
        model = GLSAR(y, X, rho=1)
        result = model.iterative_fit()

        # Extract results
        return {
            "query_id": query_id,
            "O2_coef": result.params.get("O2_std", np.nan),
            "O2_pval": result.pvalues.get("O2_std", np.nan),
            "AIC": result.aic
        }
    except Exception as e:
        print(f"Model failed for transcript {query_id}: {e}")
        return None

# Main script
if __name__ == "__main__":
    import statsmodels.api as sm

    try:
        # Load data
        logging.info("Loading data...")
        O2 = pd.read_csv('data/analysis/O2_smoother_station_130.csv')
        tpl = pd.read_csv('data/quantification/130/130_tpl.csv').melt(id_vars='target_id', var_name='sample', value_name='tpl')
        phaeo = pd.read_csv('data/annotation/taxonomy_eukprot/130/genus_bins/Phaeocystis_transcriptome_bin.csv').drop(columns='target_id')
        meta = pd.read_csv('data/samples_env.csv')

        # Preprocessing
        logging.info("Preprocessing data...")
        meta['Date'] = pd.to_datetime(meta['Date'])
        meta = meta.loc[meta['StationPrefix'] == 130]
        O2['Date'] = pd.to_datetime(O2['Date'], format='mixed').dt.floor('s')
        meta['Time'] = (meta['Date'] - meta['Date'].min()).dt.total_seconds() / 3600  # Hours since start

        ## Interpolate O2 values for meta timestamps
        ### Step 1: Add meta['Date'] timestamps to O2
        O2_augmented = pd.concat([O2, pd.DataFrame({'Date': meta['Date']})]).sort_values(by='Date')
        ### Step 2: Interpolate O2 at the added timestamps
        O2_augmented = O2_augmented.set_index('Date').interpolate(method='time').reset_index()

        ### Step 3: Extract interpolated values at meta['Date'] and merge into meta
        meta = meta.sort_values(by='Date')
        meta = pd.merge(meta, O2_augmented, on='Date', how='left')
        meta.rename(columns={'y': 'O2'}, inplace=True)

        # Merge transcript data
        data = pd.merge(tpl, phaeo, left_on='target_id', right_on='query_id', how='right')
        data['TPL_standardized'] = data['tpl'] / data.groupby('sample')['tpl'].transform('sum')
        model_data = pd.merge(data, meta, left_on='sample', right_on='Station', how='inner')

        # Standardize preO2tors
        model_data['O2_std'] = (model_data['O2'] - model_data['O2'].mean()) / model_data['O2'].std()

        # Filter transcripts expressed in at least 50% of samples
        logging.info("Filtering transcripts...")
        transcript_stats = model_data.groupby('query_id')['TPL_standardized'].agg(
            nonzero_count=lambda x: (x > 0).sum(),
            total_count='count',
        )
        transcript_stats['nonzero_fraction'] = transcript_stats['nonzero_count'] / transcript_stats['total_count']
        filtered_transcripts = transcript_stats[transcript_stats['nonzero_fraction'] >= 0.5].index
        logging.info(f"Number of transcripts retained after filtering: {len(filtered_transcripts)}")
        model_data_filtered = model_data[model_data['query_id'].isin(filtered_transcripts)]

        # Group data by transcript
        grouped = [group for _, group in model_data_filtered.groupby('query_id')]
        logging.info(f"Number of groups to process: {len(grouped)}")

        # Parallel processing for GLSAR models
        logging.info("Fitting GLSAR models...")
        with Pool(cpu_count()) as pool:
            results = pool.map(fit_glsar_with_splines, grouped)

        # Collect results
        results = [res for res in results if res is not None]
        logging.info(f"Number of successful models: {len(results)}")
        results_df = pd.DataFrame(results)

        # Multiple testing correction
        if 'O2_pval' in results_df:
            logging.info("Adjusting p-values for multiple testing...")
            results_df['O2_pval_adj'] = multipletests(results_df['O2_pval'], method='fdr_bh')[1]

            # Filter significant markers
            significant_markers = results_df[results_df['O2_pval_adj'] < 0.05]
            logging.info(f"Number of significant markers: {len(significant_markers)}")
            significant_markers.to_csv('data/analysis/significant_markers.csv', index=False)
        else:
            logging.warning("No valid p-values found in results.")

    except Exception as e:
        logging.error(f"Script failed: {e}")
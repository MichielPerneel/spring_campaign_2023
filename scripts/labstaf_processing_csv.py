import os
import pandas as pd

def parse_csv_folder(folder_path):
    """
    Parses all CSV files in a folder, extracting the relevant values under the 'rPE fit' section,
    and creates a concatenated dataframe with samples as rows and extracted values as columns.
    """
    extracted_data = []
    target_columns = ["Alpha", "Beta", "Ek", "EkBeta", "rPm", "JVPIIm", "GOPIIm", "Fo", "Fm", "Fv", "Fv/Fm", 
                      "Fv/Fmc", "F'", "Fm'", "Fq'", "Fq'/Fm'", "Fq'/Fmc'", "NPQ", "NSV", "SigmaPII_rhofit"]
    
    for file_name in sorted(os.listdir(folder_path)):
        if file_name.endswith(".csv"):
            file_path = os.path.join(folder_path, file_name)
            try:
                # Read file as CSV, keeping only first two columns
                df = pd.read_csv(file_path, delimiter=";", usecols=[0, 1], encoding="latin1", names=["Parameter", "Value"], skiprows=1, header=None)
                # Only first 20 rows are relevant
                df = df.iloc[:20]
                df = df.dropna().set_index("Parameter").T  # Transpose to have parameters as columns
                
                # Remove : from parameter names, and strip whitespace and tab characters
                df.columns = df.columns.str.replace(":", "").str.strip().str.replace("\t", "")
                
                # Ensure all target columns exist
                for col in target_columns:
                    if col not in df.columns:
                        df[col] = None
                
                df = df[target_columns]  # Keep columns in the correct order
                sample_name = os.path.splitext(file_name)[0]  # Remove .csv extension
                df.index = [sample_name]
                extracted_data.append(df)
            except Exception as e:
                print(f"Error processing {file_name}: {e}")
    # Concatenate all DataFrames
    final_df = pd.concat(extracted_data)
    
    # Convert numeric columns
    final_df[final_df.columns[1:]] = final_df[final_df.columns[1:]].apply(pd.to_numeric, errors='coerce')
    
    # Apply zero-floor correction to parameters that should not be negative
    zero_floor_columns = ["NPQ", "Fo", "Fm", "Fv", "Fv/Fm", "Fv/Fmc", "F'", "Fm'", "Fq'", "Fq'/Fm'", "Fq'/Fmc'"]
    final_df[zero_floor_columns] = final_df[zero_floor_columns].clip(lower=0)
    
    return final_df

station=130
folder_path = "data/raw/LabSTAF/csv_simple_files/{}".format(station)
df = parse_csv_folder(folder_path)

# Estimate Primary Production by multiplying the GOPIIm by 12
df["PP"] = df["GOPIIm"] * 12

# Save to CSV
df.to_csv("data/raw/LabSTAF/labstaf_csv_data_{}.csv".format(station))
print(df.head())
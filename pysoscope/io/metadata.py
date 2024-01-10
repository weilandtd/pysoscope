import pandas as pd 

def load_metadata_csv(file_path, sample_name):
    """AI is creating summary for read_metadata_csv

    Args:
        file_path ([type]): [description]
        sample_name ([type]): [description]

    Raises:
        ValueError: [description]

    Returns:
        [type]: [description]
    """
    metadata_df = pd.read_csv(file_path, header=0)
    
    if 'sample_name' not in metadata_df.columns:
        raise ValueError("The 'sample_name' column does not exist in the DataFrame.")
    
    filtered_df = metadata_df[metadata_df['sample_name'].isin(sample_name)]
    
    # Convert filtered_df into a list of dictionaries
    metadata_list = filtered_df.to_dict('records')
    
    return metadata_list

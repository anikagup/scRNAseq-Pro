import json
from preprocessing import preprocess_data, load_data

# Load configuration
with open("config.json", "r") as config_file:
    config = json.load(config_file)

input_file = config.get("input_file", "data/input_file")
file_type = config.get("file_type", "auto")
params = config.get("preprocessing_params", {})

adata = load_data(input_file, file_type)
if adata is not None:
    adata = preprocess_data(adata, params)
    adata.write("data/processed_data.h5ad")
else:
    print("No valid data loaded.")

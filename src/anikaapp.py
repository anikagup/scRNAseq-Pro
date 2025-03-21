from shiny import App, ui, render, reactive
import os
import shutil
import subprocess
import json
import scanpy as sc

# Print statement to confirm app start
print("App is starting...")

# Load configuration
# Get the absolute path to the root of the project
project_root = os.path.abspath(os.path.join(os.getcwd(), '..', '..'))
print(project_root)
config_path = os.path.join(project_root, 'scRNA-seq-Automation', 'src', 'config.json')
print(config_path)
with open(config_path, 'r') as f:
    config = json.load(f)

# Define the UI layout
app_ui = ui.page_fluid(
    ui.tags.style("body { background-color: lightblue; }"),
    ui.panel_title("Welcome to RNA Pro"),

    # File upload for multiple formats
    ui.input_file(
        "file_upload",
        "Upload scRNA-seq Data",
        multiple=False,
        accept=[".fastq", ".csv", ".h5ad", ".h5", ".loom"]
    ),

    ui.output_text("file_info"),
    ui.input_action_button("activate_button_ui", "Run analysis"),  

    # QC Parameter Inputs
    ui.h3("Modify QC Metrics"),
    ui.input_numeric("min_genes", "Min Genes per Cell", value=config["preprocessing_params"]["min_genes"], min=50, max=1000, step=50),
    ui.input_numeric("min_cells", "Min Cells per Gene", value=config["preprocessing_params"]["min_cells"], min=1, max=50, step=1),
    ui.input_numeric("target_sum", "Normalization Target Sum", value=config["preprocessing_params"]["target_sum"], min=1000, max=50000, step=1000),
    ui.input_action_button("reprocess_button", "Recalculate QC and Preprocessing"),
    
    ui.output_text("reprocess_status"),

    # UMAP Gene Selection
    ui.h3("Visualize UMAP"),
    ui.input_text("gene_input", "Enter Genes (comma-separated):", placeholder="E.g., CST3, NKG7"),
    ui.input_action_button("update_umap", "Generate UMAP"),
    ui.output_image("umap_plot"),
)

def server(input, output, session):
    upload_dir = os.path.join(os.getcwd(), 'uploads')
    os.makedirs(upload_dir, exist_ok=True)
    adata = None  # Placeholder for loaded data

    # Function to save uploaded file
    def save_uploaded_file():
        file = input.file_upload()
        if file:
            temp_path = file[0]["datapath"]
            saved_path = os.path.join(upload_dir, file[0]["name"])
            shutil.move(temp_path, saved_path)
            return saved_path
        return None

    @output
    @render.text
    def file_info():
        file = input.file_upload()
        if file is None:
            return "No file uploaded."
        
        saved_path = save_uploaded_file()
        if saved_path:
            # Update config.json
            config["input_file"] = saved_path
            with open(config_path, "w") as f:
                json.dump(config, f, indent=4)

            return f"File Saved: {saved_path}"
        return "Error saving file."
    # Function to re-run preprocessing with new QC thresholds

    @reactive.effect
    @reactive.event(input.activate_button_ui)
    def activate_analysis():
        # Run the main.py script in the 'src' directory
        script_path = os.path.join(project_root, 'scRNA-seq-Automation', 'src', 'main.py')
        try:
            subprocess.run(['python', script_path], check=True)  # Run the script
            print("main.py has been executed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error running main.py: {e}")        
 

    @output
    @render.text
    def reprocess_status():
        return "Press 'Recalculate QC' to apply new metrics."

    # Function to re-run preprocessing with new QC thresholds
    @reactive.effect
    @reactive.event(input.reprocess_button)
    def reprocess_data():
        config["preprocessing_params"]["min_genes"] = input.min_genes()
        config["preprocessing_params"]["min_cells"] = input.min_cells()
        config["preprocessing_params"]["target_sum"] = input.target_sum()

        with open(config_path, "w") as f:
            json.dump(config, f, indent=4)

        subprocess.run(["python", "src/main.py"])
        session.send_notification("info", "Reprocessing complete!")

    # Function to generate UMAP with user-defined genes
    @output
    @render.image
    @reactive.event(input.update_umap)
    def umap_plot():
        gene_list = input.gene_input().split(",")
        gene_list = [g.strip() for g in gene_list if g.strip()]

        if not gene_list:
            return None  # No valid genes entered

        # Load updated dataset
        global adata
        adata = sc.read_h5ad("data/processed_data.h5ad")

        # Check if genes exist
        valid_genes = [gene for gene in gene_list if gene in adata.var_names]
        if not valid_genes:
            return None  # No valid genes found

        # Generate UMAP
        save_path = "figures/_custom_umap.png"
        sc.pl.umap(adata, color=valid_genes, save="_custom_umap.png")
        return {"src": save_path, "height": "500px"}

app = App(app_ui, server)
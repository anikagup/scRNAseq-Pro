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
project_root = os.path.abspath(os.path.join(os.getcwd(), '..'))
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
    ui.h3("Violin QC"),
    ui.output_image("displayed_image1"),
    ui.h3("UMAP QC"),
    ui.output_image("displayed_image2"),   
    ui.h3("UMAP Top 5"),
    ui.output_image("displayed_image3"), 
    ui.h3("Ranked Genes"),
    ui.output_image("displayed_image4"), 
    ui.h3("ML UMAP"),
    ui.output_image("displayed_image5"), 

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
    """ def save_uploaded_file():
        file = input.file_upload()
        if file:
            base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))  # Move up one level
            upload_dir = os.path.join(base_dir, "UI/uploads")

            # Save the uploaded file to the uploads folder
            temp_path = file[0]["datapath"]
            saved_path = os.path.join(upload_dir, file[0]["name"])
            shutil.move(temp_path, saved_path)
            return saved_path
        return None """

    # Supported file extensions
    ALLOWED_EXTENSIONS = {
        ".h5": "10x",
        ".loom": "loom",
        ".h5ad": "h5ad",
        ".csv": "csv",
        ".txt": "txt"
    }

    def save_uploaded_file():
        file = input.file_upload()
        if not file:
            return None

        # Project root and paths
        base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
        upload_dir = os.path.join(base_dir, "uploads")  # <== NOT UI/uploads anymore
        config_path = os.path.join(base_dir, "src", "config.json")

        # Get extension
        uploaded_name = file[0]["name"]
        uploaded_ext = os.path.splitext(uploaded_name)[1]

        if uploaded_ext not in ALLOWED_EXTENSIONS:
            print(f"❌ Unsupported file type: {uploaded_ext}")
            return None

        # Save with a standard name (e.g., uploads/latest.h5ad)
        temp_path = file[0]["datapath"]
        save_path = os.path.join(upload_dir, f"latest{uploaded_ext}")
        shutil.move(temp_path, save_path)

        # Update config.json
        try:
            with open(config_path, "r") as f:
                config = json.load(f)
        except FileNotFoundError:
            config = {}

        config["input_file"] = f"uploads/latest{uploaded_ext}"
        config["file_type"] = ALLOWED_EXTENSIONS[uploaded_ext]

        with open(config_path, "w") as f:
            json.dump(config, f, indent=4)

        print(f"✅ File saved to: {save_path}")
        print(f"✅ config.json updated.")

        return save_path

    @output
    @render.text
    def file_info():
        saved_path = save_uploaded_file()
        if saved_path:
            return f"✅ File saved and config updated: {saved_path}"
        return "❌ File upload failed or unsupported format."

    # Function to clear uploads folder when app is run
    def clear_uploads_folder():
        upload_dir = os.path.join(os.getcwd(), 'uploads')
        if os.path.exists(upload_dir):
            shutil.rmtree(upload_dir)  # Delete the folder and its contents
        os.makedirs(upload_dir, exist_ok=True)  # Recreate the folder

    # Clear the uploads folder when the app is run
    clear_uploads_folder()

    # Function to clear uploads folder when app is run
    def clear_processed_data_folder():
        base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
        processed_data_dir = os.path.join(base_dir, "processed_data")
        if os.path.exists(processed_data_dir):
            shutil.rmtree(processed_data_dir)  # Delete the folder and its contents
        os.makedirs(processed_data_dir, exist_ok=True)  # Recreate the folder

    # Clear the uploads folder when the app is run
    clear_processed_data_folder()

    @reactive.effect
    @reactive.event(input.activate_button_ui)
    def activate_analysis():
        # Clearing CSV files and figures every time analysis is pressed
        datapath = os.path.join(project_root, 'scRNA-seq-Automation', 'data')
        if os.path.exists(datapath):
            shutil.rmtree(datapath) 
        # Clearing figures every time analysis is pressed
        figurepath = os.path.join(project_root, 'scRNA-seq-Automation', 'figures')
        if os.path.exists(figurepath):
            shutil.rmtree(figurepath) 
        
        # Run the main.py script in the 'src' directory
        script_path = os.path.join(project_root, 'scRNA-seq-Automation', 'src', 'main.py')
        try:
            subprocess.run(['python', script_path], check=True)  # Run the script
            print("main.py has been executed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error running main.py: {e}")    

    @output
    @render.image
    @reactive.event(input.activate_button_ui)
    def displayed_image1():
        image_path = os.path.join(project_root, 'scRNA-seq-Automation', 'figures', 'violin_qc_metrics.png')
        if os.path.exists(image_path):
            return {"src": image_path, "height": "400px"}  # Return image with height setting
        return None  # Return None if image is not found

    @output
    @render.image
    @reactive.event(input.activate_button_ui)
    def displayed_image2():
        image_path = os.path.join(project_root, 'scRNA-seq-Automation', 'figures', 'umap_qc.png')
        if os.path.exists(image_path):
            return {"src": image_path, "height": "250px"}  # Return image with height setting
        return None  # Return None if image is not found

    @output
    @render.image
    @reactive.event(input.activate_button_ui)
    def displayed_image3():
        image_path = os.path.join(project_root, 'scRNA-seq-Automation', 'figures', 'umap_top5.png')
        if os.path.exists(image_path):
            return {"src": image_path, "height": "400px"}  # Return image with height setting
        return None  # Return None if image is not found
    
    """@output
    @render.image
    @reactive.event(input.activate_button_ui)
    def displayed_image3a():
        image_path = os.path.join(project_root, 'scRNA-seq-Automation', 'figures', 'umap_custom_gene.png')
        print("🔍 Looking for:", image_path)
        print("✅ Exists?", os.path.exists(image_path))

        if os.path.exists(image_path):
            return {"src": image_path, "height": "400px"}
        return None """

    @output
    @render.image
    @reactive.event(input.activate_button_ui)
    def displayed_image4():
        image_path = os.path.join(project_root, 'scRNA-seq-Automation', 'figures', 'rank_genes_groups_leiden.png')
        if os.path.exists(image_path):
            return {"src": image_path, "height": "400px"}  # Return image with height setting
        return None  # Return None if image is not found
    
    @output
    @render.image
    @reactive.event(input.activate_button_ui)
    def displayed_image5():
        image_path = os.path.join(project_root, 'scRNA-seq-Automation', 'figures', 'ML_UMAPgit .png')
        if os.path.exists(image_path):
            return {"src": image_path, "height": "400px"}  # Return image with height setting
        return None  # Return None if image is not found

    @output
    @render.text
    def reprocess_status():
        return "Press 'Recalculate QC' to apply new metrics."

    @reactive.effect
    @reactive.event(input.reprocess_button)
    def reprocess_data():
        # Path to the figures folder
        figurepath = os.path.join(project_root, 'scRNA-seq-Automation', 'figures')

        # Clear the figures folder by deleting all its contents
        if os.path.exists(figurepath):
            shutil.rmtree(figurepath)  # Delete the folder and its contents
        os.makedirs(figurepath, exist_ok=True)  # Recreate the folder

        # Clear the displayed images on the website by setting them to None
        @output
        @render.image
        def displayed_image1():
            return None

        @output
        @render.image
        def displayed_image2():
            return None

        @output
        @render.image
        def displayed_image3():
            return None

        @output
        @render.image
        def displayed_image4():
            return None
        
        @output
        @render.image
        def displayed_image5():
            return None

        # Save the original input_file and file_type before updating other parts of the config
        original_input_file = config.get("input_file")
        original_file_type = config.get("file_type")

        # Update only the QC parameters
        config["preprocessing_params"]["min_genes"] = input.min_genes()
        config["preprocessing_params"]["min_cells"] = input.min_cells()
        config["preprocessing_params"]["target_sum"] = input.target_sum()

        # Ensure input_file and file_type remain unchanged
        if original_input_file is not None:
            config["input_file"] = original_input_file
        if original_file_type is not None:
            config["file_type"] = original_file_type

        # Write updated config back to file without overwriting input_file or file_type
        with open(config_path, "w") as f:
            json.dump(config, f, indent=4)
            print(config)

        # Re-run the analysis after the preprocessing update
        subprocess.run(["python", "src/main.py"])

        # Force the UI to update and display the new images by triggering re-render
        @output
        @render.image
        def displayed_image1():
            image_path = os.path.join(project_root, 'scRNA-seq-Automation', 'figures', 'violin_qc_metrics.png')
            if os.path.exists(image_path):
                return {"src": image_path, "height": "400px"}
            return None

        @output
        @render.image
        def displayed_image2():
            image_path = os.path.join(project_root, 'scRNA-seq-Automation', 'figures', 'umap_qc.png')
            if os.path.exists(image_path):
                return {"src": image_path, "height": "250px"}
            return None

        @output
        @render.image
        def displayed_image3():
            image_path = os.path.join(project_root, 'scRNA-seq-Automation', 'figures', 'umap_top5.png')
            if os.path.exists(image_path):
                return {"src": image_path, "height": "400px"}
            return None

        @output
        @render.image
        def displayed_image4():
            image_path = os.path.join(project_root, 'scRNA-seq-Automation', 'figures', 'rank_genes_groups_leiden.png')
            if os.path.exists(image_path):
                return {"src": image_path, "height": "400px"}
            return None
        
        @output
        @render.image
        def displayed_image5():
            image_path = os.path.join(project_root, 'scRNA-seq-Automation', 'figures', 'ML_UMAP.png')
            if os.path.exists(image_path):
                return {"src": image_path, "height": "400px"}
            return None

        # Update the UI to show the completion message by updating the 'reprocess_status' output
        @output
        @render.text
        def reprocess_status():
            return "Reprocessing complete!"

# Create the app
app = App(app_ui, server)



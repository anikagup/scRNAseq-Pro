from shiny import App, ui, render, reactive
import os
import shutil
import subprocess
import json
import scanpy as sc

print("App is starting...")

project_root = os.path.dirname(__file__)
print(project_root)
config_path = os.path.join(project_root, "../src/config.json")
print(config_path)
with open(config_path, 'r') as f:
    config = json.load(f)

def clear_folder(folder_path):
    if os.path.exists(folder_path):
        for filename in os.listdir(folder_path):
            file_path = os.path.join(folder_path, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print(f"⚠️ Failed to delete {file_path}. Reason: {e}")

def styled_button(input_id, label):
    return ui.input_action_button(input_id, label, style="background-color: #007bff; color: white; border: none; padding: 10px 20px; border-radius: 5px;")

app_ui = ui.page_fluid(
    ui.tags.style("body { background-color: white; }"),
    ui.panel_title(ui.tags.h1("Welcome to RNA Pro", style="font-weight: bold; font-size: 2em;")),
    ui.tags.br(),
    ui.tags.br(),

    ui.input_file("file_upload", "Upload scRNA-seq Data", multiple=False, accept=[".fastq", ".csv", ".h5ad", ".h5", ".loom"]),
    ui.output_text("file_info"),
    ui.tags.br(),

    styled_button("activate_button_ui", "Run analysis"),
    ui.tags.br(),
    ui.tags.br(),

    ui.h3("Modify QC Metrics", style="text-decoration: underline;"),
    ui.tags.br(),
    ui.h4("Min Genes per Cell"),
    ui.div(
        ui.tags.span("📝 This sets the minimum number of genes that must be detected in each cell for it to be included in the analysis. It helps filter out low-quality or dying cells that have little RNA content. A common threshold is 200-500 genes per cell but this may vary depending on tissue type and sequencing depth.", style="color: white;"),
        style="background-color: #28a745; padding: 10px; border-radius: 5px; box-shadow: 0px 0px 5px rgba(0, 0, 0, 0.1); width: 40%;"
    ),
    ui.input_numeric("min_genes", " ", value=config["preprocessing_params"]["min_genes"], min=50, max=1000, step=50),
    ui.h4("Min Cells per Gene"),
    ui.div(
        ui.tags.span("📝 This sets the minimum number of cells a gene must appear in to be kept in the dataset. It removes genes that are only expressed in a few cells and may represent noise. A typical threshold is 3–10 cells per gene, which balances removing artifact noise without discarding biologically relevant genes.", style="color: white;"),
        style="background-color: #28a745; padding: 10px; border-radius: 5px; box-shadow: 0px 0px 5px rgba(0, 0, 0, 0.1); width: 40%;"
    ),
    ui.input_numeric("min_cells", " ", value=config["preprocessing_params"]["min_cells"], min=1, max=50, step=1),
    ui.h4("Normalization Target Sum"),
    ui.div(
        ui.tags.span("📝 This sets the total gene expression count to which each cell is scaled, allowing fair comparison between cells with different sequencing depths. A typical value is 10,000.", style="color: white;"),
        style="background-color: #28a745; padding: 10px; border-radius: 5px; box-shadow: 0px 0px 5px rgba(0, 0, 0, 0.1); width: 40%;"
    ),
    ui.input_numeric("target_sum", " ", value=config["preprocessing_params"]["target_sum"], min=1000, max=50000, step=1000),
    styled_button("reprocess_button", "Recalculate QC and Preprocessing"),
    ui.output_text("reprocess_status"),

    ui.tags.br(),
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

    ui.tags.br(),
    styled_button("downloadData1", "Download Processed CSV"),
    ui.output_text("fileStatus_1"),

    ui.tags.br(),
    styled_button("download_deg", "Download Differential Gene Expression List"),
    ui.output_text("fileStatus_2"),

    ui.tags.br(),
    ui.h3("Visualize Custom Gene-Specific UMAP"),
    ui.input_text("gene_input", "Enter Gene:", placeholder="E.g., 'CST3' or 'NKG7' "),
    styled_button("update_umap", "Generate UMAP"),
    ui.output_text("gene_status"),
    ui.tags.br(),
    ui.output_image("displayed_image6"),
)

uploaded_path = reactive.Value(None)
def resolve_deg_path():
    # Absolute path to the DEG file in the container
    return os.path.abspath("/app/processed_data/all_degs.csv")

def server(input, output, session):
    file_ready = reactive.Value(False)

    ALLOWED_EXTENSIONS = {".h5": "10x",".loom": "loom", ".h5ad": "h5ad", ".csv": "csv", ".txt": "txt"}

    @reactive.effect
    @reactive.event(input.file_upload)
    def handle_upload():
        clear_folder(os.path.join(os.path.dirname(__file__), "../uploads"))


        # Function to clear uploads folder when app is run
        clear_folder(os.path.join(os.path.dirname(__file__), "../processed_data"))

        


        file = input.file_upload()
        if not file:
            return

        base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
        upload_dir = os.path.join(base_dir, "uploads")
        config_path = os.path.join(base_dir, "src", "config.json")
        uploaded_name = file[0]["name"]
        uploaded_ext = os.path.splitext(uploaded_name)[1]
        if uploaded_ext not in ALLOWED_EXTENSIONS:
            print(f"❌ Unsupported file type: {uploaded_ext}")
            return None
        
        temp_path = file[0]["datapath"]
        save_path = os.path.join(upload_dir, f"latest{uploaded_ext}")
        shutil.move(temp_path, save_path)

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

        uploaded_path.set(save_path)

        return save_path


    @output
    @render.text
    def file_info():
        path = uploaded_path()
        return f"✅ File saved and config updated: {path}" if path else "❌ Please upload an accepted file type."
    
    @reactive.effect
    @reactive.event(input.activate_button_ui)
    def activate_analysis():
        subprocess.run(["python", "/app/src/main.py"], check=True)
        print("main.py has been executed successfully")
        processed_csv_path = "/app/processed_data/processed_matrix.csv"
        deg_path = "/app/processed_data/all_degs.csv"
        if os.path.exists(processed_csv_path):
            print(f"✅ Found Cell Count Matrix at: {processed_csv_path}")
            file_ready.set(True)
        else:
            print(f"❌ Cell Count Matrix not found at {processed_csv_path}")
            file_ready.set(False)
        file_ready.set(True)
        if os.path.exists(deg_path):
            print(f"✅ Found DEG List at: {deg_path}")
            file_ready.set(True)
        else:
            print(f"❌ DEG list not found at {deg_path}")
            file_ready.set(False)
        file_ready.set(True)

    processed_csv_path = os.path.join("/app", "processed_data", "processed_matrix.csv")
    deg_path = os.path.join("/app", "/processed_data", "all_degs.csv")

    # Always render download button
    @output
    @render.download
    def downloadData1():
        if os.path.exists(processed_csv_path):
            return processed_csv_path
        else:
            # If file doesn't exist, prevent download by returning None
            return None
    
       
    @output
    @render.download
    def download_deg():
        path = resolve_deg_path()
        if os.path.exists(path):
            return path
        return None



    @output
    @render.image
    @reactive.event(input.activate_button_ui, input.reprocess_button)
    def displayed_image1():
        image_path = os.path.join(project_root, '..', 'figures', 'violin_qc_metrics.png')
        return {"src": image_path, "height": "350px"} if file_ready() and os.path.exists(image_path) else None
    @output
    @render.ui
    @reactive.event(input.activate_button_ui)
    def violin_qc_text():
        image_path = os.path.join(project_root, '..', 'figures', 'violin_qc_metrics.png')
        if os.path.exists(image_path):
            return ui.HTML(
                """
                <div style="margin-top: 15px; margin-bottom: 30px; padding: 15px; max-width: 1000px; background-color: #f8f9fa; border-left: 5px solid #007bff; border-radius: 5px; font-size: 15px; box-shadow: 2px 2px 8px rgba(0,0,0,0.05);">
                     📝 These violin plots display the distribution of three key scRNA-seq quality control metrics across all cells.<br><br>
                    
                    1. <strong>n_genes_by_counts</strong>: This plot shows how many genes are detected in each cell. Cells appearing at the bottom of the violin (low vertical values) likely have very few genes detected and may represent empty droplets or low-quality cells. A common threshold is filtering out cells below ~200 genes.<br><br>
                    
                    2. <strong>total_counts</strong>: This indicates the total number of RNA molecules sequenced per cell. Cells with low total counts (bottom of the violin) may be poor quality, while those with very high counts (top tail) could be doublets.<br><br>
                    
                    3. <strong>pct_counts_mt</strong>: This shows the percentage of mitochondrial RNA per cell. A long upper tail in this plot may suggest a population of stressed or dying cells. It's common practice to filter out cells with more than 5–10% mitochondrial content.
                </div>
                """
            )
        return None

    @output
    @render.image
    @reactive.event(input.activate_button_ui, input.reprocess_button)
    def displayed_image2():
        image_path = os.path.join(project_root, '..', 'figures', 'umap_qc.png')
        return {"src": image_path, "height": "250px"} if file_ready() and os.path.exists(image_path) else None

    @output
    @render.image
    @reactive.event(input.activate_button_ui, input.reprocess_button)
    def displayed_image3():
        image_path = os.path.join(project_root, '..', 'figures', 'umap_top5.png')
        return {"src": image_path, "height": "400px"} if file_ready() and os.path.exists(image_path) else None

    @output
    @render.image
    @reactive.event(input.activate_button_ui, input.reprocess_button)
    def displayed_image4():
        image_path = os.path.join(project_root, '..', 'figures', 'rank_genes_groups_leiden.png')
        return {"src": image_path, "height": "400px"} if file_ready() and os.path.exists(image_path) else None

    @output
    @render.image
    @reactive.event(input.activate_button_ui, input.reprocess_button)
    def displayed_image5():
        image_path = os.path.join(project_root, '..', 'figures', 'umapML_umap.png')
        return {"src": image_path, "height": "400px"} if file_ready() and os.path.exists(image_path) else None

    @output
    @render.image
    @reactive.event(input.update_umap)
    def displayed_image6():
        with open(config_path, "r") as f:
            config = json.load(f)
        config["visualization"]["custom_genes"] = input.gene_input()
        with open(config_path, "w") as f:
            json.dump(config, f, indent=4)
        upload_path = os.path.join("/app", 'figures', 'umap_custom_gene.png')
        print(upload_path)

        if os.path.exists(upload_path):
            os.remove(upload_path)
        
        subprocess.run(["python", "/app/src/main.py"], check=True)
        

        @output
        @render.text
        def gene_status():
            return None

        image_path = os.path.join(project_root, '..', 'figures', 'umap_custom_gene.png')
        if os.path.exists(image_path):
            return {"src": image_path, "height": "400px"}  # Return image with height setting
        else: 
            return "❌ Gene not found in the dataset for custom UMAP."
        #return None  # Return None if image is not found

    @output
    @render.text
    def reprocess_status():
        return "Press 'Recalculate QC' to apply new metrics."

    @output
    @render.text
    def fileStatus_1():
        return "✅ Processed CSV is ready for download." if file_ready() else "❌ Processed CSV not found yet. Run the analysis first."
    
    @output
    @render.text
    def fileStatus_2():
        return "✅ Differential Gene Expression list is ready for download." if file_ready() else "❌ DEG List not found yet. Run the analysis first."


    @reactive.effect
    @reactive.event(input.reprocess_button)
    def reprocess_data():
        with open(config_path, "r") as f:
            config = json.load(f)
        config["preprocessing_params"]["min_genes"] = input.min_genes()
        config["preprocessing_params"]["min_cells"] = input.min_cells()
        config["preprocessing_params"]["target_sum"] = input.target_sum()
        with open(config_path, "w") as f:
            json.dump(config, f, indent=4)

        subprocess.run(["python", "/app/src/main.py"], check=True)
        file_ready.set(True)

        @output
        @render.text
        def reprocess_status():
            return "✅ Reprocessing complete!"

app = App(app_ui, server)

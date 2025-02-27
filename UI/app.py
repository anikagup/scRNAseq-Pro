from shiny import App, ui, render, reactive
import os
import subprocess
import shutil  

# Define the UI layout
app_ui = ui.page_fluid(
    ui.tags.style("body { background-color: lightblue; }"),
    ui.panel_title("Welcome to RNA Pro"),
    ui.input_file("file_upload", "Please Upload Your fastq File", multiple=False, accept=[".fastq"]),
    ui.output_text("file_info"),
    ui.output_ui("buttons"),  
    ui.output_image("plot")  
)

def server(input, output, session):
    upload_dir = os.path.join(os.getcwd(), 'uploads')  
    os.makedirs(upload_dir, exist_ok=True)  

    selected_plot = reactive.Value("")  # Track the last clicked button

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
        return f"File Saved: {saved_path}" if saved_path else "Error saving file."

    @output
    @render.ui
    def buttons():
        file = input.file_upload()
        if file is None:
            return None  

        return ui.div(
            ui.input_action_button("button1", "UMAP"),
            ui.input_action_button("button2", "TSNE")
        )

    # Ensure reactivity is correctly triggered
    @reactive.effect
    @reactive.event(input.button1)
    def on_umap_click():
        selected_plot.set("UMAP")

    @reactive.effect
    @reactive.event(input.button2)
    def on_tsne_click():
        selected_plot.set("TSNE")

    @output
    @render.image
    def plot():
        plot_type = selected_plot.get()
        if not plot_type:
            return None  # No plot until a button is clicked

        static_dir = os.path.join(os.getcwd(), 'static')
        os.makedirs(static_dir, exist_ok=True)  

        if plot_type == "UMAP":
            img_filename = 'UMAP.png'
            img_path = os.path.join(static_dir, img_filename)
            print(f"Generating UMAP plot at {img_path}")
            subprocess.run(['python', 'pretendUMAP.py', img_path])
            return {"src": f"static/{img_filename}", "height": "400px"}

        elif plot_type == "TSNE":
            img_filename = 'TSNE.png'
            img_path = os.path.join(static_dir, img_filename)
            print(f"Generating TSNE plot at {img_path}")
            subprocess.run(['python', 'pretendTSNE.py', img_path])
            return {"src": f"static/{img_filename}", "height": "400px"}

        return None  

# Create the Shiny app
app = App(app_ui, server)







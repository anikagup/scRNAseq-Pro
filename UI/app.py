from shiny import App, ui, render, reactive
import os
import subprocess
import matplotlib.pyplot as plt

# Define the UI layout
app_ui = ui.page_fluid(
    ui.tags.style("body { background-color: lightblue; }"),
    ui.panel_title("Welcome to RNA Pro"),
    ui.input_file("file_upload", "Please Upload Your fastq File", multiple=False, accept=[".fastq"]),
    ui.output_text("file_info"),
    ui.output_ui("buttons"),  # Move buttons above the plot
    ui.output_image("plot")   # Plot now comes after the buttons
)

def server(input, output, session):
    # Display file name
    @output
    @render.text
    def file_info():
        file = input.file_upload()
        if file is None:
            return "No file uploaded."
        return f"File Name: {file[0]['name']}"
    
    # Show buttons after file upload
    @output
    @render.ui
    def buttons():
        file = input.file_upload()
        if file is None:
            return None  # No buttons if no file is uploaded

        # Create two action buttons once the file is uploaded
        return ui.div(
            ui.input_action_button("button1", "UMAP"),
            ui.input_action_button("button2", "TSNE")
        )
    
    # Plot image display (only after Button 1 or Button 2 is clicked)
    @output
    @render.image
    def plot():
        static_dir = os.path.join(os.getcwd(), 'static')
        img_filename = 'UMAP.png'
        img_path = os.path.join(static_dir, img_filename)

        # Ensure the 'static/' directory exists
        os.makedirs(static_dir, exist_ok=True)

        # If Button 1 is clicked, generate plot from "pretendUMAP.py"
        if input.button1():
            if not os.path.exists(img_path):
                print(f"Generating UMAP plot at {img_path}")
                subprocess.run(['python', 'pretendUMAP.py', img_path])

            return {"src": f"static/{img_filename}", "height": "400px"}

        # If Button 2 is clicked, generate a plot from "generate_graph.py"
        elif input.button2():
            graph_filename = 'TSNE.png'
            graph_path = os.path.join(static_dir, graph_filename)

            if not os.path.exists(graph_path):
                print(f"Generating Graph plot at {graph_path}")
                subprocess.run(['python', 'pretendTSNE.py', graph_path])

            return {"src": f"static/{graph_filename}", "height": "400px"}

        return None  # No plot if no button is clicked

# Create the Shiny app
app = App(app_ui, server)













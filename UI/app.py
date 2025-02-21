from shiny import App, ui, render, reactive
import os
import subprocess
import shutil  # Import shutil to move/copy files
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
    upload_dir = os.path.join(os.getcwd(), 'uploads')  # Define a directory to store uploads
    os.makedirs(upload_dir, exist_ok=True)  # Ensure the directory exists

    # Function to save uploaded file
    def save_uploaded_file():
        file = input.file_upload()
        if file:
            temp_path = file[0]["datapath"]  # Get temporary file path
            saved_path = os.path.join(upload_dir, file[0]["name"])  # Define target path
            shutil.move(temp_path, saved_path)  # Move file to 'uploads/' directory
            return saved_path  # Return the new saved file path
        return None

    # Display file name and save file
    @output
    @render.text
    def file_info():
        file = input.file_upload()
        if file is None:
            return "No file uploaded."
        
        saved_path = save_uploaded_file()  # Save the file
        return f"File Saved: {saved_path}" if saved_path else "Error saving file."

    # Show buttons after file upload
    @output
    @render.ui
    def buttons():
        file = input.file_upload()
        if file is None:
            return None  # No buttons if no file is uploaded

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

        os.makedirs(static_dir, exist_ok=True)  # Ensure 'static/' exists

        if input.button1():
            if not os.path.exists(img_path):
                print(f"Generating UMAP plot at {img_path}")
                subprocess.run(['python', 'pretendUMAP.py', img_path])

            return {"src": f"static/{img_filename}", "height": "400px"}

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

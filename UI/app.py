from shiny import App, ui, render, reactive
import os
import subprocess
import matplotlib.pyplot as plt

# Define the UI layout
app_ui = ui.page_fluid(
    ui.h1("File Upload, Filtering, and Plotting"),
    ui.input_file("file_upload", "Upload File", multiple=False, accept=[".fastq"]),
    ui.output_text("file_info"),
    ui.output_text("filtered_lines"),
    ui.output_image("plot")
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

    # Display filtered sequence lines
    @output
    @render.text
    def filtered_lines():
        file = input.file_upload()
        if file is None:
            return ""

        try:
            file_path = file[0]['datapath']
            filtered = []

            # Read the file line-by-line to avoid loading everything into memory
            with open(file_path, 'r') as f:
                lines = f.readlines()

                # Process the lines: keep 2nd and 4th lines in each group of 4
                for i in range(0, len(lines), 4):
                    if i + 1 < len(lines):
                        filtered.append(lines[i + 1].strip())
                    if i + 3 < len(lines):
                        filtered.append(lines[i + 3].strip())

            return '\n'.join(filtered[:6])

        except Exception as e:
            return f"Error processing file: {e}"

    # Plot image display
    @output
    @render.image
    def plot():
        static_dir = os.path.join(os.getcwd(), 'static')
        img_filename = 'plot_image.png'
        img_path = os.path.join(static_dir, img_filename)

        # Ensure the 'static/' directory exists
        os.makedirs(static_dir, exist_ok=True)

        # Check if the image file exists or generate it
        if not os.path.exists(img_path):
            print(f"Generating plot at {img_path}")
            subprocess.run(['python', 'generate_graph.py', img_path])

        # Return the relative URL path (Shiny serves 'static/' folder automatically)
        return {"src": f"static/{img_filename}", "height": "400px"}

# Create the Shiny app
app = App(app_ui, server)







        






















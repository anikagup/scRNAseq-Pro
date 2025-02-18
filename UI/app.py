from shiny import App, ui, render
import os
import subprocess

# Define the UI layout
app_ui = ui.page_fluid(
    ui.h1("Generated Graph"),
    ui.output_image("plot")
)

def server(input, output, session):
    @output
    @render.image
    def plot():
        # Get the absolute path for the 'static/' folder if present
        static_dir = os.path.join(os.getcwd(), 'static')
        img_path = os.path.join(static_dir, 'plot_image.png')

        # Ensure the 'static/' directory exists
        os.makedirs(static_dir, exist_ok=True)

        # Check if the image file exists
        if not os.path.exists(img_path):
            # If the image doesn't exist, call the external script to generate it
            print(f"Calling generate_graph.py with path: {img_path}")  # Debugging line to check the path
            subprocess.run(['python', 'generate_graph.py', img_path])  # Run the external script

        # Return the image path to display it
        return {'src': img_path, 'height': '400px'}

# Create the Shiny app
app = App(app_ui, server)







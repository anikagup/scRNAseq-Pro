from shiny import App, ui, render, reactive
import os
import shutil
import subprocess

# Define the UI layout
app_ui = ui.page_fluid(
    ui.tags.style("body { background-color: lightblue; }"),
    ui.panel_title("Welcome to RNA Pro"),
    ui.input_file("file_upload", "Please Upload Your FASTQ File", multiple=False, accept=[".fastq"]),
    ui.output_text("file_info"),
    ui.output_ui("buttons"),
    ui.output_ui("axis_inputs"),  # Inputs for axis customization
    ui.output_ui("generate_button_ui"),  # Generate plot button
    ui.output_image("plot")  # Display generated plot
)

def server(input, output, session):
    upload_dir = os.path.join(os.getcwd(), 'uploads')  
    os.makedirs(upload_dir, exist_ok=True)  

    selected_plot = reactive.Value("")  # Track selected plot type
    is_plotting = reactive.Value(False)  # Prevent multiple plot generations
    x_max_value = reactive.Value(10)  
    xmin_value = reactive.Value(0)  
    y_max_value = reactive.Value(1)  
    y_min_value = reactive.Value(-1)  
    plot_ready = reactive.Value(False)  
    uploaded_file_path = reactive.Value(None)  # Store uploaded file path

    def save_uploaded_file():
        file = input.file_upload()
        if file:
            temp_path = file[0]["datapath"]  
            saved_path = os.path.join(upload_dir, file[0]["name"])  
            shutil.move(temp_path, saved_path)  
            return saved_path  # Return the saved file path
        return None

    @output
    @render.text
    def file_info():
        file = input.file_upload()
        if file is None:
            return "No file uploaded."
        
        saved_path = save_uploaded_file()  
        if saved_path:
            uploaded_file_path.set(saved_path)  # Store file path in reactive value
            return f"File Saved: {saved_path}"
        return "Error saving file."

    @output
    @render.ui
    def buttons():
        file = input.file_upload()
        if file is None:
            return None  

        return ui.div(
            ui.h4("Choose plot type"),
            ui.input_action_button("button1", "UMAP"),
            ui.input_action_button("button2", "TSNE")
        )

    @output
    @render.ui
    def axis_inputs():
        if selected_plot.get():
            return ui.div(
                ui.input_numeric("x_max", "Enter X Max Value", value=x_max_value.get(), min=1, step=1),
                ui.input_numeric("xmin", "Enter X Min Value", value=xmin_value.get(), min=0, step=1),
                ui.input_numeric("y_max", "Enter Y Max Value", value=y_max_value.get(), min=0, step=1),
                ui.input_numeric("y_min", "Enter Y Min Value", value=y_min_value.get(), min=-10, step=1)
            )
        return None

    @output
    @render.ui
    def generate_button_ui():
        if selected_plot.get():
            return ui.input_action_button("generate_button", "Generate Plot")
        return None  

    @reactive.effect
    @reactive.event(input.button1)
    def on_umap_click():
        selected_plot.set("UMAP")
        plot_ready.set(False)  

    @reactive.effect
    @reactive.event(input.button2)
    def on_tsne_click():
        selected_plot.set("TSNE")
        plot_ready.set(False)  

    @reactive.effect
    @reactive.event(input.x_max)
    def on_x_max_change():
        x_max_value.set(input.x_max())

    @reactive.effect
    @reactive.event(input.xmin)
    def on_xmin_change():
        xmin_value.set(input.xmin())

    @reactive.effect
    @reactive.event(input.y_max)
    def on_y_max_change():
        y_max_value.set(input.y_max())

    @reactive.effect
    @reactive.event(input.y_min)
    def on_y_min_change():
        y_min_value.set(input.y_min())

    @reactive.effect
    @reactive.event(input.generate_button)
    def on_generate_button_click():
        if selected_plot.get():
            start_plot_generation(selected_plot.get())

    def start_plot_generation(plot_type):
        if is_plotting.get():
            return  

        file_path = uploaded_file_path.get()
        if not file_path:
            print("No file uploaded for processing.")
            return

        current_x_max = x_max_value.get()
        current_xmin = xmin_value.get()
        current_y_max = y_max_value.get()
        current_y_min = y_min_value.get()
        
        is_plotting.set(True)
        plot_ready.set(False)  

        generate_plot(plot_type, file_path, current_x_max, current_xmin, current_y_max, current_y_min)

    def generate_plot(plot_type, file_path, x_max, xmin, y_max, y_min):
        static_dir = os.path.join(os.getcwd(), 'static')
        os.makedirs(static_dir, exist_ok=True)

        img_filename = f'{plot_type}.png'
        img_path = os.path.join(static_dir, img_filename)

        print(f"Generating {plot_type} plot for {file_path} at {img_path} with X Max {x_max}, X Min {xmin}, Y Max {y_max}, Y Min {y_min}")

        if plot_type == "UMAP":
            subprocess.run(['python', 'pretendUMAP.py', img_path, file_path, str(x_max), str(xmin), str(y_max), str(y_min)])
        elif plot_type == "TSNE":
            subprocess.run(['python', 'pretendTSNE.py', img_path, file_path, str(x_max), str(xmin), str(y_max), str(y_min)])

        plot_ready.set(True)  
        is_plotting.set(False)

    @output
    @render.image
    def plot():
        plot_type = selected_plot.get()

        if not plot_ready.get():
            return None  

        static_dir = os.path.join(os.getcwd(), 'static')
        img_filename = f'{plot_type}.png'
        img_path = os.path.join(static_dir, img_filename)

        print(f"Image path: {img_path}")

        if os.path.exists(img_path):
            return {"src": f"static/{img_filename}", "height": "400px"}

        return None  

app = App(app_ui, server)

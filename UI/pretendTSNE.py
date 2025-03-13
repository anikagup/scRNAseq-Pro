import sys
import matplotlib.pyplot as plt
import numpy as np

def generate_tsne_plot(img_path, x_max, xmin, y_max, y_min):
    # Simulate TSNE plot generation using the input values

    # Generate sample data for the plot (for demonstration purposes)
    x = np.linspace(xmin, x_max, 100)
    y = np.cos(x) * y_max  # Apply y_max to scale the y-values
    y = np.clip(y, y_min, y_max)  # Apply y_min to clip the values to a specific range

    fig, ax = plt.subplots()
    ax.plot(x, y)
    ax.set_title("Pretend this is a TSNE XD")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True)

    # Save the plot to the specified path
    fig.savefig(img_path)
    plt.close(fig)

if __name__ == "__main__":
    # The arguments are passed from the Shiny app
    img_path = sys.argv[1]
    x_max = float(sys.argv[2])
    xmin = float(sys.argv[3])
    y_max = float(sys.argv[4])
    y_min = float(sys.argv[5])

    generate_tsne_plot(img_path, x_max, xmin, y_max, y_min)







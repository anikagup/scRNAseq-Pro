import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def generate_plot(img_path):
    # Ensure the directory exists
    os.makedirs(os.path.dirname(img_path), exist_ok=True)

    # Generate the plot
    x = np.linspace(0, 10, 100)
    y = np.sin(x)

    fig, ax = plt.subplots()
    ax.plot(x, y)
    ax.set_title("y = sin(x)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True)

    # Save the plot to the specified file
    fig.savefig(img_path)
    plt.close(fig)

if __name__ == "__main__":
    img_path = sys.argv[1]  # Path passed from Shiny app
    print(f"Saving image to {img_path}")  # Debugging line to check the path
    generate_plot(img_path)



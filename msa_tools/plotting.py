import matplotlib.pyplot as plt
import numpy as np


def plot_sliding(results, output=None, title="Sliding window"): #NEEDS TESTING
    """
    Plot sliding window metric.

    results: list of (start, end, value)
    """
    positions = [(s + e) / 2 for s, e, _ in results]
    values = [v for _, _, v in results]

    plt.figure(figsize=(10, 4))
    plt.plot(positions, values)
    plt.xlabel("Genome position")
    plt.ylabel("Metric value")
    plt.title(title)
    plt.tight_layout()

    if output:
        plt.savefig(output, dpi=300)
    else:
        plt.show()

    plt.close()


def plot_distance_matrix(names, matrix, output=None, title="Distance matrix", fontsize=10):
    """
    Plot heatmap of distance matrix
    """
    plt.figure(figsize=(6, 5))
    im = plt.imshow(matrix)
    plt.colorbar(im)

    plt.xticks(range(len(names)), names, rotation=90, fontsize=fontsize)
    plt.yticks(range(len(names)), names, fontsize=fontsize)

    plt.title(title)
    plt.tight_layout()

    if output:
        plt.savefig(output, dpi=300)
    else:
        plt.show()

    plt.close()
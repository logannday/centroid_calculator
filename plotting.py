import matplotlib.pyplot as plt
import seaborn as sns

def plot_centroid_vs_rmsd(df, title="Centroid Distance vs Interface RMSD", save_path=None):
    """
    Generate a scatter plot of centroid_distance vs interface_rmsd.

    Parameters:
        df (pd.DataFrame): DataFrame containing 'centroid_distance' and 'interface_rmsd' columns.
        title (str): Plot title.
        save_path (str or None): Path to save the figure (e.g., 'plot.png'). If None, does not save.

    Returns:
        None
    """
    if 'centroid_distance' not in df.columns or 'interface_rmsd' not in df.columns:
        raise ValueError("DataFrame must contain 'centroid_distance' and 'interface_rmsd' columns.")

    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=df, x='centroid_distance', y='interface_rmsd', s=60, edgecolor='w', alpha=0.7)
    plt.title(title)
    plt.xlabel("Centroid Distance")
    plt.ylabel("Interface RMSD")
    plt.grid(True)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300)
    
    plt.show()

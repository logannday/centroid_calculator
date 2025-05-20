import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os



def plot_per_residue(df, output_dir="plots", pdb_id="protein"):
    """
    Generate separate scatter plots of centroid_distance vs interface_rmsd for each unique residue.

    Parameters:
        df (pd.DataFrame): DataFrame containing 'residue', 'centroid_distance', and 'interface_rmsd' columns.
        output_dir (str): Directory where plots will be saved.
        pdb_id (str): Identifier used in filenames to distinguish plots.

    Returns:
        None
    """
    required_cols = {"residue", "centroid_distance", "interface_rmsd"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"DataFrame must contain columns: {required_cols}")

    os.makedirs(output_dir, exist_ok=True)

    # Plot one figure per unique residue
    for residue, group in df.groupby("residue"):
        plt.figure(figsize=(8, 6))
        sns.scatterplot(
            data=group,
            x="centroid_distance",
            y="interface_rmsd",
            s=60,
            edgecolor="w",
            alpha=0.7,
        )
        plt.title(f"{pdb_id.upper()} - {residue}: Centroid vs RMSD")
        plt.xlabel("Centroid Distance")
        plt.ylabel("Interface RMSD")
        plt.grid(True)
        plt.tight_layout()

        filename = f"{pdb_id}_{residue}_plot.png".replace(" ", "_")
        plt.savefig(os.path.join(output_dir, filename), dpi=300)
        plt.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--filepath', type=str, required=True, help="Path to input CSV file")
    parser.add_argument('--pdb_id', type=str, required=True, help="PDB ID of the protein")
    parser.add_argument('--output_dir', type=str, default="plots", help="Directory to save plots")
    args = parser.parse_args()

    df = pd.read_csv(args.filepath)
    plot_per_residue(df, output_dir=args.output_dir, pdb_id=args.pdb_id)

if __name__ == "__main__":
    main()


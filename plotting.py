import argparse, os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def plot_centroid_vs_rmsd(df, pdb_id):
    """
    Generate a scatter plot of centroid_distance vs interface_rmsd.

    Parameters:
        df (pd.DataFrame): DataFrame containing 'centroid_distance' and 'interface_rmsd' columns.
        title (str): Plot title.
        save_path (str or None): Path to save the figure (e.g., 'plot.png'). If None, does not save.

    Returns:
        None
    """
    if "centroid_distance" not in df.columns or "interface_rmsd" not in df.columns:
        raise ValueError(
            "DataFrame must contain 'centroid_distance' and 'interface_rmsd' columns."
        )

    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        data=df,
        x="centroid_distance",
        y="interface_rmsd",
        s=60,
        edgecolor="w",
        alpha=0.7,
    )
    plt.title("Distance from centroid to inserted residue vs Interface RMSD")
    plt.xlabel("Centroid Distance")
    plt.ylabel("Interface RMSD")
    plt.grid(True)
    plt.tight_layout()

    plt.savefig(f"plots/{pdb_id}/Centroid_vs_rmsd", dpi=300)

    plt.show()
    import pandas as pd


def plot_rmsd_and_angle_by_location(df, pdb_id):
    """
    Takes in a DataFrame with columns 'rmsd', 'angle', and 'at_near_away',
    and creates boxplots for rmsd and angle grouped by at_near_away.
    """

    # Ensure required columns exist
    required_cols = {'interface_rmsd', 'angle', 'at_near_away'}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"DataFrame must contain columns: {required_cols}")
    
    # Set style
    sns.set_theme(style="whitegrid")

    # Plot RMSD vs at/near/away
    plt.figure(figsize=(12, 5))

    plt.subplot(2, 2, 1)
    sns.boxplot(data=df, x='at_near_away', y='interface_rmsd', palette="Set2")
    plt.title("RMSD by Mutation Location")
    plt.xlabel("Mutation Location")
    plt.ylabel("RMSD Angstroms")

    # Plot Angle vs at/near/away
    plt.subplot(2, 2, 2)
    sns.boxplot(data=df, x='at_near_away', y='angle', palette="Set3")
    plt.title("Angle by Mutation Location")
    plt.xlabel("Mutation Location")
    plt.ylabel("Angle (degrees)")

    # Plot Angle vs at/near/away
    plt.subplot(2, 2, 3)
    sns.boxplot(data=df, x='at_near_away', y='delta_SASA', palette="Set3")
    plt.title("∆SASA by Mutation Location")
    plt.xlabel("Mutation Location")
    plt.ylabel("∆SASA Angstroms")
    plt.title("∆SASA by Mutation Location")
    plt.xlabel("Mutation Location")

    plt.tight_layout()
    plt.savefig(f"plots/{pdb_id}/away_near_at", dpi=300)

pdb_id = None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--filepath', type = str, help = "path to csv")
    parser.add_argument('--pdb_id', type = str, help = "pdb id of protien")
    args = parser.parse_args()

    pdb_id = args.pdb_id
    os.makedirs(f"./plots/{pdb_id}", exist_ok=True)
    filepath = args.filepath
    df = pd.read_csv(filepath)
    # plot_centroid_vs_rmsd(df, pdb_id)
    plot_rmsd_and_angle_by_location(df, pdb_id)

if __name__ == "__main__":
    main()

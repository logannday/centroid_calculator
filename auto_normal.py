from pymol2 import PyMOL
import numpy as np
from pymol.cgo import BEGIN, COLOR, TRIANGLES, END, VERTEX, ALPHA


def draw_plane(
    cmd,
    center,
    normal,
    plane_width=20.0,
    alpha_val=0.5,
    color_val=(0.8, 0.3, 0.8),
    plane_name="Plane",
):
    """
    Draw a square plane in PyMOL using CGO

    Args:
        cmd: PyMOL cmd object
        center (np.array): center point of the plane
        normal (np.array): normal vector (doesn't have to be normalized)
        plane_width (float): width of the square plane
        alpha_val (float): transparency (0 to 1)
        color_val (tuple): RGB color (0 to 1)
        plane_name (str): name of the CGO object
    """
    normal = normal / np.linalg.norm(normal)

    # Choose vector not parallel to normal
    if abs(normal[2]) < 0.9:
        temp = np.array([0, 0, 1])
    else:
        temp = np.array([0, 1, 0])

    u = np.cross(normal, temp)
    u /= np.linalg.norm(u)
    v = np.cross(normal, u)
    v /= np.linalg.norm(v)

    half = plane_width / 2
    corner1 = center + half * u + half * v
    corner2 = center - half * u + half * v
    corner3 = center - half * u - half * v
    corner4 = center + half * u - half * v

    plane = [
        ALPHA,
        alpha_val,
        BEGIN,
        TRIANGLES,
        COLOR,
        *color_val,
        VERTEX,
        *corner1,
        VERTEX,
        *corner2,
        VERTEX,
        *corner3,
        VERTEX,
        *corner3,
        VERTEX,
        *corner4,
        VERTEX,
        *corner1,
        END,
    ]

    cmd.load_cgo(plane, plane_name)


def plane_fig(pdb_fp, center1, normal1, center2, normal2, fig_name):
    with PyMOL() as pymol:
        cmd = pymol.cmd
        cmd.load(f"/research/jagodzinski/interface_mutations/complete_runs/MUTANT_PDBS/5m2o/{pdb_fp}")
        cmd.bg_color("white")
        cmd.color("deepsalmon", "chain A")
        cmd.color("deepteal", "chain B")

        # Plane 1 parameters
        draw_plane(
            cmd,
            center1,
            normal1,
            plane_width=20.0,
            alpha_val=0.5,
            color_val=(0.8, 0.3, 0.8),
            plane_name="Plane_A",
        )

        draw_plane(
            cmd,
            center2,
            normal2,
            plane_width=15.0,
            alpha_val=0.4,
            color_val=(0.2, 0.6, 0.4),
            plane_name="Plane_B",
        )

        view_matrix = (
            -0.037086401134729385,
            -0.7610383629798889,
            -0.6476132869720459,
            -0.3061339557170868,
            0.6255485415458679,
            -0.7175810933113098,
            0.9512314200401306,
            0.1716453582048416,
            -0.25619155168533325,
            0.0,
            0.0,
            -168.339599609375,
            -0.8642916679382324,
            -5.857709884643555,
            35.0895881652832,
            152.60862731933594,
            184.07057189941406,
            -20.0,
        )

        cmd.set_view(view_matrix)
        cmd.zoom()
        # Ray trace at 4K resolution
        cmd.ray(3840, 2160)

        # Save PNG with 4K resolution and high dpi (dpi affects print size, not pixel dimensions)
        cmd.png(f"./plane_figs/{fig_name}", dpi=300)


if __name__ == "__main__":
    center1 = np.array([1.1541645e-02, -6.7973328e00, 4.3206039e01])
    normal1 = np.array([-0.14236136, -0.46168193, 0.8755473])
    center2 = center1 + np.array([5, 0, 0])
    normal2 = np.array([0.3, -0.7, 0.65])

    plane_fig("mut.pdb", center1, normal1, center2, normal2)

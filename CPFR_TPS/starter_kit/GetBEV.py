#!/usr/bin/env python3
########################################################################
########################################################################
# Nov 2025 A.Burattini
########################################################################
########################################################################

import argparse
import numpy as np
from mhd_io import *
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from matplotlib.patches import Polygon

import numpy as np

def rotate_slit(coords, angle_deg):
    """
    Ruota un vettore 2D (o una lista N×2 di vettori) 
    di angle_deg gradi in senso antiorario.
    Se angle_deg è negativo ⇒ rotazione oraria.
    """
    theta = np.deg2rad(angle_deg)

    R = np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta),  np.cos(theta)]
    ])

    coords = np.asarray(coords)

    # singolo punto
    if coords.ndim == 1:
        return R @ coords

    # lista di punti
    return coords @ R.T


def plot_filled_ptv_contour(ptv_coords, ax=None, facecolor="red", edgecolor="black", alpha=0.3, show_points=False):
    """
    ptv_coords: array di shape (N, 2), coordinate [x, y] in cm
    ax: oggetto matplotlib.axes; se None ne crea uno nuovo
    facecolor: colore di riempimento
    edgecolor: colore del bordo
    alpha: trasparenza del riempimento
    show_points: se True mostra anche i punti originali
    """
    ptv_coords = np.asarray(ptv_coords)

    if ptv_coords.shape[1] != 2:
        raise ValueError("ptv_coords deve avere shape (N, 2) con coordinate [x, y].")

    # Calcola l'involucro convesso (convex hull)
    hull = ConvexHull(ptv_coords)

    # Prendi i vertici nel giusto ordine
    hull_points = ptv_coords[hull.vertices]

    # Crea axes se necessario
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))

    # Riempie il contorno
    ax.fill(hull_points[:, 0], hull_points[:, 1],
            facecolor=facecolor, edgecolor=edgecolor, alpha=alpha)

    # (Opzionale) mostra anche i punti
    if show_points:
        ax.scatter(ptv_coords[:, 0], ptv_coords[:, 1], s=5, color="k")

    return ax


parser = argparse.ArgumentParser(
    description='Get the size of a dose grid keeping the same spacing of the CT'
)
parser.add_argument("-ptv",   help="CT map",   required=True)
parser.add_argument("-slitsize", help="size of the slits", type=float, nargs=2, required=True)
parser.add_argument("-W",  help="horizontal field size", type=float, required=True)
parser.add_argument("-H",  help="vertical field size", type=float, required=True)
parser.add_argument("-angle",  help="shape angle", type=float, required=True)
parser.add_argument("-showsave", help="0 to show or 1 to save image, else for none", type=int, required=True)

args = parser.parse_args()

ptv=args.ptv
shortslit=np.min(np.array(args.slitsize))
longslit=np.max(np.array(args.slitsize))
W=args.W
H=args.H
angle=args.angle
showsave=args.showsave
[nn, hs, x0, Map]=unpackVoxels(mhd_read(ptv))

print(nn, "\n", hs, "\n",  x0, "\n")

ptv_proj=np.max(np.array(Map), axis=2)
#plt.imshow(ptv_proj, cmap="gray", interpolation="nearest")
#plt.show()

ptv_idx=np.array(list(zip(*np.where(ptv_proj))))
#print(ptv_idx[:,0])
ptv_coords=ptv_idx.astype(float)
ptv_coords[:,0]=ptv_idx[:,0]*hs[0]+x0[0]
ptv_coords[:,1]=ptv_idx[:,1]*hs[1]+x0[1]
#print(ptv_coords)

ax = plot_filled_ptv_contour(ptv_coords, facecolor="red", edgecolor="black", alpha=0.4, show_points=False)
ax.set_aspect("equal")
ax.set_xlabel("x [cm]")
ax.set_ylabel("y [cm]")
ax.set_title("PTV contour")
#ax.plot([1, 1], [-1, 1], color="blue", linewidth=2)

upslit=[(-longslit/2, H/2), (-longslit/2, H/2+shortslit), (longslit/2, H/2+shortslit), (longslit/2, H/2)]
downslit=[(-longslit/2, -H/2), (-longslit/2, -H/2-shortslit), (longslit/2, -H/2-shortslit), (longslit/2, -H/2)]
leftslit=[(-W/2, -longslit/2), (-W/2-shortslit, -longslit/2), (-W/2-shortslit, longslit/2), (-W/2, longslit/2)]
rightslit=[(W/2, -longslit/2), (W/2+shortslit, -longslit/2), (W/2+shortslit, longslit/2), (W/2, longslit/2)]


verts = [upslit, downslit, leftslit, rightslit]
for v in verts:
  rot_v=rotate_slit(v, -angle)
  poly = Polygon(rot_v, closed=True, facecolor="grey", edgecolor="black", alpha=0.5)
  ax.add_patch(poly)

p=np.array(upslit + downslit + leftslit + rightslit)
text_coords=[np.max(p[:,0])+1, np.max(p[:,1])+1]
plt.text(text_coords[0], text_coords[1],
    f"SHAPER ANGLE: {angle}\u00B0 \nFIELD WIDTH: {W}cm \nFIELD HEIGHT: {H}cm",      # testo
    fontsize=8,
    bbox=dict(
        facecolor="white",
        edgecolor="black",
        boxstyle="square,pad=0.5"
    )
)


if showsave==0:
  plt.show()
elif showsave==1:
  plt.savefig("BEV.jpg", dpi=200, bbox_inches="tight")
  print("BEV.jpg saved")
else:
  print("image is neither shown or saved")

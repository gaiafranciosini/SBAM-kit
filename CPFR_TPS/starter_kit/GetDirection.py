import SimpleITK as sitk
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from skimage import measure
import argparse
from mhd_io import *
import scipy.ndimage as ndi

def loadCT(CT_path):
    [nn,hs,xmin,Map] = unpackVoxels(mhd_read(CT_path))
    voxels = mhd_read(CT_path)
    return voxels, nn, hs, xmin, Map

def generate_closed_body_mask(ct_array, hu_threshold=-500, closing_radius=3, plot=False):
    binary = ct_array > hu_threshold
    structure = ndi.generate_binary_structure(3, 2)
    closed = ndi.binary_closing(binary, structure=structure, iterations=closing_radius)
    filled = np.zeros_like(closed)
    for i in range(ct_array.shape[2]):
        filled[:,:,i] = ndi.binary_fill_holes(closed[:,:,i])

    if plot:
        mid_z = ct_array.shape[2] // 2
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        axes[0].imshow(binary[:,:,mid_z], cmap='gray')
        axes[0].set_title('Maschera iniziale (HU > soglia)')
        axes[1].imshow(closed[:,:,mid_z], cmap='gray')
        axes[1].set_title('Dopo binary_closing')
        axes[2].imshow(filled[:,:,mid_z], cmap='gray')
        axes[2].set_title('Dopo fill_holes (finale)')
        for ax in axes:
            ax.axis('off')
        plt.tight_layout()
        plt.show()

    return filled.astype(np.uint8)

def get_skin(body_mask, spacing, shell_thickness_cm=0.1, plot=False):
    inner_dist = ndi.distance_transform_edt(body_mask, sampling=spacing)
    shell_mask = np.logical_and(body_mask == 1, inner_dist <= shell_thickness_cm)
    shell_mask = shell_mask.astype(np.uint8)

    if plot:
        mid_z = body_mask.shape[2] // 2
        plt.imshow(shell_mask[:,:,mid_z], cmap='gray')
        plt.title('Shell di {} cm sul contorno del corpo'.format(shell_thickness_cm))
        plt.axis('off')
        plt.show()

    return shell_mask

def extract_voxels_in_sphere(skin, rep_coords, spacing, origin, radius):
    sphere_voxels = []
    for x in range(skin.shape[0]):
        for y in range(skin.shape[1]):
            for z in range(skin.shape[2]):
                if skin[x, y, z] > 0:
                    voxel_coords = np.array([x, y, z]) * spacing + origin
#                    print(voxel_coords)
                    distance = np.linalg.norm(voxel_coords - np.array(rep_coords))
#                    print(distance)
                    if distance <= radius:
                        sphere_voxels.append(voxel_coords)
    return np.array(sphere_voxels)

def set_verse(body_mask_array, origin, spacing, centroid, front, max_distance_mm=500, step_mm=1):
    spacing=spacing*10
    origin=origin*10
    centroid_mm = centroid * 10
    steps = int(max_distance_mm / step_mm)

    def distance_to_exit(direction):
        for i in range(1, steps + 1):
            point_mm = centroid_mm + i * step_mm * direction
            idx = np.round((point_mm - origin) / spacing).astype(int)
            z, y, x = idx[2], idx[1], idx[0]
            if not (0 <= z < body_mask_array.shape[2] and 0 <= y < body_mask_array.shape[1] and 0 <= x < body_mask_array.shape[0]):
                return i * step_mm
            if body_mask_array[x, y, z] == 0:
                return i * step_mm
        return max_distance_mm + 1

    dist_forward = distance_to_exit(front)
    dist_backward = distance_to_exit(-front)

    if(dist_forward < dist_backward):
        return front
    else:
        return -front

def get_centroid(array, spacing, origin):
    array=np.transpose(array, (2,1,0))
    labels = (array > 0).astype(np.uint8)
    props = measure.regionprops(labels)
    if not props:
        raise ValueError(f"Nessuna regione trovata")
    centroid_voxel = props[0].centroid  # (z, y, x)
    z, y, x = centroid_voxel
    cx = origin[0] + x * spacing[0]
    cy = origin[1] + y * spacing[1]
    cz = origin[2] + z * spacing[2]
    return np.array([cx, cy, cz])

def get_ptv_coords(ptv_array, spacing, origin):
    ptv_idx = np.argwhere(ptv_array > 0)
    phys_coords = np.array([[origin[0] + x[0]*spacing[0],origin[1] + x[1]*spacing[1], origin[2] + x[2]*spacing[2]] for x in ptv_idx]) 
    return phys_coords

def fit_plane(points, centroid):
    centered_points = points - centroid
    _, _, vh = np.linalg.svd(centered_points)
    normal = vh[2, :]
    return normal

def compute_beam_vectors(normal):
    front = normal / np.linalg.norm(normal)
    arbitrary = np.array([1, 0, 0]) if abs(front[0]) < 0.9 else np.array([0, 1, 0])
    up = np.cross(front, arbitrary)
    up /= np.linalg.norm(up)
    left = np.cross(up, front)
    left /= np.linalg.norm(left)
    return front, up, left

def plot_3d(ptv_coords, sphere_coords, skin_coords, direction, rep_coords):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(ptv_coords[:, 0], ptv_coords[:, 1], ptv_coords[:, 2], s=1, label='PTV')
    ax.scatter(sphere_coords[:, 0], sphere_coords[:, 1], sphere_coords[:, 2], s=1, label='Sfera')
    ax.scatter(skin_coords[:, 0], skin_coords[:, 1], skin_coords[:, 2], s=1, alpha=0.3, label='Skin')
    centroid = rep_coords
    ax.quiver(centroid[0], centroid[1], centroid[2], direction[0], direction[1], direction[2], length=20, color='r', linewidth=2, label='Direzione')
#    ax.quiver(centroid[0], centroid[1], centroid[2], -direction[0], -direction[1], -direction[2], length=20, color='r', linewidth=2, label='min Direzione')
    ax.set_title('Visualizzazione 3D PTV e Sfera con Direzione')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.legend()
    plt.tight_layout()
    plt.show()

def main():
    parser = argparse.ArgumentParser(description='Get the beam direction from the CT repery')
    parser.add_argument('ct_path', help='Percorso della CT')
    parser.add_argument('-PTV', required=True, help='Percorso del PTV')
    parser.add_argument('-marker', required=True, help='Percorso file marker')
    parser.add_argument('-plot', action='store_true', help='Visualizza il grafico 3D')
    args = parser.parse_args()

    voxels, nn, spacing, origin, Map = loadCT(args.ct_path)
#    mid_z = Map.shape[2] // 2
#    plt.imshow(Map[:,:,mid_z], cmap='gray')
    Map=np.transpose(np.array(sitk.GetArrayFromImage(sitk.ReadImage(args.ct_path))), (2,1,0))
    ptv_array=np.transpose(np.array(sitk.GetArrayFromImage(sitk.ReadImage(args.PTV))), (2,1,0))
#    print(ptv_array.shape)
    marker_array=np.transpose(np.array(sitk.GetArrayFromImage(sitk.ReadImage(args.marker))), (2,1,0))
    rep_coords=get_centroid(marker_array, spacing, origin)
    print("MARKER coords [cm]:", rep_coords)
    iso_ptv=get_centroid(ptv_array, spacing, origin)
    ptv_coords=get_ptv_coords(ptv_array, spacing, origin)
    body_mask_array=generate_closed_body_mask(Map, hu_threshold=-500, closing_radius=3)
    skin_mask=get_skin(body_mask_array, spacing, shell_thickness_cm=0.2, plot=False)

#    mid_z = skin_mask.shape[2] // 2
#    plt.imshow(skin_mask[:,:,mid_z], cmap='gray')
#    plt.show()

    print(skin_mask.shape)
    points_in_sphere=extract_voxels_in_sphere(skin_mask, rep_coords, spacing, origin, 2)
    print(points_in_sphere.shape)
    skin_coords=get_ptv_coords(skin_mask, spacing, origin)
#    ptv_and_sphere=np.append(ptv_coords, points_in_sphere, axis=0)
#    normal=[1,0,0]
    normal=fit_plane(points_in_sphere, rep_coords) 
    front, up, left = compute_beam_vectors(normal)
    front=set_verse(body_mask_array, origin, spacing, iso_ptv, front, max_distance_mm=250, step_mm=3) 
    front, up, left = compute_beam_vectors(front)

    print("--DIRECTION ORTHOGONAL TO SURFACE--")
    print("BEAM DIRECTION:  "+str(front[0])+" "+str(front[1])+" "+str(front[2]))
#    print("UP", up)
#    print("LEFT", left)

    if args.plot:
        plot_3d(ptv_coords, points_in_sphere, skin_coords, front, rep_coords)

if __name__ == "__main__":
    main()

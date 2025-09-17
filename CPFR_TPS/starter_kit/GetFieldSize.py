import numpy as np
import sys
import SimpleITK as sitk
import matplotlib.pyplot as plt
from skimage.measure import find_contours
import argparse
from mhd_io import *

def load_ptv(ptv_path):
    voxels = mhd_read(ptv_path)
    [dim, spacing, offset, Map] = unpackVoxels(voxels)
    return Map, spacing, offset

def calculate_isocenter(ptv_array, spacing, origin):
    coords = np.array(np.nonzero(ptv_array)).T
    physical_coords = coords * spacing + origin
    isocenter = physical_coords.mean(axis=0)
    return isocenter

def find_contour_projection(ptv_array):
    projection = np.max(ptv_array, axis=2)
    contours = find_contours(projection, level=0.5)
    if contours:
        return projection, contours[0]
    return projection, None

def rotate_contour(contour, angle, center):
    angle_rad = np.radians(angle)
    rotation_matrix = np.array([
        [np.cos(angle_rad), -np.sin(angle_rad)],
        [np.sin(angle_rad), np.cos(angle_rad)]
    ])

    centered_contour = contour - center
    rotated_contour = np.dot(centered_contour, rotation_matrix.T) + center
    return rotated_contour

def calculate_rectangle(contour, spacing, origin, isocenter):
    physical_contour = contour * spacing[:2]+origin[:2]
    min_coords = physical_contour.min(axis=0)
    max_coords = physical_contour.max(axis=0)

    bottom_left = min_coords
    top_right = max_coords
    width = max_coords[0] - min_coords[0]
    height = max_coords[1] - min_coords[1]
    print("isocenter: ", isocenter)
    print("bottom_left: ", bottom_left)
    print("top_right: ", top_right)

    relative_coords = {
        "left": isocenter[0] - bottom_left[0],
        "right": top_right[0] - isocenter[0],
        "bottom": isocenter[1] - bottom_left[1],
        "top": top_right[1] - isocenter[1],
    }

    rectangle = {
        "bottom_left": bottom_left,
        "top_right": top_right,
        "width": width,
        "height": height,
        "relative_coords": relative_coords,
        "area": width * height,
    }
    
    return rectangle

def find_best_rectangle(projection, contour, spacing, origin, isocenter, iso_idx, rotation_angles):
    ptv_area = np.sum(projection) * spacing[0] * spacing[1]
    rectangles = []

    for angle in rotation_angles:
        rotated_contour = rotate_contour(contour, angle, iso_idx[:2])
        rectangle = calculate_rectangle(rotated_contour, spacing, origin, isocenter)
        rectangle["angle"] = angle
        rectangle["area_diff"] = abs(rectangle["area"] - ptv_area)
        rectangles.append(rectangle)

    best_rectangle = min(rectangles, key=lambda x: x["area_diff"])
    return rectangles, best_rectangle

def plot_results(projection, contour, rectangles, best_rectangle, isocenter, origin, spacing, rotation_angles):
    isocenter_voxel = (isocenter[:2]-origin[:2]) / spacing[:2]

    plt.figure(figsize=(10, 10), facecolor="white")
    #plt.imshow(projection.T, cmap="binary", origin="lower")

    # Disegna il contorno del PTV originale
    if contour is not None:
        plt.fill(contour[:, 0], contour[:, 1], color='red', alpha=0.2, label="PTV Region")
        min_x, max_x = contour[:, 0].min(), contour[:, 0].max()
        min_y, max_y = contour[:, 1].min(), contour[:, 1].max()
        margin = 10  # voxel
        plt.xlim(min_x - margin, max_x + margin)
        plt.ylim(min_y - margin, max_y + margin)
#        plt.xticks(np.linspace((min_x-margin)*spacing[0]+origin[0], (max_x+margin)*spacing[0]+origin[0], 6))

    plt.plot(isocenter_voxel[0], isocenter_voxel[1], "yx", label="Isocenter", markersize=10)

    # Disegna i rettangoli per ogni angolo
    for rect in rectangles:
        bl = (rect["bottom_left"]-origin[:2]) / spacing[:2]
        tr = (rect["top_right"]-origin[:2]) / spacing[:2]
        rect_x = [bl[0], tr[0], tr[0], bl[0], bl[0]]
        rect_y = [bl[1], bl[1], tr[1], tr[1], bl[1]]
        plt.plot(rect_x, rect_y, "g--", alpha=0.7, label=f"Angle {rect['angle']}\u00b0")

    # Disegna il miglior rettangolo
    best_bl = best_rectangle["bottom_left"] / spacing[:2]
    best_tr = best_rectangle["top_right"] / spacing[:2]
    best_rect_x = [best_bl[0], best_tr[0], best_tr[0], best_bl[0], best_bl[0]]
    best_rect_y = [best_bl[1], best_bl[1], best_tr[1], best_tr[1], best_bl[1]]
    plt.plot(best_rect_x, best_rect_y, "r-", linewidth=2, label=f"Best Rectangle (Angle {best_rectangle['angle']}\u00b0)")

    for angle in rotation_angles:
#        rotated_contour = rotate_contour(contour, angle, isocenter[:2] / spacing[:2])
        rotated_contour = rotate_contour(contour, angle, isocenter_voxel[:2])
        if angle == best_rectangle["angle"]:
            plt.plot(rotated_contour[:, 0], rotated_contour[:, 1], "r-", linewidth=2, label=f"Best Contour (Angle {angle}\u00b0)")
        else:
            plt.plot(rotated_contour[:, 0], rotated_contour[:, 1], "c--", alpha=0.5)

  
    plt.legend()
    plt.title("PTV Contour and Fitted Rectangles")
    plt.xlabel("X (voxels)")
    plt.ylabel("Y (voxels)")
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Find the best rectangle fitting a PTV projection.")
    parser.add_argument("ptv_path", help="Path to the PTV file.")
    parser.add_argument("-rot", nargs="+", type=int, default=[0], help="Rotation angles (degrees).")
    parser.add_argument("-plot", action="store_true", help="Enable visualization.")
    args = parser.parse_args()

    invalid_angles = [angle for angle in args.rot if angle < 0 or angle > 90]
    if invalid_angles:
        print(f"Error: Invalid rotation angles {invalid_angles}. Angles must be between 0 and 90 degrees.")
        sys.exit(1)  # Termina il programma con un codice di errore

    ptv_array, spacing, origin = load_ptv(args.ptv_path)
    print("spacing[cm]:", spacing)
    print("origin[cm]:", origin)
    print(ptv_array.shape)
    isocenter = calculate_isocenter(ptv_array, spacing, origin)
    iso_idx = (isocenter-origin)/spacing
    print(f"Isocenter [cm]: X={isocenter[0]:.2f}, Y={isocenter[1]:.2f}, Z={isocenter[2]:.2f}")
    print("Isocenter [index]:", iso_idx)
    projection, contour = find_contour_projection(ptv_array)
    if contour is None:
        print("No contour found for the PTV projection.")
        return
    contour=np.array(contour)[::-1]
    print(contour.shape)
    rectangles, best_rectangle = find_best_rectangle(projection, contour, spacing, origin, isocenter, iso_idx, args.rot)

    for R in rectangles:
        print(f"  Angle: {R['angle']}\u00b0")
        print(f"  Width: {R['width']:.2f} cm")
        print(f"  Height: {R['height']:.2f} cm")
        print(f"  Left: {R['relative_coords']['left']:.2f} mm")
        print(f"  Right: {R['relative_coords']['right']:.2f} mm")
        print(f"  Bottom: {R['relative_coords']['bottom']:.2f} mm")
        print(f"  Top: {R['relative_coords']['top']:.2f} mm")
        print(f"  Area Difference: {R['area_diff']:.2f} mm²")
        print("")
    print("***************")
    print("BEST RECTANGLE:")
    print(f"  Angle: {best_rectangle['angle']}\u00b0")
    print(f"  Width: {best_rectangle['width']:.2f} cm")
    print(f"  Height: {best_rectangle['height']:.2f} cm")
    print(f"  Left: {best_rectangle['relative_coords']['left']:.2f} mm")
    print(f"  Right: {best_rectangle['relative_coords']['right']:.2f} mm")
    print(f"  Bottom: {best_rectangle['relative_coords']['bottom']:.2f} mm")
    print(f"  Top: {best_rectangle['relative_coords']['top']:.2f} mm")
    print(f"  Area Difference: {best_rectangle['area_diff']:.2f} mm²")

    if args.plot:
        plot_results(projection, contour, rectangles, best_rectangle, isocenter, origin, spacing, args.rot)

if __name__ == "__main__":
    main()

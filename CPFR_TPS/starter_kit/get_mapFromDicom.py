import SimpleITK as sitk
import os
import argparse

def dicom_to_mhd(dicom_folder, output_filename):
    reader = sitk.ImageSeriesReader()
    dicom_files = reader.GetGDCMSeriesFileNames(dicom_folder)
    reader.SetFileNames(dicom_files)
    
    print(f"Found {len(dicom_files)} DICOM files in {dicom_folder}.")

    image = reader.Execute()
    
    sitk.WriteImage(image, output_filename)
    print(f"Converted DICOM to MHD: {output_filename}")

def main():
    parser = argparse.ArgumentParser(description="Convert a DICOM folder to MHD/RAW files")
    parser.add_argument("dicom_folder", help="Path to the folder containing DICOM images")
    parser.add_argument("output_filename", help="Path to the output MHD file (e.g., 'output.mhd')")
    args = parser.parse_args()
    
    dicom_to_mhd(args.dicom_folder, args.output_filename)

if __name__ == "__main__":
    main()

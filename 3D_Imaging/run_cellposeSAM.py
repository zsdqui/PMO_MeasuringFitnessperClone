# Runs Cellpose-SAM on 3D-tiff files (nuclear channel)
import os, glob, re
import numpy as np
import pandas as pd
from tifffile import imread
import argparse, traceback
from skimage.color import label2rgb

try:
    from cellpose import models
    from cellpose import io, utils
    import torch

    print("Cellpose library imported successfully.")
    CELLPOSE_AVAILABLE = True
except ImportError:
    print("🛑 Cellpose library not found. Please install it (pip install cellpose).")
    print("🛑 Prediction step will be skipped.")
    CELLPOSE_AVAILABLE = False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir",help= "Input directory for 3D Tiffs" )
    parser.add_argument("-c", "--colorize", help="colorize output", optional=True )
    parser.add_argument("-o", "--output", help="Output directory for 3D masks")
    args = parser.parse_args()

    image_base_dir = args.input_dir
    volume_subdirs = []
    try:
        if not os.path.isdir(image_base_dir):
             raise FileNotFoundError(f"Image base directory not found: {image_base_dir}")
        volume_subdirs = sorted([
            d for d in os.listdir(image_base_dir)
            if os.path.isdir(os.path.join(image_base_dir, d))
        ])
        print(f"Found {len(volume_subdirs)} potential volume directories in {image_base_dir}.")
        if not volume_subdirs:
            print("Warning: No subdirectories found. Check the image_base_dir path and contents.")

    except FileNotFoundError as e:
        print(f"ERROR: {e}")
    except Exception as e:
        print(f"ERROR: An unexpected error occurred while listing directories: {e}")
    if not volume_subdirs:  # Check if variable from Cell 1 exists and is not empty
        print(
            "No volume subdirectories found (list 'volume_subdirs' is empty or was not created in Cell 1). Exiting batch processing."
        )
    # Also check if Cellpose is ready if we intend to use it
    elif not CELLPOSE_AVAILABLE or cellpose_model is None:
        print(
            "Cellpose is not available or model failed to initialize. Cannot perform predictions. Exiting batch processing."
        )
    else:
        processed_volume_count = 0
        # Set how many volumes to process (e.g., 5 for testing, None for all)
        max_volumes_to_process = 1  # Set to None to process all found volumes

        print(
            f"\nStarting batch processing loop for up to {max_volumes_to_process if max_volumes_to_process is not None else len(volume_subdirs)} volumes..."
        )
        print(
            f"Using Cellpose Model: {CELLPOSE_MODEL_TYPE}, Diameter: {CELLPOSE_DIAMETER}, Anisotropy: {CELLPOSE_ANISOTROPY}, Device: {CELLPOSE_DEVICE}"
        )

        for volume_index, volume_name in enumerate(volume_subdirs):
            # Check processing limit
            if (
                max_volumes_to_process is not None
                and processed_volume_count >= max_volumes_to_process
            ):
                print(
                    f"\nReached processing limit ({max_volumes_to_process}). Stopping."
                )
                break

            print(
                f"\n--- Processing Volume {processed_volume_count + 1}/{max_volumes_to_process or len(volume_subdirs)}: {volume_name} ---"
            )
            volume_img_dir_path = os.path.join(image_base_dir, volume_name)

            # Reset variables for this volume
            image = None
            img_h, img_w, img_d = None, None, None  # Initialize dimensions
            start_time_volume_total = (
                time.time()
            )  # Start timer for the whole volume processing

            try:
                # STEP 1: Assemble Image 
                print("(Step 1) Assembling Image...")
                start_time_img = time.time()
                # for 2d files- image, img_h, img_w, img_d = assemble_3d_image(volume_img_dir_path)
                image = imread(
                    os.path.join(sample_img_dir_path, "nucleus.p.tif")
                ).transpose(1, 2, 0)
                [img_h, img_w, img_d] = np.shape(sample_image)
                if image is None:
                    print(f"  Skipping volume {volume_name}: Image assembly failed.")
                    # No need to increment processed_volume_count here, finally block handles it
                    continue  # Skip to the next volume in the finally block
                img_assembly_time = time.time() - start_time_img
                print(
                    f"  Image assembly took {img_assembly_time:.2f} seconds. Image Shape (H,W,D): ({img_h}, {img_w}, {img_d})"
                )
            except Exception as e:
                pass

            try:
                print("  Waiting for Cellpose inference result...")
                masks, _, _ = model.eval(
                    image,
                    do_3d=True,
                    batch_size=16,
                    z_axis=2,
                    channel_axis=None,
                    flow3D_smooth=1,
                )

                print("  Cellpose inference finished.")
                if masks is None:
                    pred_error = "Cellpose inference function returned None"
            except Exception as e_pred:
                print(
                    f"!! Error during threaded Cellpose inference for {volume_name}: {e_pred}"
                )
                traceback.print_exc()
                # Handle failure, maybe set masks_volume to None
                masks_volume = None  # Indicate prediction failure
                pred_error = str(e_pred)
                print(f"  (Cellpose inference encountered issue: {pred_error})")
            print(f" Inference took {time.time() - start_time_cellpose:.2f} seconds.")

            if args.colorize:
                masks = label2rgb(masks, bg_label=0)

        output_path = 
        tifffile.imwrite(
            output_path,
            (masks * 255).astype(np.uint8),
            imagej=True
        )
    print(f"Colorized image saved to: {output_path}")


if __name__ == "__main__":
    main()

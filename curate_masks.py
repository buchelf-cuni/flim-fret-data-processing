"""Manual curation of segmentation masks in napari.

Run `python curate_masks.py` and select, in this order:

1. the folder with the intensity images (`<sample>_export`),
2. one or more folders with label images (`<sample>_ib_masks`,
   `<sample>_cell_masks`); press Cancel when all are selected,
3. the output folder (the curated masks go to `<output>/<label folder name>/`).

Each image opens with its label layers, which are edited with napari's label
tools. A label image belongs to an image when its name starts with the image
name followed by '_' (img.tif -> img_ib_labels.tif).

Keys:  s  save the labels and go to the next image
       n  next image without saving
       p  previous image
"""
import os

import napari
import numpy as np
from qtpy.QtWidgets import QApplication, QFileDialog
from skimage.io import imread, imsave


current = 0          # index of the image being shown
pairs = []           # (image path, [label paths])
output_folder = None

def select_directory(dialogTitle, start_dir=None, required=True):
    """Ask for a folder; returns None if cancelled and not `required`."""
    QApplication.instance() or QApplication([])
    directory_path = QFileDialog.getExistingDirectory(
        None, dialogTitle, start_dir or "", QFileDialog.ShowDirsOnly
    )
    if not directory_path:
        if required:
            raise ValueError("No directory selected.")
        return None
    return directory_path


def select_label_directories(start_dir=None):
    """Ask for one label folder per mask type (e.g. _ib_masks and _cell_masks)."""
    folders = []
    while True:
        title = ("Select folder containing labels." if not folders
                 else "Select another label folder (Cancel when done).")
        folder = select_directory(title, start_dir, required=not folders)
        if folder is None:
            break
        if folder not in folders:
            folders.append(folder)
    return folders


def get_image_label_pairs(image_fldr, label_fldrs):
    """[(image path, [label paths])] for every image that has at least one label image."""
    images = sorted(f for f in os.listdir(image_fldr) if f.lower().endswith(".tif"))
    labels_by_folder = [
        (fldr, sorted(f for f in os.listdir(fldr) if f.lower().endswith(".tif")))
        for fldr in label_fldrs
    ]

    image_label_pairs = []

    for image in images:
        image_stem = os.path.splitext(image)[0]
        image_path = os.path.join(image_fldr, image)

        label_paths = []
        for label_fldr, labels in labels_by_folder:
            for label in labels:
                if label.startswith(image_stem + "_"):
                    label_paths.append(os.path.join(label_fldr, label))

        if label_paths:
            image_label_pairs.append((image_path, label_paths))
    return image_label_pairs


def report_no_pairs(image_fldr, label_fldrs):
    """Explain why nothing matched instead of just saying that nothing did."""
    print("No matching image-label pairs found.")

    def tifs(folder):
        return sorted(f for f in os.listdir(folder) if f.lower().endswith(".tif"))

    images = tifs(image_fldr)
    print(f"  {len(images)} .tif in image folder {image_fldr}")
    for name in images[:3]:
        print(f"    {name}")
    if not images:
        subdirs = sorted(d for d in os.listdir(image_fldr)
                         if os.path.isdir(os.path.join(image_fldr, d)))
        if subdirs:
            print("  (this folder only holds subfolders; the script does not recurse "
                  f"- pick one of: {', '.join(subdirs[:5])}{' ...' if len(subdirs) > 5 else ''})")

    for label_fldr in label_fldrs:
        labels = tifs(label_fldr)
        print(f"  {len(labels)} .tif in label folder {label_fldr}")
        for name in labels[:3]:
            print(f"    {name}")

    print("  A label matches when its name starts with the image name (minus .tif) plus '_',")
    print("  e.g. img.tif -> img_ib_labels.tif. An identically named label never matches.")


def label_layer_names(image_stem, label_paths):
    """Name each labels layer after the part of its filename following the image stem."""
    names = []
    for label_path in label_paths:
        stem = os.path.splitext(os.path.basename(label_path))[0]
        name = stem[len(image_stem) + 1:] or stem
        if name in names:
            name = f"{name}_{len(names)}"
        names.append(name)
    return names


def bind_navigation_keys(viewer, layer):
    layer.bind_key("s", lambda _layer: save_and_next_image(viewer), overwrite=True)
    layer.bind_key("n", lambda _layer: next_image(viewer), overwrite=True)
    layer.bind_key("p", lambda _layer: previous_image(viewer), overwrite=True)


def load_pair(viewer, i):
    viewer.layers.clear()

    image_path, label_paths = pairs[i]
    image_stem = os.path.splitext(os.path.basename(image_path))[0]
    image = imread(image_path)

    viewer.add_image(image, name=os.path.basename(image_path))

    names = label_layer_names(image_stem, label_paths)
    for label_path, name in zip(label_paths, names):
        labels_layer = viewer.add_labels(imread(label_path), name=name)
        bind_navigation_keys(viewer, labels_layer)

    print(f"Loaded {i + 1}/{len(pairs)}: {os.path.basename(image_path)} + {', '.join(names)}")


def save_current_mask(viewer):
    image_path, label_paths = pairs[current]
    image_stem = os.path.splitext(os.path.basename(image_path))[0]

    for label_path, name in zip(label_paths, label_layer_names(image_stem, label_paths)):
        if name not in viewer.layers:
            print(f"Layer '{name}' is gone, not saving it.")
            continue
        # mirror the source folder so the mask types stay separated in the output
        save_dir = os.path.join(output_folder, os.path.basename(os.path.dirname(label_path)))
        os.makedirs(save_dir, exist_ok=True)
        save_path = os.path.join(save_dir, os.path.basename(label_path))
        imsave(save_path, viewer.layers[name].data.astype(np.uint16), check_contrast=False)
        print(f"Saved {save_path}")


def advance_image(viewer, step):
    global current

    if not pairs:
        print("No matching image-label pairs found.")
        return

    new_current = current + step

    if 0 <= new_current < len(pairs):
        current = new_current
        load_pair(viewer, current)
    elif new_current >= len(pairs):
        current = len(pairs)
        viewer.layers.clear()
        print("Finished.")
    else:
        print("Already at the first pair.")


def save_and_next_image(viewer):
    if not pairs:
        print("No matching image-label pairs found.")
        return
    if current >= len(pairs):
        print("Finished.")
        return

    save_current_mask(viewer)
    advance_image(viewer, 1)


def next_image(viewer):
    advance_image(viewer, 1)


def previous_image(viewer):
    advance_image(viewer, -1)


def main():
    global current, pairs, output_folder

    images_folder = select_directory("Select folder containing images.")
    # start the later dialogs beside the image folder, not inside it
    sibling_dir = os.path.dirname(images_folder)
    labels_folders = select_label_directories(sibling_dir)
    output_folder = select_directory("Select output folder.", sibling_dir)

    print(f"Images: {images_folder}")
    for labels_folder in labels_folders:
        print(f"Labels: {labels_folder}")
    print(f"Output: {output_folder}")

    pairs = get_image_label_pairs(images_folder, labels_folders)
    current = 0

    if not pairs:
        report_no_pairs(images_folder, labels_folders)
        return

    viewer = napari.Viewer()

    viewer.bind_key("s", save_and_next_image)
    viewer.bind_key("n", next_image)
    viewer.bind_key("p", previous_image)

    load_pair(viewer, current)

    napari.run()


if __name__ == "__main__":
    main()

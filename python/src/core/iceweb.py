#!/usr/bin/env python
#
# iceweb: create a standalone HTML webpage of an image
# folder containing either pdf or png images in (sub)subfolders
# 
# Usage:
#   python -m core.iceweb <folder_path>
#   then a standalone 'html' folder is created
# 
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


import argparse
import html
import os
import pathlib
import re
import shutil
import time
from datetime import datetime
from urllib.parse import quote  # For URL encoding

from pdf2image import convert_from_path
from termcolor import cprint

from core import resource
from core.io.files import ensure_dir

HTML_TEMPLATE = resource("templates/iceweb.html")


# Format the iceweb completion message with elapsed wall-clock time
def format_done_message(elapsed_seconds):
    return f"[iceweb: done in {elapsed_seconds:.1f} sec]"


def sanitize_filename(filename):
    """
    Function to sanitize filenames and directory names (remove or replace problematic characters)
    """
    return re.sub(
        r"[^\w\-_\. ]", "_", filename
    )  # Replace anything that's not alphanumeric, _, -, ., or space with _


def convert_pdf_to_png(pdf_path, output_folder):
    """
    Function to convert PDF to PNG and save it in the output folder
    """
    print(f"Converting PDF to PNG: {pdf_path}")

    # Get the original modification timestamp of the PDF file
    modification_time = os.path.getmtime(pdf_path)
    formatted_modification_time = time.strftime(
        "%Y-%m-%d %H:%M:%S", time.localtime(modification_time)
    )

    images = convert_from_path(pdf_path)
    png_paths = []
    stem = sanitize_filename(os.path.splitext(os.path.basename(pdf_path))[0])
    ensure_dir(output_folder, exist_ok=True)  # Ensure the directory is created
    for page, image in enumerate(images, start=1):
        suffix = "" if len(images) == 1 else f"__page_{page:03d}"
        output_file = os.path.join(output_folder, f"{stem}{suffix}.png")
        image.save(output_file, "PNG")
        png_paths.append(output_file)

    return png_paths, formatted_modification_time


def copy_png_to_output(png_path, output_subfolder):
    """
    Function to copy PNG files to the output folder, preserving folder structure
    """
    sanitized_filename = sanitize_filename(os.path.basename(png_path))
    ensure_dir(output_subfolder, exist_ok=True)

    output_file = os.path.join(output_subfolder, sanitized_filename)
    shutil.copy2(png_path, output_file)  # Copy with metadata preservation

    # Get the original modification timestamp of the PNG file
    modification_time = os.path.getmtime(png_path)
    formatted_modification_time = time.strftime(
        "%Y-%m-%d %H:%M:%S", time.localtime(modification_time)
    )

    return output_file, formatted_modification_time


# Create one HTML gallery from the collected image hierarchy
def create_html(png_files, output_file, image_size, base_folder_name):
    """
    Function to create HTML page with hierarchical folder structure
    """
    generation_timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    print(f"Generating HTML file: {output_file}")
    content = []

    current_root = None
    current_subfolder = None
    current_deeper_subfolder = None
    inside_grid = False

    for folder_name, files in png_files.items():
        folder_hierarchy = folder_name.split(os.sep)
        root_folder = folder_hierarchy[0]

        if root_folder == ".":
            root_folder = base_folder_name

        if root_folder != current_root:
            if inside_grid:
                content.append("</div>")
            if current_root is not None:
                content.append("<hr>")
            content.append(f"<h2>{html.escape(root_folder)}</h2><div class='grid-container'>")
            inside_grid = True
            current_root = root_folder
            current_subfolder = None

        if len(folder_hierarchy) > 1:
            subfolder = folder_hierarchy[1]
            if subfolder != current_subfolder:
                if current_subfolder is not None:
                    content.append("</div><div class='grid-container'>")
                content.append(f"<h3>{html.escape(subfolder)}</h3>")
                current_subfolder = subfolder
                current_deeper_subfolder = None

        if len(folder_hierarchy) > 2:
            deeper_subfolder = "/".join(folder_hierarchy[2:])
            if deeper_subfolder != current_deeper_subfolder:
                if current_deeper_subfolder is not None:
                    content.append("</div><div class='grid-container'>")
                content.append(f"<h4>{html.escape(deeper_subfolder)}</h4>")
                current_deeper_subfolder = deeper_subfolder

        for file_name, png_file, original_timestamp in files:
            relative_path = os.path.relpath(png_file, os.path.dirname(output_file))
            relative_path = quote(relative_path.replace(os.sep, "/"))

            # Add image with the original file modification timestamp
            title = html.escape(file_name)
            content.append(
                f"""
            <div class="grid-item">
                <div class="title">{title}</div>
                <a href="{relative_path}" target="_blank">
                    <img src="{relative_path}" alt="{title}">
                </a>
                <div class="timestamp">Last modified: {original_timestamp}</div>
            </div>
            """
            )

    if inside_grid:
        content.append("</div>")

    document = HTML_TEMPLATE.read_text(encoding="utf-8")
    document = document.replace("__ICEWEB_IMAGE_SIZE__", str(int(image_size)))
    document = document.replace("__ICEWEB_GENERATED_AT__", generation_timestamp)
    document = document.replace("__ICEWEB_CONTENT__", "\n".join(content))
    ensure_dir(pathlib.Path(output_file).parent)
    pathlib.Path(output_file).write_text(document, encoding="utf-8")
    cprint(f"HTML file generated: {output_file}", "green")


def process_folders(base_folder, output_folder, html_output, image_size):
    """
    Main function to process folders recursively
    """
    png_files = {}

    ensure_dir(output_folder)

    base_folder_name = os.path.basename(base_folder)

    # Walk through all subfolders and process PDFs and PNGs, skip the "html" output folder
    for root, dirs, files in os.walk(base_folder):
        dirs[:] = [d for d in dirs if d != "html"]

        relative_folder_path = os.path.relpath(root, base_folder)
        output_subfolder = os.path.join(output_folder, relative_folder_path)
        png_files[relative_folder_path] = []

        cprint(f"Processing folder: {relative_folder_path}", "yellow")

        png_file_basenames = {os.path.splitext(f)[0] for f in files if f.endswith(".png")}

        for file in files:
            if file.endswith(".pdf"):
                if os.path.splitext(file)[0] in png_file_basenames:
                    continue

                pdf_path = os.path.join(root, file)
                png_paths, pdf_timestamp = convert_pdf_to_png(pdf_path, output_subfolder)

                for png_file in png_paths:
                    png_files[relative_folder_path].append(
                        (os.path.splitext(os.path.basename(png_file))[0], png_file, pdf_timestamp)
                    )

            elif file.endswith(".png"):
                png_path = os.path.join(root, file)
                copied_png_path, png_timestamp = copy_png_to_output(png_path, output_subfolder)
                png_files[relative_folder_path].append(
                    (os.path.splitext(file)[0], copied_png_path, png_timestamp)
                )

    create_html(png_files, html_output, image_size, base_folder_name)


# Run the iceweb conversion and HTML generation workflow
def main():
    start_time = time.perf_counter()

    # Use argparse to parse command-line arguments
    parser = argparse.ArgumentParser(description="Convert PDFs to PNGs and generate an HTML grid.")
    parser.add_argument(
        "base_folder", type=str, help="Base folder containing PDF and PNG files and subfolders."
    )
    parser.add_argument(
        "--image_size",
        type=int,
        default=400,
        help="Size of the images to display (default: 400px).",
    )

    args = parser.parse_args()

    base_folder = args.base_folder  # Base folder provided as argument
    image_size = args.image_size  # Image size specified as argument
    output_folder = os.path.join(base_folder, "html")  # Folder to save PNGs
    html_output = os.path.join(output_folder, "main.html")  # HTML output file

    cprint("Starting PDF to PNG conversion and HTML generation...", "green")
    process_folders(base_folder, output_folder, html_output, image_size)
    elapsed_seconds = time.perf_counter() - start_time
    cprint(format_done_message(elapsed_seconds), "green")


if __name__ == "__main__":
    main()

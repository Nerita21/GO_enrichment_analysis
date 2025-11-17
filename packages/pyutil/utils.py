import requests 
import os
import gzip
import shutil
from pathlib import Path

def download_mesh():
    url = "https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/desc2025.xml"

    filename = os.path.basename(url)

    response = requests.get(url)

    if response.status_code == 200:
        with open(filename, 'wb') as file:
            file.write(response.content)
        print(f"File downloaded and saved as {filename}")
    else:
        print(f"Failed to download file. Status code: {response.status_code}")

def unzip(file_path_input: str, file_path_output: str | None = None) -> None:
    """Unzips a .gz file and writes the decompressed file to the specified output path."""
    
    input_path = Path(file_path_input)

    if input_path.suffix != ".gz":
        raise ValueError(f"Input file must have a .gz extension, got {input_path.suffix}")

    # If no output path is provided, remove the .gz suffix and save in the same location
    output_path = Path(file_path_output) if file_path_output else input_path.with_suffix("")

    # Unzip the gzipped file into a plain file
    with gzip.open(input_path, "rb") as f_in:
        with open(output_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    print(f"Unzipped file written to: {output_path}")
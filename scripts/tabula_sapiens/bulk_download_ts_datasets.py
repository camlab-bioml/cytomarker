import requests
import os
from tqdm import *
import zipfile

tissue_links = {'Bladder': 34701958, 'Blood': 34701964, 'Bone_Marrow': 34701967,
                'Eye': 34701970, 'Fat': 34701973, 'Heart': 34701976, 'Kidney': 34701979,
                'Large_Intestine': 34701982, 'Liver': 34701985,
                'Lung': 34701991, 'Lymph_Node': 34701994,
                'Mammary': 34701997, 'Muscle': 34702000, 'Pancreas': 34702003,
                'Prostate': 34702006, 'Salivary_Gland': 34702009, 'Skin': 34702012,
                'Small_Intestine': 34702015, 'Spleen': 34702018, 'Thymus': 34702027,
                }

dest_dir = "/home/matt/cytosel/tabula_sapiens/datasets/"

for key, value in tissue_links.items():
    print(key, value)
    url = "https://figshare.com/ndownloader/files/" + str(value)
    print(url)
    content = requests.get(url, stream=True)
    dest_path = os.path.join(dest_dir, key + ".h5ad.zip")

    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        if not os.path.isfile(dest_path):
            with open(dest_path, 'wb') as f:
                pbar = tqdm(total=int(r.headers['Content-Length']))
                for chunk in r.iter_content(chunk_size=8192):
                    if chunk:  # filter out keep-alive new chunks
                        f.write(chunk)
                        pbar.update(len(chunk))

    with zipfile.ZipFile(dest_path, 'r') as zip_ref:
        zip_ref.extractall(dest_dir)

#!/usr/bin/env python

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
                'Tongue': 34702030, 'Trachea': 34702033, 'Uterus': 34702036, 'Vasculature': 34702039,
                'Endothelial': 34702042, 'Epithelial': 34702048, 'Germ_Line': 34702051, 'Stromal': 34702054,
                'Immune': 34702069
                }

dest_dir = "/home/matt/cytosel/tabula_sapiens/datasets/"

for key, value in tissue_links.items():
    url = "https://figshare.com/ndownloader/files/" + str(value)
    dest_path = os.path.join(dest_dir, key + ".h5ad.zip")
    dest_file = os.path.join(dest_dir, "TS_" + key + ".h5ad")

    # check if the uncompressed version is present before downloading
    # uncompressed files have the format TS_{tissue}.h5ad such as TS_Liver.h5ad

    if not os.path.isfile(dest_file):
        print(key, value)
        print(url)

        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(dest_path, 'wb') as f:
                pbar = tqdm(total=int(r.headers['Content-Length']))
                for chunk in r.iter_content(chunk_size=8192):
                    if chunk:  # filter out keep-alive new chunks
                        f.write(chunk)
                        pbar.update(len(chunk))

            with zipfile.ZipFile(dest_path, 'r') as zip_ref:
                zip_ref.extractall(dest_dir)

import numpy as np
import Bio
import math
import pandas as pd
import requests
import io
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.ticker import MultipleLocator, FixedLocator
from PIL import Image,ImageFont
from pathlib import Path
from helpers import get_integers_around
from mpl_toolkits.mplot3d import Axes3D
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

def get_mmcif(uniprot):
    response = requests.get(f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot}").json()
    mmcif_url = response[0]['cifUrl']

    response = requests.get(mmcif_url)

    return response.text

def mmcif_to_df(filename):
    sonic_dict = MMCIF2Dict(filename)
    x_coords = [float(i)for i in sonic_dict['_atom_site.Cartn_x']]
    y_coords = [float(i)for i in sonic_dict['_atom_site.Cartn_y']]
    z_coords = [float(i)for i in sonic_dict['_atom_site.Cartn_z']]

    x_coords = [i - min(x_coords) for i in x_coords]
    y_coords = [i - min(y_coords) for i in y_coords]
    z_coords = [i - min(z_coords) for i in z_coords]


    aminos = sonic_dict['_atom_site.label_comp_id']

    df = pd.DataFrame({
        'x' : x_coords,
        'y': y_coords,
        'z' : z_coords,
        'amino_acid' : aminos
    })
    return df

def create_array(coords_df, size=20, buffersize=1):
    maxx = coords_df['x'].max()
    maxy = coords_df['y'].max()
    maxz = coords_df['z'].max()

    # Get scaling factor for everything
    maxbottom = max(maxx,maxy)
    factor = size/maxbottom

    coords_df['x'] = coords_df['x']*factor
    coords_df['y'] = coords_df['y']*factor
    #divide by two for the z coordinate, 
    # because lego brings are 1.2x taller than they are wide
    coords_df['z'] = (coords_df['z']*factor) / 1.2

    # get Height
    height = math.ceil(coords_df['z'].max())

    print(coords_df.head())

    array_3d = np.zeros((size,size,height))

    print(array_3d.shape)

    for i in range(coords_df.shape[0]):
        x = list(coords_df['x'])[i]
        y = list(coords_df['y'])[i]
        z = list(coords_df['z'])[i]

        x_ints = get_integers_around(x,buffersize)
        y_ints = get_integers_around(y,buffersize)
        z_ints = get_integers_around(z,buffersize)

        for i in x_ints:
            for j in y_ints:
                for k in z_ints:
                    if i in range(0,size) and j in range(0,size) and k in range(0,height):
                        array_3d[i][j][k] = 1
                        
    return array_3d

def visualize_array(array):
    data = array
    axes=list(array.shape)
    alpha = 0.9

    colors = np.empty(axes + [4], dtype=np.float32)
    colors[:] = [1, 0, 0, alpha]  # red

    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.voxels(data, facecolors=colors, edgecolors='grey')
    plt.show()

def make_slices(array,uniprot):
    folder_path = Path(f'{uniprot}/slices')
    folder_path.mkdir(parents=True,exist_ok=True)

    for i in range(array.shape[2]):
        image = array[:,:,i]
        cmap = ListedColormap(['red','white','green'])
        norm = plt.Normalize(vmin=-1.5, vmax=1.5)
        plt.matshow(image,cmap=cmap, norm=norm)
        plt.xticks(np.arange(-0.5, image.shape[1], 1), labels=np.arange(0, image.shape[1]+1, 1))
        plt.yticks(np.arange(-0.5, image.shape[0], 1), labels=np.arange(0, image.shape[0]+1, 1))
        plt.grid(which='major', color='black', linestyle='-', linewidth=1)
        plt.grid(which='minor', color='gray', linestyle='--', linewidth=0.5)

        
        plt.savefig(f'{uniprot}/slices/{i}.png')
        plt.clf()

def mmcif_to_visualize(filename,size=20,buffersize=1):
    df = mmcif_to_df(filename)
    array = create_array(df,size,buffersize)
    visualize_array(array)

def uniprot_to_visualize(uniprot,size=20,buffersize=1):
    mmcif = get_mmcif(uniprot)

    df = mmcif_to_df(io.StringIO(mmcif))
    array = create_array(df,size,buffersize)
    visualize_array(array)

def array_to_overhang(array):
    overhang_array = np.copy(array)

    for x in range(array.shape[0]):
        for y in range(array.shape[1]):
            for z in range(array.shape[2]):
                if z > 0 and z < array.shape[2] - 1:
                    if overhang_array[x][y][z] == 1 and array[x][y][z-1] == 0 and array[x][y][z+1] == 0:
                        overhang_array[x][y][z] *= -1

                if z == array.shape[2] - 1:
                    if overhang_array[x][y][z] == 1 and array[x][y][z-1] == 0:
                        overhang_array[x][y][z] *= -1
    return overhang_array



if __name__ == "__main__":
    human_myosin = "P12882"
    # iga = "P01876"
    # mmcif = get_mmcif(iga)

    # df = mmcif_to_df(io.StringIO(mmcif))
    # array = create_array(df)
    # overhangarray = array_to_overhang(array)

    # make_slices(overhangarray,iga)

    uniprot_to_visualize(human_myosin)

    #resp = mmcif_to_visualize("pdb_files/1iga_updated.cif",size=100,buffersize=2)


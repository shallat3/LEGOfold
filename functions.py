import numpy as np
import Bio
import math
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3d
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

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
    # because lego brings are 2x taller than they are wide
    coords_df['z'] = (coords_df['z']*factor) / 2

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



def get_integers_around(x,buffer):
    lower = x - buffer/2
    upper = x + buffer/2

    return list(range(math.ceil(lower),math.ceil(upper)))



if __name__ == "__main__":
    df = mmcif_to_df('pdb_files/6pjv.cif')
    array = create_array(df)
    print(array.sum())

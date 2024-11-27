import numpy as np
import requests
import numpy as np
import Bio
import math
import pandas as pd
import requests
import io
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.colors import ListedColormap
from matplotlib.ticker import MultipleLocator, FixedLocator
from helpers import get_integers_around
from PIL import Image,ImageFont
from pathlib import Path
from mpl_toolkits.mplot3d import Axes3D
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from LegoPiece import LegoPiece

class LegoProtein:
    def __init__(self, size=20,uniprot = None,filename=None, buffersize=1):
        self.size = size
        self.uniprot = uniprot
        self.buffersize=buffersize
        self.piece_dimension_list = [(12,1),(1,12),(8,1),(1,8),(2,6),(6,2),(2,4),(4,2),(6,1),(1,6),(4,1),(1,4),(3,2),(2,3),(2,2),(1,3),(3,1),(2,1),(1,2),(1,1)]
        if filename:
            self.filename = filename
        elif uniprot:
            self.filename = self.__get_mmcif()

        
        self.pieces = []

        
        self.coords_df = self.__mmcif_to_df()
        self.array = self.__create_array()

        self.x_dim = self.array.shape[0]
        self.y_dim = self.array.shape[1]
        self.z_dim = self.array.shape[2]

        self.overhangarray = self.__array_to_overhang()

        self.piecearray = self.__placepieces()


    def __placepieces(self):
        filled = self.array.copy()

        for piece in self.piece_dimension_list:

            for z in range(self.z_dim):
                for x in range(self.x_dim):
                    for y in range(self.y_dim):
                        if self.__can_place(filled,piece,x,y,z):
                            legopiece = LegoPiece(piece,(x,y,z))
                            filled[x:x+piece[0],y:y+piece[1],z] = legopiece
                            self.__add_connections(filled,legopiece)
                            self.pieces.append(legopiece)

        assert((filled != 1).all())
        return filled

    def __add_connections(self,filled,legopiece):
        z = legopiece.get_coordinates()[2]
        

        if z != self.z_dim - 1:
            for x in range(legopiece.get_coordinates()[0],legopiece.get_coordinates()[0]+legopiece.get_dimensions()[0]):
                for y in range(legopiece.get_coordinates()[1],legopiece.get_coordinates()[1]+legopiece.get_dimensions()[1]):
                    
                    above = filled[x,y,z+1]
                    
                    if  isinstance(above,LegoPiece):
                        above.add_connection(legopiece)
                        legopiece.add_connection(above)
        if z != 0:
            for x in range(legopiece.get_coordinates()[0],legopiece.get_coordinates()[0]+legopiece.get_dimensions()[0]):
                for y in range(legopiece.get_coordinates()[1],legopiece.get_coordinates()[1]+legopiece.get_dimensions()[1]):
                    below = filled[x,y,z-1]
                    if isinstance(below,LegoPiece):
                        below.add_connection(legopiece)
                        legopiece.add_connection(below)
    
    def __can_place(self,filled,dimensions,x,y,z):
        try:
            if (filled[x:(x+dimensions[0]),y:(y+dimensions[1]),z] == 1).all() and x+dimensions[0] - 1<filled.shape[0] and y+dimensions[1] - 1 < filled.shape[1]:
                return True
            else:
                return False
        except:
            return False

    def __get_mmcif(self):
        response = requests.get(f"https://alphafold.ebi.ac.uk/api/prediction/{self.uniprot}").json()
        mmcif_url = response[0]['cifUrl']

        response = requests.get(mmcif_url)

        return io.StringIO(response.text)

    def __mmcif_to_df(self):
        sonic_dict = MMCIF2Dict(self.filename)
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

    def __create_array(self):
        maxx = self.coords_df['x'].max()
        maxy = self.coords_df['y'].max()
        maxz = self.coords_df['z'].max()

        # Get scaling factor for everything
        maxbottom = max(maxx,maxy)
        factor = self.size/maxbottom

        self.coords_df['x'] = self.coords_df['x']*factor
        self.coords_df['y'] = self.coords_df['y']*factor
        #divide by 1.2 for the z coordinate, 
        # because lego brings are 1.2x taller than they are wide
        self.coords_df['z'] = (self.coords_df['z']*factor) / 1.2

        # get Height
        height = math.ceil(self.coords_df['z'].max())

        print(self.coords_df.head())

        array_3d = np.zeros((self.size,self.size,height),dtype=object)

        print(array_3d.shape)

        for i in range(self.coords_df.shape[0]):
            x = list(self.coords_df['x'])[i]
            y = list(self.coords_df['y'])[i]
            z = list(self.coords_df['z'])[i]

            x_ints = get_integers_around(x,self.buffersize)
            y_ints = get_integers_around(y,self.buffersize)
            z_ints = get_integers_around(z,self.buffersize)

            for i in x_ints:
                for j in y_ints:
                    for k in z_ints:
                        if i in range(0,self.size) and j in range(0,self.size) and k in range(0,height):
                            array_3d[i][j][k] = 1
                            
        return array_3d

    def visualize_array(self):
        data = self.array
        axes=list(self.array.shape)
        alpha = 0.9

        colors = np.empty(axes + [4], dtype=np.float32)
        colors[:] = [1, 0, 0, alpha]  # red

        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')
        ax.voxels(data, facecolors=colors, edgecolors='grey')
        plt.show()

    def make_slices(self):
        folder_path = Path(f'{self.uniprot}/slices')
        folder_path.mkdir(parents=True,exist_ok=True)

        for i in range(self.array.shape[2]):
            image = self.array[:,:,i]
            cmap = ListedColormap(['red','white','green'])
            norm = plt.Normalize(vmin=-1.5, vmax=1.5)
            plt.matshow(image,cmap=cmap, norm=norm)
            plt.xticks(np.arange(-0.5, image.shape[1], 1), labels=np.arange(0, image.shape[1]+1, 1))
            plt.yticks(np.arange(-0.5, image.shape[0], 1), labels=np.arange(0, image.shape[0]+1, 1))
            plt.grid(which='major', color='black', linestyle='-', linewidth=1)
            plt.grid(which='minor', color='gray', linestyle='--', linewidth=0.5)

            
            plt.savefig(f'{self.uniprot}/slices/{i}.png')
            plt.clf()

    def connection_graph(self):
        lego_graph = nx.Graph()

        for piece in self.pieces:
            lego_graph.add_node(piece.piece_id())

        for piece in self.pieces:
            for conn in piece.connections():
                
                lego_graph.add_edge(piece.piece_id(), conn.piece_id())

        self.graph = lego_graph
        return lego_graph

    def __array_to_overhang(self):
        overhang_array = np.copy(self.array)

        for x in range(self.array.shape[0]):
            for y in range(self.array.shape[1]):
                for z in range(self.array.shape[2]):
                    if z > 0 and z < self.array.shape[2] - 1:
                        if overhang_array[x][y][z] == 1 and self.array[x][y][z-1] == 0 and self.array[x][y][z+1] == 0:
                            overhang_array[x][y][z] *= -1

                    if z == self.array.shape[2] - 1:
                        if overhang_array[x][y][z] == 1 and self.array[x][y][z-1] == 0:
                            overhang_array[x][y][z] *= -1
        return overhang_array
    
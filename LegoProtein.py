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
from matplotlib.colors import ListedColormap,LinearSegmentedColormap, TwoSlopeNorm,BoundaryNorm,Normalize
from matplotlib.ticker import MultipleLocator, FixedLocator
from helpers import get_integers_around
from PIL import Image,ImageFont
from pathlib import Path
import shutil
from mpl_toolkits.mplot3d import Axes3D
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from LegoPiece import LegoPiece

class LegoProtein:
    def __init__(self, size=20,uniprot = None,filename=None, buffersize=1):
        self.size = size
        self.uniprot = uniprot
        self.buffersize=buffersize
        self.piece_dimensions_overhang = [(1,2),(2,1),(1,3),(3,1),(4,1),(1,4),(2,2)]
        self.piece_dimension_list = [(12,1),(1,12),(8,1),(1,8),(2,6),(6,2),(2,4),(4,2),(1,6),(6,1),(4,1),(1,4),(3,2),(2,3),(2,2),(1,3),(3,1),(2,1),(1,2),(1,1)]
        self.piece_dimension_list_odd = [(1,12),(12,1),(1,8),(8,1),(2,6),(6,2),(2,4),(4,2),(6,1),(1,6),(1,4),(4,1),(2,3),(3,2),(2,2),(3,1),(1,3),(1,2),(2,1),(1,1)]
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
        self.pieces_visual_array = self.overhangarray.copy()
        self.pieces_visual_no_overhang = self.overhangarray.copy()

        self.piecearray = self.__placepieces()


    def __placepieces(self):
        filled = self.overhangarray.copy().astype(object)
        print(np.sum(filled==-1))
        for piece in self.piece_dimensions_overhang:
            for z in range(self.z_dim):
                for x in range(self.x_dim):
                    for y in range(self.y_dim):
                        if self.__can_place(filled,piece,x,y,z) and (filled[x:x+piece[0],y:y+piece[1],z] == -1).any():
                            legopiece = LegoPiece(piece,(x,y,z))
                            self.pieces_visual_array[x:x+piece[0],y:y+piece[1],z] = legopiece.color
                            filled[x:x+piece[0],y:y+piece[1],z] = legopiece
                            self.__add_connections(filled,legopiece)
                            self.pieces.append(legopiece)

        print(np.sum(filled==-1))

        filled[np.where(filled==-1)] = 0

        for z in range(self.z_dim):
            piece_list = []
            
            if (z % 2 == 0):
                piece_list = self.piece_dimension_list
            else:
                piece_list = self.piece_dimension_list_odd

            for piece in piece_list:
                for x in range(self.x_dim):
                    for y in range(self.y_dim):
                        if self.__can_place(filled,piece,x,y,z):
                            legopiece = LegoPiece(piece,(x,y,z))
                            filled[x:x+piece[0],y:y+piece[1],z] = legopiece
                            self.pieces_visual_array[x:x+piece[0],y:y+piece[1],z] = legopiece.color
                            self.__add_connections(filled,legopiece)
                            self.pieces.append(legopiece)

        assert((filled != 1).all() and (filled != -1).all())
        self.pieces_visual_no_overhang = self.pieces_visual_array
        self.pieces_visual_no_overhang[np.where(self.pieces_visual_no_overhang == -1)] = 0
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
            if ((np.abs(filled[x:(x+dimensions[0]),y:(y+dimensions[1]),z]) == 1).all() 
            and x+dimensions[0] - 1<filled.shape[0] and y+dimensions[1] - 1 < filled.shape[1] and (filled[x:(x+dimensions[0]),y:(y+dimensions[1]),z] == 1).any()):
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

        array_3d = np.zeros((self.size,self.size,height))

        print(array_3d.shape)

        for i in range(self.coords_df.shape[0]):
            x = list(self.coords_df['x'])[i]
            y = list(self.coords_df['y'])[i]
            z = list(self.coords_df['z'])[i]

            x_ints = get_integers_around(x,self.buffersize)
            y_ints = get_integers_around(y,self.buffersize)
            z_ints = get_integers_around(z,self.buffersize)

            for x2 in x_ints:
                for y2 in y_ints:
                    for z2 in z_ints:
                        if x2 in range(0,self.size) and y2 in range(0,self.size) and z2 in range(0,height):
                            array_3d[x2][y2][z2] = 1
                            
        return array_3d

    def visualize_array(self):
        data = self.array
        axes=list(self.array.shape)
        alpha = 0.9

        colors = np.empty(axes + [4], dtype=np.float32)
        colors[:] = [1, 0, 0, alpha]  # red

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111,projection='3d')
        ax.voxels(data, facecolors=colors, edgecolors='grey')
        plt.show()

    def visualize_piece_array(self):
        data = self.pieces_visual_array
        voxels = data != 0

        discrete_colors = ['red', 'white', 'yellow']  # Colors for discrete ranges
        discrete_cmap = ListedColormap(discrete_colors)

        continuous_cmap = LinearSegmentedColormap.from_list('smooth', ['yellow', 'green', 'blue'])

        # Step 3: Define boundaries and normalization
        boundaries = [-0.5,0.5]  # Discrete range boundaries
        norm_discrete = BoundaryNorm(boundaries, len(discrete_colors), extend='both')
        norm_continuous = Normalize(vmin=0, vmax=256)  # Linear gradient in the range [0, 2]

        # Start with discrete colors
        facecolors = discrete_cmap(norm_discrete(data))

        # Overlay the continuous region
        continuous_mask = (data > 0.5)  # Region for continuous mapping
        facecolors[continuous_mask] = continuous_cmap(norm_continuous(data[continuous_mask]))

        # Make non-visible voxels transparent
        voxels = data != 0  # Define visible voxels
        facecolors[~voxels] = [0, 0, 0, 0]  # Transparent for invisible voxels

        # Step 5: Plot the voxels
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')
        ax.voxels(voxels, facecolors=facecolors, edgecolor='k')

        plt.show()

    def visualize_piece_array_no_overhang(self):
        data = self.pieces_visual_no_overhang
        voxels = data != 0

        discrete_colors = ['red', 'white', 'yellow']  # Colors for discrete ranges
        discrete_cmap = ListedColormap(discrete_colors)

        continuous_cmap = LinearSegmentedColormap.from_list('smooth', ['yellow', 'green', 'blue'])

        # Step 3: Define boundaries and normalization
        boundaries = [-0.5,0.5]  # Discrete range boundaries
        norm_discrete = BoundaryNorm(boundaries, len(discrete_colors), extend='both')
        norm_continuous = Normalize(vmin=0, vmax=256)  # Linear gradient in the range [0, 2]

        # Start with discrete colors
        facecolors = discrete_cmap(norm_discrete(data))

        # Overlay the continuous region
        continuous_mask = (data > 0.5)  # Region for continuous mapping
        facecolors[continuous_mask] = continuous_cmap(norm_continuous(data[continuous_mask]))

        # Make non-visible voxels transparent
        voxels = data != 0  # Define visible voxels
        facecolors[~voxels] = [0, 0, 0, 0]  # Transparent for invisible voxels

        # Step 5: Plot the voxels
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')
        ax.voxels(voxels, facecolors=facecolors, edgecolor='k')

        plt.show()

    def make_overhang_slices(self):
        folder_path = Path(f'{self.uniprot}/slices')
        if folder_path.exists():
            shutil.rmtree(f'{self.uniprot}/slices')
        folder_path.mkdir(parents=True,exist_ok=True)

        colors = ['red','white','green']  # Colors for each range
        cmap = ListedColormap(colors)

        # Step 2: Define boundaries for discrete ranges
        boundaries = [-0.5,0.5]  # The boundaries of the bins
        norm = BoundaryNorm(boundaries, len(colors), extend='both')

        for i in range(self.overhangarray.shape[2]):
            image = self.overhangarray[:,:,i]
            
            plt.matshow(image,cmap=cmap, norm=norm)
            plt.xticks(np.arange(-0.5, image.shape[0], 1), labels=np.arange(0, image.shape[0]+1, 1))
            plt.yticks(np.arange(-0.5, image.shape[1], 1), labels=np.arange(0, image.shape[1]+1, 1))
            plt.grid(which='major', color='black', linestyle='-', linewidth=1)
            plt.grid(which='minor', color='gray', linestyle='--', linewidth=0.5)

            plt.xlabel('x')
            plt.ylabel('y')
            plt.savefig(f'{self.uniprot}/slices/{i}.png')
            plt.clf()

    def make_lego_slices(self):
        self.pieces_visual_array[np.where(np.abs(self.pieces_visual_array) < 0.5)] = 0
        
        folder_path = Path(f'{self.uniprot}/legoslices')
        if folder_path.exists():
            shutil.rmtree(f'{self.uniprot}/legoslices')
        folder_path.mkdir(parents=True,exist_ok=True)

        discrete_colors = ['red', 'white','green','yellow']  # Colors for discrete ranges
        discrete_cmap = ListedColormap(discrete_colors)

        # Step 2: Define the continuous colormap (smooth transition for positive values)
        continuous_cmap = LinearSegmentedColormap.from_list('smooth', ['blue', 'yellow'])

        # Step 3: Define boundaries and normalizations
        boundaries = [-0.5, 0.5,1.5]  # Edges for discrete ranges
        norm_discrete = BoundaryNorm(boundaries, len(discrete_colors), extend='both')
        norm_continuous = plt.Normalize(2, 255) 

    
        for i in range(self.pieces_visual_array.shape[2]):
            image = self.pieces_visual_array[:,:,i]
            
            if i == 8:
                print(image)
            
            # Plot
            plt.matshow(image, cmap=discrete_cmap, norm=norm_discrete)
            mask = image > 1.5  # Mask for the continuous section
            masked_data = np.ma.masked_where(~mask, image)
            ax = plt.gca()
            mappable_continuous = ax.matshow(masked_data, cmap=continuous_cmap, norm=norm_continuous)

            plt.xticks(np.arange(-0.5, image.shape[1], 1), labels=np.arange(0, image.shape[1]+1, 1))
            plt.yticks(np.arange(-0.5, image.shape[0], 1), labels=np.arange(0, image.shape[0]+1, 1))
            plt.grid(which='major', color='black', linestyle='-', linewidth=1)
            plt.grid(which='minor', color='gray', linestyle='--', linewidth=0.5)

            plt.xlabel('x')
            plt.ylabel('y')
            plt.savefig(f'{self.uniprot}/legoslices/{i}.png')
            plt.clf()

    def make_lego_slices_no_overhang(self):
        self.pieces_visual_no_overhang[np.where(np.abs(self.pieces_visual_no_overhang) < 0.5)] = 0
        folder_path = Path(f'{self.uniprot}/legoslices_no_overhang')
        if folder_path.exists():
            shutil.rmtree(f'{self.uniprot}/legoslices_no_overhang')
        folder_path.mkdir(parents=True,exist_ok=True)
        discrete_colors = ['red', 'white','green','yellow']  # Colors for discrete ranges
        discrete_cmap = ListedColormap(discrete_colors)

        # Step 2: Define the continuous colormap (smooth transition for positive values)
        continuous_cmap = LinearSegmentedColormap.from_list('smooth', ['blue', 'yellow'])

        # Step 3: Define boundaries and normalizations
        boundaries = [-0.5, 0.5,1.5]  # Edges for discrete ranges
        norm_discrete = BoundaryNorm(boundaries, len(discrete_colors), extend='both')
        norm_continuous = plt.Normalize(2, 255) 

    
        for i in range(self.pieces_visual_no_overhang.shape[2]):
            image = self.pieces_visual_no_overhang[:,:,i]
            
            if i == 8:
                print(image)
            
            # Plot
            plt.matshow(image, cmap=discrete_cmap, norm=norm_discrete)
            mask = image > 1.5  # Mask for the continuous section
            masked_data = np.ma.masked_where(~mask, image)
            ax = plt.gca()
            mappable_continuous = ax.matshow(masked_data, cmap=continuous_cmap, norm=norm_continuous)

            plt.xticks(np.arange(-0.5, image.shape[1], 1), labels=np.arange(0, image.shape[1]+1, 1))
            plt.yticks(np.arange(-0.5, image.shape[0], 1), labels=np.arange(0, image.shape[0]+1, 1))
            plt.grid(which='major', color='black', linestyle='-', linewidth=1)
            plt.grid(which='minor', color='gray', linestyle='--', linewidth=0.5)

            plt.xlabel('x')
            plt.ylabel('y')
            plt.savefig(f'{self.uniprot}/legoslices_no_overhang/{i}.png')
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
    
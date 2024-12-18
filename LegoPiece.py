
import numpy.random as rand

class LegoPiece:
    def __init__(self, dimensions, coordinates):
        self.length = dimensions[0]
        self.width = dimensions[1]
        self.coordinates = coordinates
        self.connected_set = set()
        self.pieceid=f"{dimensions}at{coordinates}"

        self.color = rand.randint(2,255)

    def add_connection(self,piece):
        self.connected_set.add(piece)

    def connections(self):
        return self.connected_set
    
    def get_coordinates(self):
        return self.coordinates
    
    def get_dimensions(self):
        return (self.length,self.width)
    
    def piece_id(self):
        return self.pieceid
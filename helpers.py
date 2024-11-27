import math


def get_integers_around(x,buffer):
    lower = x - buffer/2
    upper = x + buffer/2

    return list(range(math.ceil(lower),math.ceil(upper)))


def can_place(array,dimensions,x,y,z):
    try:
        for i in range(dimensions[0]):
            for j in range(dimensions[1]):
                if array[i,j,z] != 1:
                    return False
        return True
    except:
        return False
    
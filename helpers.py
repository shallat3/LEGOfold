import math


def get_integers_around(x,buffer):
    lower = x - buffer/2
    upper = x + buffer/2

    return list(range(math.ceil(lower),math.ceil(upper)))
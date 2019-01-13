import numpy as np

def mask(h,w, center, radius):
    # x = positions[:, 0]
    # y = positions[:, 1]
    # w = x.max() - x.min()
    # h = y.max() - y.min()
    if center is None: # use the middle of the slab
            center = [int(w/2), int(h/2)]
            
    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    m = dist_from_center <= radius
    return m
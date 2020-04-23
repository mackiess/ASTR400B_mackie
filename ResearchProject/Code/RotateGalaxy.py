import numpy as np

# Function will rotate M31 to view it edge on
# From Lab 7 Part B
def EdgeRotate(posInit, velInit):
    # Inputs:
    #    posInit = initial [x, y, z] array
    #    velInit = initial [vx, vy, vz] array
    # Returns:
    #    posNew  = [x, y, z] array with coordinates rotated to view M31 edge on 
    #    velNew  = [vx, vy, vz] array with coordinates rotated to view M31 edge on 
    
    # compute the angular momentum
    L = np.sum(np.cross(posInit,velInit), axis=0)
    # normalize the vector
    L_norm = L/np.sqrt(np.sum(L**2))

    # Set up rotation matrix to map L_norm to z unit vector (disk in xy-plane)
    
    # z unit vector
    z_norm = np.array([0, 0, 1])
    
    # cross product between L and z
    vv = np.cross(L_norm, z_norm)
    s = np.sqrt(np.sum(vv**2))
    
    # dot product between L and z 
    c = np.dot(L_norm, z_norm)
    
    # rotation matrix
    I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    v_x = np.array([[0, -vv[2], vv[1]], [vv[2], 0, -vv[0]], [-vv[1], vv[0], 0]])
    R = I + v_x + np.dot(v_x, v_x)*(1 - c)/s**2

    # Rotate coordinate system
    posNew = np.dot(R, posInit.T).T
    velNew = np.dot(R, velInit.T).T
    
    return posNew, velNew
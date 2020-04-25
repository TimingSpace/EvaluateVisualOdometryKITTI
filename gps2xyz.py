#ref https://stackoverflow.com/questions/18759601/converting-lla-to-xyz
# the code for calculation are copied from https://stackoverflow.com/questions/18759601/converting-lla-to-xyz
# 

import glob
import numpy as np
# gps from lat lon
# gps to lat lon
import math

R = 6378137
f_inv = 298.257224
f = 1.0 / f_inv
e2 = 1 - (1 - f) * (1 - f)



def gps_to_ecef(latitude, longitude, altitude):
    # (lat, lon) in WSG-84 degrees
    # h in meters
    cosLat = math.cos(latitude * math.pi / 180)
    sinLat = math.sin(latitude * math.pi / 180)

    cosLong = math.cos(longitude * math.pi / 180)
    sinLong = math.sin(longitude * math.pi / 180)

    c = 1 / math.sqrt(cosLat * cosLat + (1 - f) * (1 - f) * sinLat * sinLat)
    s = (1 - f) * (1 - f) * c

    x = (R*c + altitude) * cosLat * cosLong
    y = (R*c + altitude) * cosLat * sinLong
    z = (R*s + altitude) * sinLat

    return x, y, z

# ecef2enu
def ecef_to_enu(x, y, z, latRef, longRef, altRef):

    cosLatRef = math.cos(latRef * math.pi / 180)
    sinLatRef = math.sin(latRef * math.pi / 180)

    cosLongRef = math.cos(longRef * math.pi / 180)
    sinLongRef = math.sin(longRef * math.pi / 180)

    cRef = 1 / math.sqrt(cosLatRef * cosLatRef + (1 - f) * (1 - f) * sinLatRef * sinLatRef)

    x0 = (R*cRef + altRef) * cosLatRef * cosLongRef
    y0 = (R*cRef + altRef) * cosLatRef * sinLongRef
    z0 = (R*cRef*(1-e2) + altRef) * sinLatRef

    xEast = (-(x-x0) * sinLongRef) + ((y-y0)*cosLongRef)

    yNorth = (-cosLongRef*sinLatRef*(x-x0)) - (sinLatRef*sinLongRef*(y-y0)) + (cosLatRef*(z-z0))

    zUp = (cosLatRef*cosLongRef*(x-x0)) + (cosLatRef*sinLongRef*(y-y0)) + (sinLatRef*(z-z0))

    return xEast, yNorth, zUp

def geodetic_to_enu(lat, lon, h, lat_ref, lon_ref, h_ref):
    x, y, z = gps_to_ecef(lat, lon, h)
    return ecef_to_enu(x, y, z, lat_ref, lon_ref, h_ref)


def gps2xyz(gps_from,gps_to):
    return geodetic_to_enu(gps_to[0],gps_to[1],gps_to[2],gps_from[0],gps_from[1],gps_from[2])   

def get_rotation(yaw):
    cy   = np.cos(yaw)
    sy   = np.sin(yaw)
    rotation_matrix=np.matrix([[cy,-sy,0],[sy,cy,0],[0,0,1]])
    return rotation_matrix
def main():
    import sys
    import matplotlib.pyplot as plt
    gps_files = sorted(glob.glob(sys.argv[1]+'/*'))
    gps_ref = np.loadtxt(gps_files[0])
    gps_ref_lla = gps_ref[0:3]
    yaw_ref     = gps_ref[5]
    rotation_matrix_ref =get_rotation(yaw_ref)
    rotation_matrix_zx=np.matrix([[0,-1,0],[1,0,0],[0,0,1]])
    rotation_matrix_I = rotation_matrix_ref.I*rotation_matrix_zx
    xyzs=[]
    pose = np.matrix(np.eye(4))
    poses = []
    coor= np.matrix(np.eye(4))
    coor[0:3,0:3] = np.array([[1,0,0],[0,0,-1],[0,1,0]])
    for i in range(0,len(gps_files)):

        gps_cur = np.loadtxt(gps_files[i])
        gps_cur_lla = gps_cur[0:3]
        yaw         = gps_cur[5]
        x,y,z =  gps2xyz(gps_ref_lla,gps_cur_lla)
        rotation = get_rotation(yaw)
        rotation_ = rotation_matrix_ref.I*rotation
        pose[0:3,0:3] = rotation_
        xyz =np.array([[x],[y],[z]])
        xyz_r = rotation_matrix_I@xyz
        pose[0:3,3] = xyz_r
        pose = coor@pose@coor.I
        xyzs.append([x,y,z])
        pose_line = np.array(pose[0:3,:]).reshape(-1).copy()
        poses.append(pose_line)
    xyzs = np.array(xyzs)
    poses= np.array(poses)
    xyzs_r = (rotation_matrix_I@xyzs.transpose()).T
    np.savetxt('pose_from_gt_raw_02.txt',poses)


if __name__ == '__main__':
    main()

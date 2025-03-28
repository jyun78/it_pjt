import numpy as np
from edit_pdb import read_pdb,save_pdb,update_pdb_atoms
import math

## rotate around axis (x, y, z) by angle theta (degree)
def getQuaternion(v_axis, theta):

    x, y, z = v_axis
    norm = math.sqrt(x*x + y*y + z*z)
    x /= norm
    y /= norm
    z /= norm

    q0 = math.cos(math.radians(theta * 0.5))
    q1 = x * math.sin(math.radians(theta * 0.5))
    q2 = y * math.sin(math.radians(theta * 0.5))
    q3 = z * math.sin(math.radians(theta * 0.5))

    return np.array([q0, q1, q2, q3])

## Returns the rotation of a vector by a given quaternion.
## Shamelessly from Greg's quatvec function in functors.h in his MC code
# input quat and vec are arrays of four numbers: [q0, q1, q2, q3]
def quatVec(quat, vec):

    v0 = (quat[0]**2 + quat[1]**2 - quat[2]**2 - quat[3]**2) * vec[0] \
         + (2 * quat[1] * quat[2] - 2 * quat[0] * quat[3]) * vec[1] \
         + (2 * quat[3] * quat[1] + 2 * quat[0] * quat[2]) * vec[2]
    
    v1 = (2 * quat[1] * quat[2] + 2 * quat[0] * quat[3]) * vec[0] \
         + (quat[0]**2 - quat[1]**2 + quat[2]**2 - quat[3]**2) * vec[1] \
         + (2 * quat[2] * quat[3] - 2 * quat[0] * quat[1]) * vec[2]
    
    v2 = (2 * quat[1] * quat[3] - 2 * quat[0] * quat[2]) * vec[0] \
         + (2 * quat[0] * quat[1] + 2 * quat[2] * quat[3]) * vec[1] \
         + (quat[0]**2 - quat[1]**2 - quat[2]**2 + quat[3]**2) * vec[2]
    
    return np.array([v0, v1, v2])

## Rotates a list of vectors according to the quaternion quat
# input quat and vec are arrays of four numbers: [q0, q1, q2, q3]
def rotateFrame(quat, veclist):

    vecarr = np.asarray(veclist)
    new_veclist = [quatVec(quat, vec) for vec in vecarr]
    return np.asarray(new_veclist)

## Returns the product of two quaternions.
# input quats are arrays of four numbers: [q0, q1, q2, q3]
def quatMultiply(quatA, quatB):

    qAs, qAv = quatA[0], np.array(quatA[1:])
    qBs, qBv = quatB[0], np.array(quatB[1:])
    
    qABs = qAs * qBs - np.dot(qAv, qBv)
    qABv = qAs * qBv + qBs * qAv + np.cross(qAv, qBv)
    
    return np.array([qABs, qABv[0], qABv[1], qABv[2]])

## rdb파일 읽고, center of mass 계산
## input: file_path(str), output: ceter of mass(tuple: (x,y,z))
def calculate_center_of_mass(pdb_atoms):
    if not pdb_atoms:
        print("No atoms found in the PDB file.")
        return None

    # 모든 x, y, z 좌표를 리스트로 저장
    x_coords = np.array([atom['x'] for atom in pdb_atoms])
    y_coords = np.array([atom['y'] for atom in pdb_atoms])
    z_coords = np.array([atom['z'] for atom in pdb_atoms])

    # 중심 좌표 계산 (평균)
    center_x = np.mean(x_coords)
    center_y = np.mean(y_coords)
    center_z = np.mean(z_coords)
    return center_x, center_y, center_z

def rotate_pdb(pdb_atoms, rotation_axis, rotation_angle, custom_center=None):
    original_center_x, original_center_y, original_center_z = calculate_center_of_mass(pdb_atoms)
    original_center = np.array([original_center_x, original_center_y, original_center_z])

    # Define atom_coords regardless of custom_center
    atom_coords = np.array([[atom['x'], atom['y'], atom['z']] for atom in pdb_atoms])

    if custom_center is None:
        center = original_center
    else:
        center = np.array(custom_center)
        translation = center - original_center
        atom_coords += translation  

    centered_coords = atom_coords - center
    quat = getQuaternion(rotation_axis, rotation_angle)
    rotated_coords = rotateFrame(quat, centered_coords)
    final_coords = rotated_coords + center
    rot_atoms = update_pdb_atoms(pdb_atoms, final_coords)
    
    return rot_atoms
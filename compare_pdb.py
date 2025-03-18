import numpy as np
from edit_pdb import read_pdb
from scipy.spatial.transform import Rotation as R

## 두 원자 좌표간 RMSD (Root Mean Square Deviation) 단, 같은 길이만 가능
def calculate_rmsd(pdb1, pdb2):

    if len(pdb1) != len(pdb2):
        raise ValueError("The two PDB files must have the same number of atoms.")

    coords1 = np.array([[atom['x'], atom['y'], atom['z']] for atom in pdb1])
    coords2 = np.array([[atom['x'], atom['y'], atom['z']] for atom in pdb2])

    rmsd = np.sqrt(np.mean(np.sum((coords1 - coords2) ** 2, axis=1)))
    return rmsd

## 두 원자 사이 각도 계산 (SVD사용)
def calculate_rotation_angle(pdb1, pdb2):

    if len(pdb1) != len(pdb2):
        raise ValueError("The two PDB files must have the same number of atoms.")

    coords1 = np.array([[atom['x'], atom['y'], atom['z']] for atom in pdb1])
    coords2 = np.array([[atom['x'], atom['y'], atom['z']] for atom in pdb2])

    centroid1 = np.mean(coords1, axis=0)
    centroid2 = np.mean(coords2, axis=0)
    centered_coords1 = coords1 - centroid1
    centered_coords2 = coords2 - centroid2

    H = np.dot(centered_coords1.T, centered_coords2)
    U, S, Vt = np.linalg.svd(H)
    rotation_matrix = np.dot(Vt.T, U.T)

    # Convert rotation matrix to quaternion
    r = R.from_matrix(rotation_matrix)
    quat = r.as_quat()  # [x, y, z, w]
    rotation_angle = np.degrees(2 * np.arccos(quat[3]))  #  degrees

    return rotation_angle

if __name__ == "__main__":
    file1 = "RK718_full.pdb"  # 원본 PDB 파일
    file2 = "rotated_example_x.pdb"  # 회전된 PDB 파일

    pdb1 = read_pdb(file1)
    pdb2 = read_pdb(file2)

    rmsd_value = calculate_rmsd(pdb1, pdb2)
    rotation_angle = calculate_rotation_angle(pdb1, pdb2)

    print(f"RMSD between the two structures: {rmsd_value:.3f} Å")
    print(f"Rotation angle between the two structures: {rotation_angle:.3f} degrees")

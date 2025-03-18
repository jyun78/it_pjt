import numpy as np
from cal_quternion import rotate_pdb
from edit_pdb import save_pdb

if __name__ == "__main__":
    file_path = "RK718_full.pdb"  # 회전할 PDB 파일 경로
    output_file = "rotated_example_x10.pdb"  # 저장할 파일 경로
    rotation_axis = [0, 1, 0]  # 어떤 축으로 회전할지
    rotation_angle = -90  # 회전 각도 (degree)
    rot_atoms=rotate_pdb(file_path, rotation_axis, rotation_angle,(10,0,0))
    save_pdb(output_file, rot_atoms)

import os,re,glob
from datetime import datetime

def read_pdb(file_path):

    pdb_data = []
    
    try:
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith("ATOM"):
                    atom_data = {
                        'atom_serial': int(line[6:11].strip()),
                        'atom_name': line[12:16].strip(),
                        'residue_name': line[17:20].strip(),
                        'chain': line[21].strip(),
                        'residue_seq': int(line[22:26].strip()),
                        'x': float(line[30:38].strip()),
                        'y': float(line[38:46].strip()),
                        'z': float(line[46:54].strip())
                    }
                    pdb_data.append(atom_data)
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
    
    return pdb_data


def save_pdb(polytype, pdb_atoms):
    # 현재 작업 디렉토리 내 "output" 폴더 생성 (없으면)
    output_dir = os.path.join(os.getcwd(), "output")
    os.makedirs(output_dir, exist_ok=True)
    
    # 오늘 날짜 문자열 (YYYYMMDD 형식)
    today_str = datetime.now().strftime("%Y%m%d")
    
    # polytype 기반의 패턴으로 검색
    pattern = os.path.join(output_dir, f"{polytype}_{today_str}_*.pdb")
    existing_files = glob.glob(pattern)
    
    # 정규식도 polytype 반영해서 설정
    regex = re.compile(rf"{re.escape(polytype)}_{today_str}_(\d+)\.pdb")
    max_index = 0
    for file in existing_files:
        basename = os.path.basename(file)
        m = regex.match(basename)
        if m:
            num = int(m.group(1))
            if num > max_index:
                max_index = num

    next_index = max_index + 1
    output_file = os.path.join(output_dir, f"{polytype}_{today_str}_{next_index}.pdb")

    with open(output_file, "w") as f:
        for atom in pdb_atoms:
            f.write(f"ATOM  {atom['atom_serial']:5d} {atom['atom_name']:^4} {atom['residue_name']:3} {atom['chain']:1}"
                    f"{atom['residue_seq']:4d}    {atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}\n")

    print(f"PDB file saved as: {output_file}")

def update_pdb_atoms(pdb_atoms, new_coords):
    for atom, new_coord in zip(pdb_atoms, new_coords):
        atom['x'], atom['y'], atom['z'] = new_coord  # 좌표만 업데이트
    return pdb_atoms
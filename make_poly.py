import numpy as np
from edit_pdb import read_pdb,save_pdb
from cal_quternion import *
import copy

def get_polycor(polytype):
    """
    정다면체의 종류에 따라 꼭지점 좌표(단위 벡터)를 반환합니다.
    
    polytype: 문자열 (대소문자 구분 없이)
        지원 종류:
         - "tetrahedron" 또는 "정사면체"
         - "octahedron" 또는 "정팔면체"
         - "icosahedron" 또는 "정이십면체"
         - "dodecahedron" 또는 "정십이면체"  (여기서는 dodecahedron으로 처리)
         
    반환값: 각 꼭지점 좌표가 [x, y, z]인 리스트.
    """
    polytype = polytype.lower()
    
    if polytype in ['tetrahedron', '정사면체']:
        # 정사면체의 꼭지점 (예시 좌표, 단위벡터로 정규화)
        vertices = np.array([
            [1, 1, 1],
            [1, -1, -1],
            [-1, 1, -1],
            [-1, -1, 1]
        ], dtype=float)
        # 각 좌표를 단위벡터로 정규화
        vertices = np.array([v / np.linalg.norm(v) for v in vertices])
        return vertices.tolist()
    
    elif polytype in ['octahedron', '정팔면체']:
        # 정팔면체의 꼭지점 (이미 단위벡터)
        vertices = np.array([
            [1, 0, 0],
            [-1, 0, 0],
            [0, 1, 0],
            [0, -1, 0],
            [0, 0, 1],
            [0, 0, -1]
        ], dtype=float)
        return vertices.tolist()
    
    elif polytype in ['icosahedron', '정이십면체']:
        # 정이십면체: 황금비(phi)를 이용
        phi = (1 + np.sqrt(5)) / 2
        vertices = []
        # (0, ±1, ±phi)
        for sign1 in [1, -1]:
            for sign2 in [1, -1]:
                vertices.append([0, sign1, sign2 * phi])
        # (±1, ±phi, 0)
        for sign1 in [1, -1]:
            for sign2 in [1, -1]:
                vertices.append([sign1, sign2 * phi, 0])
        # (±phi, 0, ±1)
        for sign1 in [1, -1]:
            for sign2 in [1, -1]:
                vertices.append([sign1 * phi, 0, sign2])
        # 단위벡터로 정규화
        vertices = np.array(vertices, dtype=float)
        vertices = np.array([v / np.linalg.norm(v) for v in vertices])
        return vertices.tolist()
    
    elif polytype in ['dodecahedron', '정십이면체']:
        # 정십이면체: dodecahedron (12면체가 아니라 20개의 꼭지점을 가짐)
        phi = (1 + np.sqrt(5)) / 2
        vertices = []
        # 8개: (±1, ±1, ±1)
        for i in [-1, 1]:
            for j in [-1, 1]:
                for k in [-1, 1]:
                    vertices.append([i, j, k])
        # 12개: (0, ±1/phi, ±phi), (±1/phi, ±phi, 0), (±phi, 0, ±1/phi)
        for j in [-1, 1]:
            for k in [-1, 1]:
                vertices.append([0, j/phi, k*phi])
        for i in [-1, 1]:
            for k in [-1, 1]:
                vertices.append([i/phi, k*phi, 0])
        for i in [-1, 1]:
            for j in [-1, 1]:
                vertices.append([i*phi, 0, j/phi])
        # 단위벡터로 정규화
        vertices = np.array(vertices, dtype=float)
        vertices = np.array([v / np.linalg.norm(v) for v in vertices])
        return vertices.tolist()
    
    else:
        raise ValueError("Unknown polyhedron type: {}".format(polytype))

def scale_poly(poly_coords, scale):
    """
    poly_coords에 있는 각 좌표(단위벡터)를 scale 만큼 스케일링합니다.
    
    poly_coords: 각 꼭지점 좌표 ([x, y, z])의 리스트.
    scale: 스케일 인자 (예: 120)
    
    반환값: 스케일이 적용된 좌표 목록 (각 좌표는 numpy 배열)
    """
    scaled_coords = [np.array(coord) * scale for coord in poly_coords]
    return scaled_coords

def separate_three(pdb_atoms):
    return [atom for atom in pdb_atoms if atom.get('chain') == 'A']

def rotate_vector(vec, axis, angle_deg):
    """
    vec (3D 벡터)를 axis(축)로 angle_deg(도 단위)만큼 회전(Rodrigues 공식을 사용).
    방향벡터이므로 평행이동 없음.
    """
    # 1) 라디안 변환 및 축 정규화
    angle_rad = np.radians(angle_deg)
    axis = axis / (np.linalg.norm(axis) + 1e-12)

    # 2) 로드리게스(Rodrigues) 공식
    #    v_rot = v*cosθ + (k×v)*sinθ + k*(k·v)*(1-cosθ)
    v = vec
    k = axis
    cos_t = np.cos(angle_rad)
    sin_t = np.sin(angle_rad)

    v_rot = (v * cos_t
             + np.cross(k, v) * sin_t
             + k * np.dot(k, v) * (1.0 - cos_t))
    return v_rot

# #def trans2polyhedron(pdb_atoms, polytype):
    
#     # (1) 정사면체 꼭지점 좌표 & 스케일 계산
#     poly_cor = get_polycor(polytype)
#     scale_cor = scale_poly(poly_cor, 120)
    
#     #  단백질이 이미 z축 주위로 3겹(120°) 대칭을 이루고 있음. 이때 하나는 x축과 얼라인 되어있음음
#     protein_axis = np.array([0.0, 0.0, 1.0])

#     for i, vertex in enumerate(scale_cor):   
#         next_i = (i + 1) % len(scale_cor)  # 인접한 꼭짓점 (마지막이면 첫 번째로 순환)
#         next_vertex = scale_cor[next_i]
#         rotation_axis = np.cross(protein_axis, vertex)
#         axis_norm = np.linalg.norm(rotation_axis)
        
#         if axis_norm < 1e-7:
#             # 평행한 경우 (동일 또는 반대 방향)
#             dot_val = np.dot(protein_axis, vertex)
#             if dot_val < 0:
#                 # 반대 방향 -> 180도 회전
#                 rotation_axis = np.array([1.0, 0.0, 0.0])
#                 rotation_angle = 180.0
#             else:
#                 # 동일 방향 -> 회전 불필요
#                 rotation_axis = np.array([1.0, 0.0, 0.0])
#                 rotation_angle = 0.0
#         else:
#             # 교차축과 회전각 계산
#             rotation_axis /= axis_norm
#             dot_val = np.dot(protein_axis, vertex)
#             cos_val = dot_val / (np.linalg.norm(protein_axis) * np.linalg.norm(vertex))
#             cos_val = max(min(cos_val, 1.0), -1.0)
#             rotation_angle = np.degrees(np.arccos(cos_val))
        
#         # 원본 복사
#         atoms_copy = copy.deepcopy(pdb_atoms)
        
#         # (가) 첫 회전 + 무게중심 → vertex 이동
#         rotated_atoms = rotate_pdb(
#             atoms_copy,
#             rotation_axis,
#             rotation_angle,
#             custom_center=vertex
#         )
#         new_axis = rotate_vector(protein_axis, rotation_axis, rotation_angle)

#         separates = separate_three(rotated_atoms)
#         #print(separates)
#         x = calculate_center_of_mass(separates)
#         theta_rad,final_angle = rotation_needed(vertex,x,vertex,next_vertex)
#         atoms_copy_2 = copy.deepcopy(rotated_atoms)
#         final_atoms = rotate_pdb(
#             atoms_copy_2,
#             new_axis,
#             final_angle
#         )
#         print(new_axis,final_angle,x)

#         save_pdb(polytype,final_atoms)
    
def trans2polyhedron(pdb_atoms, polytype, vertex_index=0):
    """
    pdb_atoms: 원본 PDB 원자 리스트 (각 원소는 딕셔너리)
    polytype: 정다면체 종류 (예: "tetrahedron")
    vertex_index: 사용할 꼭지점 번호 (예: 0,1,2,...)
    
    반환: 해당 꼭지점에 대해 변환된 pdb_atoms (복제본)
    """
    # (1) 정사면체 꼭지점 좌표 & 스케일 계산
    poly_cor = get_polycor(polytype)
    scale_cor = scale_poly(poly_cor, 120)
    
    # 단백질은 이미 z축 주위로 3겹(120°) 대칭을 이루고 있으며, 
    # 이때 하나는 x축과 얼라인되어 있다고 가정.
    protein_axis = np.array([0.0, 0.0, 1.0])
    
    # 입력된 vertex_index가 유효한지 체크
    if vertex_index < 0 or vertex_index >= len(scale_cor):
        raise ValueError("vertex_index가 올바르지 않습니다.")
    
    vertex = scale_cor[vertex_index]
    next_index = (vertex_index + 1) % len(scale_cor)  # 인접 꼭지점 (순환)
    next_vertex = scale_cor[next_index]
    
    # === (A) 첫 회전: z축 ↔ vertex 정렬 ===
    rotation_axis = np.cross(protein_axis, vertex)
    axis_norm = np.linalg.norm(rotation_axis)
    if axis_norm < 1e-7:
        dot_val = np.dot(protein_axis, vertex)
        if dot_val < 0:
            rotation_axis = np.array([1.0, 0.0, 0.0])
            rotation_angle = 180.0
        else:
            rotation_axis = np.array([1.0, 0.0, 0.0])
            rotation_angle = 0.0
    else:
        rotation_axis /= axis_norm
        dot_val = np.dot(protein_axis, vertex)
        cos_val = dot_val / (np.linalg.norm(protein_axis) * np.linalg.norm(vertex))
        cos_val = max(min(cos_val, 1.0), -1.0)
        rotation_angle = np.degrees(np.arccos(cos_val))
    
    atoms_copy = copy.deepcopy(pdb_atoms)
    # 첫 회전 + 무게중심을 vertex로 이동
    rotated_atoms = rotate_pdb(
        atoms_copy,
        rotation_axis,
        rotation_angle,
        custom_center=vertex
    )
    # protein_axis를 첫 회전으로 돌린 결과 (새로운 기준 축)
    new_axis = rotate_vector(protein_axis, rotation_axis, rotation_angle)
    
    # 분리된 3개의 영역(예: 3개의 헬릭스 부분)을 추출하고, 그 중심을 계산
    separates = separate_three(rotated_atoms)
    x = calculate_center_of_mass(separates)
    
    # vertex와 인접 꼭지점(next_vertex)을 이용해 추가 회전각 계산 (roll 보정)
    theta_rad, final_angle = rotation_needed(vertex, x, vertex, next_vertex)
    
    atoms_copy_2 = copy.deepcopy(rotated_atoms)
    final_atoms = rotate_pdb(
        atoms_copy_2,
        new_axis,
        final_angle,
        custom_center=vertex  # 두 번째 회전에서도 vertex를 중심으로 회전하도록 함
    )
    print("새로운 축:", new_axis, "추가 회전각:", final_angle, "분리 영역 중심:", x)
    
    save_pdb(polytype, final_atoms)

        
def rotation_needed(V, P, A, B):
    """
    V에서 P로 향하는 벡터와, A에서 B로 향하는 엣지 벡터 사이의 각도를 계산합니다.
    이 각도가 바로 V→P를 엣지 A→B 방향으로 일치시키기 위해 필요한 회전각입니다.
    
    Parameters
    ----------
    V : array-like, shape (3,)
        정사면체의 한 꼭지점 (회전의 기준점)
    P : array-like, shape (3,)
        V와 잇는 임의의 점
    A : array-like, shape (3,)
        엣지의 시작점
    B : array-like, shape (3,)
        엣지의 끝점
    
    Returns
    -------
    theta_rad : float
        회전각 (라디안 단위)
    theta_deg : float
        회전각 (도 단위)
    """
    # V→P 벡터 계산
    v = np.array(P, dtype=float) - np.array(V, dtype=float)
    # 엣지 A→B 벡터 계산
    e = np.array(B, dtype=float) - np.array(A, dtype=float)
    
    norm_v = np.linalg.norm(v)
    norm_e = np.linalg.norm(e)
    
    if norm_v == 0 or norm_e == 0:
        raise ValueError("입력 벡터 중 0벡터가 있습니다. 올바른 좌표를 입력하세요.")
    
    # 단위 벡터로 정규화
    v_unit = v / norm_v
    e_unit = e / norm_e
    
    # 두 단위 벡터 사이의 내적으로 코사인값 구하기
    dot_val = np.dot(v_unit, e_unit)
    # 수치오차로 인한 문제 방지를 위해 -1 ~ 1 범위로 제한
    dot_val = np.clip(dot_val, -1.0, 1.0)
    
    # 두 벡터 사이의 각도 (라디안)
    theta_rad = np.arccos(dot_val)
    # 도(degree) 단위로 변환
    theta_deg = np.degrees(theta_rad)
    
    return theta_rad, theta_deg

if __name__ == "__main__":
    file_path = "RK718_full.pdb"  # 원본 PDB 파일
    pdb_atoms = read_pdb(file_path)
    polytype='tetrahedron'
    trans2polyhedron(pdb_atoms, polytype,0)
    ##final angle이 제대로 안구해지고 있음. 



    

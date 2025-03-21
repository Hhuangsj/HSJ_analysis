import os
import MDAnalysis as mda
from MDAnalysis.coordinates.XTC import XTCWriter
import glob
import shutil

# 主目录
root_dir = "Coordinates_splite_pdb_p2"
output_xtc = "output.xtc"
failed_dir = "failed_pdbs"  # 存放失败的 PDB
os.makedirs(failed_dir, exist_ok=True)

# 找到所有子文件夹中的 PDB 文件
pdb_files = sorted(glob.glob(os.path.join(root_dir, "*", "*.pdb")))

if not pdb_files:
    print("未找到 PDB 文件")
else:
    # 以第一个 PDB 为拓扑
    u = mda.Universe(pdb_files[0])
    reference_n_atoms = u.atoms.n_atoms  # 记录参考原子数
    print(f"参考原子数: {reference_n_atoms}")

    # 轨迹写入
    with XTCWriter(output_xtc, n_atoms=reference_n_atoms) as xtc:
        for pdb in pdb_files:
            try:
                u = mda.Universe(pdb)
                if u.atoms.n_atoms != reference_n_atoms:
                    print(f"❌ {pdb} 原子数不匹配 ({u.atoms.n_atoms} vs {reference_n_atoms})，移动到 {failed_dir}")
                    shutil.move(pdb, os.path.join(failed_dir, os.path.basename(pdb)))
                    continue

                print(f"✅ Processing {pdb}")
                xtc.write(u.atoms)  # 逐帧写入 XTC
            except Exception as e:
                print(f"❌ 处理 {pdb} 时出错: {e}，移动到 {failed_dir}")
                shutil.move(pdb, os.path.join(failed_dir, os.path.basename(pdb)))

    print(f"✅ 轨迹文件已保存为 {output_xtc}")
    print(f"❗ 失败的 PDB 文件已移动到 {failed_dir}")

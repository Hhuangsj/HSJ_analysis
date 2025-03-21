import sys

def sort_pdb_by_residue(input_pdb, output_pdb):
    """按残基编号和氨基酸名称对 PDB 文件中的 ATOM 记录排序"""
    with open(input_pdb, "r") as f:
        lines = f.readlines()

    # 分离 ATOM 记录和其他记录
    atom_lines = [line for line in lines if line.startswith("ATOM")]
    other_lines = [line for line in lines if not line.startswith("ATOM")]

    # 按残基编号（resid）、氨基酸名称（resname）、原子名称（atom name）排序
    atom_lines.sort(key=lambda x: (int(x[22:26].strip()), x[17:20].strip(), x[12:16].strip()))

    # 写入排序后的 PDB 文件
    with open(output_pdb, "w") as f:
        f.writelines(other_lines)  # 先写入非 ATOM 记录
        f.writelines(atom_lines)   # 再写入排序后的 ATOM 记录

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("用法: python sort_pdb.py input.pdb output.pdb")
        sys.exit(1)

    input_pdb = sys.argv[1]
    output_pdb = sys.argv[2]

    sort_pdb_by_residue(input_pdb, output_pdb)
    print(f"排序完成，结果已保存至 {output_pdb}")


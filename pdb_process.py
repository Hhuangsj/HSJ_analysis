import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
import shutil
import logging
import tempfile
from tqdm import tqdm
import psutil

# 原始和目标主文件夹
input_root = "Coordinates_splite_pdb"
output_root = "Coordinates_splite_pdb_p1"
log_dir = "log"

# 确保输出目录和 log 目录存在
os.makedirs(output_root, exist_ok=True)
os.makedirs(log_dir, exist_ok=True)

# 配置日志
logging.basicConfig(filename='error.log', level=logging.ERROR,
                    format='%(asctime)s %(levelname)s: %(message)s')

def run_command(cmd, timeout=60):
    """运行子进程命令，并设置超时时间"""
    try:
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        logging.error(f"Timeout executing command: {' '.join(cmd)}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error executing command: {' '.join(cmd)}: {e}")

def process_pdb(input_pdb_path):
    """处理单个 PDB 文件"""
    base_name = os.path.splitext(os.path.basename(input_pdb_path))[0]

    # 使用临时目录管理中间文件
    with tempfile.TemporaryDirectory(dir=log_dir) as temp_dir:
        tmp_pdb = os.path.join(temp_dir, f"{base_name}_tmp.pdb")
        noH_pdb = os.path.join(temp_dir, f"{base_name}_noH.pdb")
        reH_pdb = os.path.join(temp_dir, f"{base_name}_reH.pdb")

        # 执行处理步骤
        run_command(["python", "sort_pdb.py", input_pdb_path, tmp_pdb])
        run_command(["pdb4amber", "-i", tmp_pdb, "-o", noH_pdb, "-y"])
        run_command(["pdb4amber", "-i", noH_pdb, "-o", reH_pdb, "--reduce"])

        # 移动最终文件到 output_root
        shutil.move(reH_pdb, os.path.join(output_root, os.path.basename(reH_pdb)))

def get_max_workers():
    """动态调整并行度，根据 CPU 负载决定"""
    cpu_load = psutil.cpu_percent(interval=1)
    if cpu_load > 80:  # 如果 CPU 负载超过 80%，减少并行度
        return max(1, os.cpu_count() // 4)
    else:
        return min(32, os.cpu_count() * 2)  # 限制最大并行度为 32

def main():
    # 获取所有 PDB 文件路径
    pdb_files = [os.path.join(root, file)
                 for root, _, files in os.walk(input_root)
                 for file in files if file.endswith(".pdb")]

    # 动态调整并行度
    max_workers = get_max_workers()
    print(f"Using {max_workers} workers")

    # 使用多进程并行处理
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_pdb, pdb_file): pdb_file for pdb_file in pdb_files}
        for future in tqdm(as_completed(futures), total=len(pdb_files)):
            try:
                future.result()  # 等待任务完成
            except Exception as e:
                logging.error(f"Error processing {futures[future]}: {e}")

if __name__ == "__main__":
    main()

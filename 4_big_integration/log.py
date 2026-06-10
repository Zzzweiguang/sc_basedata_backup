import subprocess
import os
import time

def run_script_in_background(script_path, log_file=None):
    """在后台运行Python脚本"""
    if not os.path.exists(script_path):
        print(f"错误: 脚本文件 {script_path} 不存在")
        return None
    
    # 直接运行Python脚本（假设已在正确环境中）
    cmd = ["python", script_path]
    
    # 设置输出处理
    stdout = subprocess.PIPE if log_file is None else open(log_file, 'w')
    stderr = subprocess.STDOUT if log_file is None else stdout
    
    # 启动进程
    process = subprocess.Popen(cmd, stdout=stdout, stderr=stderr, start_new_session=True)
    
    print(f"已启动脚本 {script_path}，进程ID: {process.pid}")
    return process

if __name__ == "__main__":
    # 配置参数 - 根据你的需求修改这些值
    script1_path = "/home/lixiangyu/zr/Annotate/ANNOTATE_new/4_big_integration/4_1_adata_X_sparse.py"   # 第一个要运行的脚本
    script2_path = "/home/lixiangyu/zr/Annotate/ANNOTATE_new/4_big_integration/4_2_adata_layers_sparse.py"  # 第二个要运行的脚本
    
    # 日志文件路径
    log1_path = "/home/lixiangyu/zr/Annotate/ANNOTATE_new/4_big_integration/log/4_1_adata_X_sparse.log"
    log2_path = "/home/lixiangyu/zr/Annotate/ANNOTATE_new/4_big_integration/log/4_2_adata_layers_sparse.log"
    
    print("开始运行第一个脚本...")
    process1 = run_script_in_background(script1_path, log1_path)
    
    if process1:
        # 等待第一个脚本完成
        print(f"等待脚本 {script1_path} 完成...")
        while True:
            return_code = process1.poll()
            if return_code is not None:
                print(f"脚本 {script1_path} 已完成，退出码: {return_code}")
                break
            time.sleep(1)
        
        # 运行第二个脚本
        print("开始运行第二个脚本...")
        process2 = run_script_in_background(script2_path, log2_path)
        
        if process2:
            print(f"脚本 {script2_path} 已在后台启动，进程ID: {process2.pid}")
            print("所有脚本已按顺序启动")
    
import matlab.engine
import os
import concurrent.futures
import time
from tqdm import tqdm # 进度条库

# ================= 配置区域 =================
# 1. 脚本目录
SCRIPT_DIR = r'D:\work\Task_Convergence\03_Terms_Calculation' 

# 2. 最大并行数（建议：内存 16G 设为 4，32G 设为 8）
MAX_WORKERS = 4 

# 3. 指定脚本名（留空则运行目录下所有 .m 文件）
TARGET_SCRIPTS = [] 
# ===========================================

def run_matlab_task(script_name):
    """运行单个脚本并返回详细状态"""
    result = {"name": script_name, "status": "Pending", "error": "", "duration": 0}
    eng = None
    try:
        # 启动引擎（-nodisplay 速度最快）
        eng = matlab.engine.start_matlab("-nodisplay -nosplash")
        eng.cd(SCRIPT_DIR, nargout=0)
        
        start_t = time.time()
        # 执行脚本，nargout=0 表示不要求返回值
        eng.eval(script_name, nargout=0)
        
        result["duration"] = time.time() - start_t
        result["status"] = "Success"
    except Exception as e:
        result["status"] = "Failed"
        # 捕获完整的 MATLAB 报错信息
        result["error"] = str(e)
    finally:
        if eng:
            eng.quit()
    return result

if __name__ == "__main__":
    # 自动检索脚本
    if not TARGET_SCRIPTS:
        TARGET_SCRIPTS = [f.replace('.m', '') for f in os.listdir(SCRIPT_DIR) if f.endswith('.m')]

    total_tasks = len(TARGET_SCRIPTS)
    print(f"🚀 准备并行处理 {total_tasks} 个 MATLAB 脚本...")
    print(f"⚙️  并行池大小: {MAX_WORKERS}")
    print("-" * 50)

    results = []
    
    # 使用 ThreadPoolExecutor（对于启动多进程引擎更稳定）
    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        # 提交所有任务
        future_to_script = {executor.submit(run_matlab_task, name): name for name in TARGET_SCRIPTS}
        
        # 使用 tqdm 显示进度条
        with tqdm(total=total_tasks, desc="执行进度", unit="脚本") as pbar:
            for future in concurrent.futures.as_completed(future_to_script):
                res = future.result()
                results.append(res)
                # 每完成一个，更新进度条
                pbar.update(1)
                if res["status"] == "Failed":
                    pbar.write(f"⚠️  检测到错误: {res['name']}")

    # --- 最终总结报告 ---
    print("\n" + "="*20 + " 运行总结报告 " + "="*20)
    success_count = sum(1 for r in results if r["status"] == "Success")
    failed_count = total_tasks - success_count
    
    for r in results:
        if r["status"] == "Success":
            print(f"✅ {r['name'].ljust(20)} | 耗时: {r['duration']:6.2f}s")
        else:
            print(f"❌ {r['name'].ljust(20)} | 状态: 失败")
            print(f"   👉 报错信息: {r['error']}\n")

    print("-" * 54)
    print(f"完成统计: 成功 {success_count} / 失败 {failed_count}")
    print("=" * 54)
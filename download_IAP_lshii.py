import os
import requests
from bs4 import BeautifulSoup
from urllib.parse import urljoin
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import random

# ================= 配置区域 =================

# 1. 保存路径
BASE_SAVE_DIR = r"D:\work\Ishii_v7.3.1_Data_05_24"

# 2. 下载目标
TARGETS = {
    "Temperature": "https://climate.mri-jma.go.jp/pub/ocean/ts/v7.3.1/temp/nc/",
    #"Salinity":    "https://climate.mri-jma.go.jp/pub/ocean/ts/v7.3.1/sal/nc/"
}

# 3. 并行下载数量 (建议 4-8，太多会被封IP或反而变慢)
MAX_WORKERS = 5

# 4. 年份过滤
START_YEAR = 2005
END_YEAR   = 2024
TARGET_YEARS = [str(y) for y in range(START_YEAR, END_YEAR + 1)]
FILE_EXTS = (".nc", ".gz")

# ===========================================

def get_session():
    """创建一个带有自动重试功能的 Session"""
    session = requests.Session()
    retry = Retry(
        total=5,                # 最多重试5次
        backoff_factor=1,       # 每次重试间隔时间 (1s, 2s, 4s...)
        status_forcelist=[500, 502, 503, 504], # 针对这些错误码重试
        allowed_methods=["HEAD", "GET", "OPTIONS"]
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    # 随机 User-Agent 防止被屏蔽
    session.headers.update({
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/115.0.0.0 Safari/537.36"
    })
    return session

def is_target_file(filename):
    if not filename.endswith(FILE_EXTS):
        return False
    for year in TARGET_YEARS:
        if year in filename:
            return True
    return False

def download_single_file(url, save_dir):
    """单个文件的下载逻辑"""
    file_name = url.split('/')[-1]
    save_path = os.path.join(save_dir, file_name)
    
    session = get_session()
    
    try:
        # 1. 获取文件大小 (HEAD 请求)
        head = session.head(url, timeout=30)
        total_size = int(head.headers.get('content-length', 0))
        
        # 2. 检查本地文件
        if os.path.exists(save_path):
            local_size = os.path.getsize(save_path)
            if local_size == total_size and total_size > 0:
                return f"Skipped (Exist): {file_name}"
        
        # 3. 开始下载 (GET 请求)
        response = session.get(url, stream=True, timeout=60)
        response.raise_for_status()
        
        block_size = 1024 * 1024 # 1MB
        # 设置 position 参数以避免多线程进度条错乱
        pbar = tqdm(total=total_size, unit='iB', unit_scale=True, desc=file_name, leave=False)
        
        with open(save_path, 'wb') as file:
            for data in response.iter_content(block_size):
                pbar.update(len(data))
                file.write(data)
        pbar.close()
        
        return f"✅ Success: {file_name}"

    except Exception as e:
        if os.path.exists(save_path):
            os.remove(save_path) # 删除下载失败的残留文件
        return f"❌ Failed: {file_name} | Error: {str(e)}"

def process_category(category, url):
    save_dir = os.path.join(BASE_SAVE_DIR, category)
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
        
    print(f"\n>>> 分析目录: {category} ...")
    
    session = get_session()
    try:
        # 获取文件列表
        response = session.get(url, timeout=30)
        soup = BeautifulSoup(response.text, 'html.parser')
        links = soup.find_all('a')
        
        # 筛选链接
        tasks = []
        for link in links:
            href = link.get('href')
            if href and is_target_file(href):
                full_url = urljoin(url, href)
                tasks.append(full_url)
                
        
        print(f"    需下载 {len(tasks)} 个文件。启动 {MAX_WORKERS} 线程下载...")
        
        # 多线程下载池
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            # 提交所有任务
            future_to_url = {executor.submit(download_single_file, u, save_dir): u for u in tasks}
            
            # 打印结果
            for future in as_completed(future_to_url):
                print(future.result())
                
    except Exception as e:
        print(f"!!! 获取目录失败: {e}")

def main():
    if not os.path.exists(BASE_SAVE_DIR):
        os.makedirs(BASE_SAVE_DIR)
        
    print(f"=== Ishii 高速下载器 (线程数: {MAX_WORKERS}) ===")
    
    for category, url in TARGETS.items():
        process_category(category, url)
        
    print("\n=== 全部完成 ===")

if __name__ == "__main__":
    main()
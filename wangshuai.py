import pandas as pd
import requests
import time

# ==========================================
# 1. 配置您的专属高德 API Key (Web服务)
# ==========================================
AMAP_API_KEY = "7c5c954346a4004cfbde3f016dbf8645"

def get_amap_coordinate(city_name):
    """
    调用高德地理编码 API 获取城市中心经纬度
    """
    url = "https://restapi.amap.com/v3/geocode/geo"
    params = {
        'address': city_name,
        'key': AMAP_API_KEY,
        'output': 'json'
    }
    
    try:
        response = requests.get(url, params=params, timeout=5)
        result = response.json()
        
        # 状态码 1 表示成功，且 count > 0 表示有匹配结果
        if result.get('status') == '1' and int(result.get('count', 0)) > 0:
            location = result['geocodes'][0]['location']
            lng, lat = location.split(',')
            return float(lng), float(lat)
        else:
            print(f"  [警告] 未找到坐标，API返回: {result.get('info')}")
            return None, None
            
    except Exception as e:
        print(f"  [错误] 请求发生异常: {e}")
        return None, None

def main():
    input_csv = "city_template.csv"
    output_csv = "city_coordinates_amap.csv"
    
    print(f"开始读取文件: {input_csv} ...")
    try:
        df = pd.read_csv(input_csv)
    except FileNotFoundError:
        print(f"找不到文件 {input_csv}，请确保它与此脚本在同一目录下。")
        return

    if "城市" not in df.columns:
        print("CSV 文件中找不到“城市”这一列，请检查表头名称。")
        return

    total_cities = len(df)
    print(f"共加载了 {total_cities} 个城市。开始向高德 API 请求坐标...\n")

    lngs = []
    lats = []

    for index, row in df.iterrows():
        # 清理可能存在的空格
        city_name = str(row['城市']).strip()
        
        # 优化查询词：部分单纯的名字（如"上海"）查询不如"上海市"精准
        query_name = city_name
        if not city_name.endswith(('市', '州', '区', '县', '盟', '地区')):
            query_name = city_name + '市'

        print(f"[{index + 1}/{total_cities}] 正在获取: {city_name} (查询词: {query_name})...", end="")
        
        lng, lat = get_amap_coordinate(query_name)
        
        if lng and lat:
            print(f" 成功 -> 经度: {lng}, 纬度: {lat}")
        else:
            # 如果加“市”失败（比如某些特殊自治州名称），尝试用原始名称再查一次
            if query_name != city_name:
                print(f"\n  尝试使用原名 '{city_name}' 重新查询...", end="")
                lng, lat = get_amap_coordinate(city_name)
                if lng and lat:
                    print(f" 成功 -> 经度: {lng}, 纬度: {lat}")
                else:
                    print(" 最终失败。")
            else:
                print(" 失败。")
                
        lngs.append(lng)
        lats.append(lat)
        
        # 【重要安全机制】
        # 高德个人开发者免费并发限制通常为 50次/秒。
        # 加上 0.1 秒的延迟，既能几分钟内跑完，又绝对安全不会被封禁或报错
        time.sleep(0.1) 

    # 将获取到的经纬度添加为新的两列
    df['经度(GCJ-02)'] = lngs
    df['纬度(GCJ-02)'] = lats

    print(f"\n处理完成！准备保存结果到: {output_csv} ...")
    # encoding='utf-8-sig' 确保生成的 csv 在 Excel 中直接打开不会中文乱码
    df.to_csv(output_csv, index=False, encoding='utf-8-sig')
    print("保存成功！您现在可以打开生成的 csv 文件查看结果。")

if __name__ == "__main__":
    main()
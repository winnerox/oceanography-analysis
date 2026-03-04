# 09_Utils - 工具函数

## 功能说明

本文件夹包含项目中使用的通用工具函数。

## 文件列表

| 文件名 | 功能 |
|--------|------|
| `draw_jitter_unit.m` | 绘制带抖动散点的箱线图 |

## 函数详解

### draw_jitter_unit

绘制单个带抖动散点的箱线图，用于展示数据分布。

**调用格式：**

```matlab
draw_jitter_unit(ax, center, data, w_box, w_scatter, sz_scatter, color_dots, color_box, alpha_s)
```

**参数说明：**

| 参数 | 类型 | 说明 |
|------|------|------|
| `ax` | handle | 坐标轴句柄 |
| `center` | scalar | X轴中心位置 |
| `data` | vector | 数据向量 |
| `w_box` | scalar | 箱体宽度 |
| `w_scatter` | scalar | 散点抖动范围 |
| `sz_scatter` | scalar | 散点大小 |
| `color_dots` | [R,G,B] | 散点颜色 |
| `color_box` | [R,G,B] | 箱体填充颜色 |
| `alpha_s` | scalar | 散点透明度 (0-1) |

**显示内容：**

- 箱体: IQR (25%-75%)
- 中线: 中位数 (红色)
- 点: 平均值 (黑色)
- 须线: 2.5%-97.5% 范围
- 散点: 显示离群值或采样点

## 使用示例

```matlab
figure;
ax = gca;
data = randn(1000, 1);
draw_jitter_unit(ax, 1, data, 0.5, 0.3, 10, [0.5 0.5 0.5], [0.95 0.95 1], 0.4);
```

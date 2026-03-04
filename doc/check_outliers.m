load('D:\work\Task_Convergence\Trend_Results\EN4_Average_Trends.mat', 'trend_TSLA');
EN4_data = trend_TSLA;
load('D:\work\Task_Convergence\Trend_Results\IAP_Average_Trends.mat', 'trend_TSLA');
IAP_data = trend_TSLA;
load('D:\work\Task_Convergence\Trend_Results\Ishii_Average_Trends.mat', 'trend_TSLA');
Ishii_data = trend_TSLA;

fprintf('===== 验证各数据集离群点规模 =====\n');
for n = 1:3
    e = EN4_data(:,:,n); e = e(~isnan(e(:)));
    i = IAP_data(:,:,n); i = i(~isnan(i(:)));
    s = Ishii_data(:,:,n); s = s(~isnan(s(:)));
    
    all_data = [e; i; s];
    low_w = prctile(all_data, 2.5);
    up_w = prctile(all_data, 97.5);
    
    e_out = sum(e < low_w | e > up_w);
    i_out = sum(i < low_w | i > up_w);
    s_out = sum(s < low_w | s > up_w);
    
    total = e_out + i_out + s_out;
    fprintf('\n=== 第 %d 阶 (T%d) ===\n', n, n);
    fprintf('全局 2.5%% ~ 97.5%%: [%.5f, %.5f]\n', low_w, up_w);
    fprintf('EN4   超出数量: %d (%.1f%%)\n', e_out, e_out/total*100);
    fprintf('IAP   超出数量: %d (%.1f%%)\n', i_out, i_out/total*100);
    fprintf('Ishii 超出数量: %d (%.1f%%)\n', s_out, s_out/total*100);
end

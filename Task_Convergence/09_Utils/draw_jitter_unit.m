function draw_jitter_unit(ax, center, data, w_box, w_scatter, sz_scatter, color_dots, color_box, alpha_s)
% draw_jitter_unit: Draws a single boxplot with jittered scatter points
% 
% Inputs:
%   ax: Axes handle
%   center: X-axis position center
%   data: Vector of data points
%   w_box: Width of the box
%   w_scatter: Width of the jitter area
%   sz_scatter: Size of scatter dots
%   color_dots: Color of scatter dots
%   color_box: Color of the box fill
%   alpha_s: Transparency of scatter dots

    data = data(~isnan(data)); 
    if isempty(data), return; end
    
    q1 = prctile(data, 25); q3 = prctile(data, 75);
    med_val = median(data); mean_val = mean(data);
    low_w = prctile(data, 2.5); up_w  = prctile(data, 97.5);
    
    % Identify Outliers for Scatter (Outside 2.5% - 97.5%)
    % Or simply plot all points? The user asked to "show outliers".
    % Let's plot ALL points but with transparency, or just outliers?
    % The user said "show outliers for different datasets". 
    % Standard boxplot convention: Plot outliers individually.
    % Here we plot ALL points to show distribution density ("Jitter Plot").
    % But to optimize performance, if N is huge, we might sample.
    
    if length(data) > 3000
        idx_sub = randperm(length(data), 3000);
        data_plot = data(idx_sub);
    else
        data_plot = data;
    end
    
    x_jit = center + (rand(size(data_plot)) - 0.5) * w_scatter;
    
    try 
        scatter(ax, x_jit, data_plot, sz_scatter, color_dots, 'filled', ...
                'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', alpha_s); 
    catch
    end
    
    % Draw Box components on top
    % Whiskers
    plot(ax, [center, center], [q1, low_w], 'k-', 'LineWidth', 1.2);
    plot(ax, [center, center], [q3, up_w],  'k-', 'LineWidth', 1.2);
    
    % Caps
    cap_w = w_box * 0.3;
    plot(ax, [center-cap_w, center+cap_w], [low_w, low_w], 'k-', 'LineWidth', 1.2);
    plot(ax, [center-cap_w, center+cap_w], [up_w, up_w],   'k-', 'LineWidth', 1.2);
    
    % Box
    x_L = center - w_box/2; x_R = center + w_box/2;
    patch(ax, [x_L, x_R, x_R, x_L], [q1, q1, q3, q3], color_box, ...
          'FaceAlpha', 0.8, 'EdgeColor', 'k', 'LineWidth', 1.2);
    
    % Median & Mean
    plot(ax, [x_L, x_R], [med_val, med_val], 'r-', 'LineWidth', 2);
    plot(ax, center, mean_val, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
end

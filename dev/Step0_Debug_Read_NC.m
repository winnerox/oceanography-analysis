%% Step0_Debug_Read_NC.m
% =========================================================================
% 数据诊断脚本 (Data Diagnostics)
% =========================================================================
% 功能：
%   直接读取 IAP 和 Ishii 的各一个样本文件。
%   打印：
%     1. 文件是否存在
%     2. 变量列表 (ncinfo)
%     3. 核心变量 (温度/盐度) 的维度、Min/Max/Mean
%   目的：
%     彻底搞清变量名（是 temp 还是 theta?）、维度顺序（Lon,Lat,Depth,Time?）
%     以及是否存在异常填充值 (FillValue)。
% =========================================================================

clear; clc;

% 定义我们要探测的目标文件 (根据之前用户提供的信息)
TargetFiles = { ...
    'IAP_Temp',  'D:\work\IAP_05_24\TEMP\IAPv4_Temp_monthly_1_6000m_year_2005_month_01.nc'; ...
    'IAP_Salt',  'D:\work\IAP_05_24\SALT\IAPv2_Salinity_monthly_1_6000m_year_2005_month_01.nc'; ...
    'Ishii_Temp','D:\work\Ishii_05_24\Temperature\temp.2005.nc'; ...
    'Ishii_Salt','D:\work\Ishii_05_24\Salinity\sal.2005.nc' ...
};

fprintf('>>> 开始数据诊断 <<<\n');

for k = 1:size(TargetFiles, 1)
    Label = TargetFiles{k, 1};
    FilePath = TargetFiles{k, 2};
    
    fprintf('\n=======================================================\n');
    fprintf('  检查文件 [%s]\n', Label);
    fprintf('  路径: %s\n', FilePath);
    fprintf('=======================================================\n');
    
    if ~exist(FilePath, 'file')
        fprintf('❌ 文件不存在！请检查路径。\n');
        % 尝试列出该目录下的文件，帮用户找找
        [Folder, ~, ~] = fileparts(FilePath);
        if exist(Folder, 'dir')
            fprintf('   -> 目录存在，列出前 5 个文件:\n');
            d = dir(fullfile(Folder, '*.nc'));
            for i = 1:min(5, length(d))
                fprintf('      %s\n', d(i).name);
            end
        else
            fprintf('   -> 连目录都不存在: %s\n', Folder);
        end
        continue;
    end
    
    % 1. 读取信息
    try
        Info = ncinfo(FilePath);
        fprintf('✅ 读取成功。包含变量:\n');
        VarNames = {Info.Variables.Name};
        fprintf('   %s ', VarNames{:}); fprintf('\n');
        
        % 2. 猜测核心变量
        TargetVar = '';
        if contains(Label, 'Temp')
            Possible = {'temperature', 'theta', 'temp', 't', 'pt'};
        else
            Possible = {'salinity', 'salt', 'sal', 's'};
        end
        
        for p = Possible
            if any(strcmpi(p{1}, VarNames))
                TargetVar = p{1};
                % 修正大小写
                TargetVar = VarNames{strcmpi(p{1}, VarNames)}; 
                break;
            end
        end
        
        if isempty(TargetVar)
            fprintf('⚠️ 没猜到核心变量，请人工核对上面的变量列表。\n');
            continue;
        end
        
        fprintf('-> 选中核心变量: [%s]\n', TargetVar);
        
        % 3. 读取少量数据进行统计
        % 为了快，只读第一帧 (Time=1)
        VarInfo = Info.Variables(strcmp(VarNames, TargetVar));
        Size = VarInfo.Size;
        fprintf('   维度: %s\n', mat2str(Size));
        
        Start = ones(1, length(Size));
        Count = Size;
        if length(Size) >= 4, Count(4) = 1; end % 只读第一个时间步
        if length(Size) >= 3 && Size(3) > 100, Count(3) = 1; end % 如果深度层太多，只读第一层(仅供测试)
        
        % 尝试完整读取第一帧 (不切深度，为了看数值范围)
        Count = Size; 
        if length(Size) >= 4, Count(4) = 1; end 
        
        Data = ncread(FilePath, TargetVar, Start, Count);
        Data = double(Data);
        
        % 统计
        TotalPoints = numel(Data);
        NaNPoints = sum(isnan(Data(:)));
        ValidData = Data(~isnan(Data));
        
        fprintf('   数据统计 (Time=1):\n');
        fprintf('     Total Points: %d\n', TotalPoints);
        fprintf('     NaN Points  : %d (%.1f%%)\n', NaNPoints, NaNPoints/TotalPoints*100);
        
        if ~isempty(ValidData)
            v_min = min(ValidData(:));
            v_max = max(ValidData(:));
            v_mean = mean(ValidData(:));
            v_med = median(ValidData(:));
            fprintf('     Min : %.4e\n', v_min);
            fprintf('     Max : %.4e\n', v_max);
            fprintf('     Mean: %.4e\n', v_mean);
            fprintf('     Med : %.4e\n', v_med);
            
            % 检查 FillValue
            if v_max > 1000 || v_min < -100
                fprintf('⚠️ 警告: 发现异常大/小的值，可能是 FillValue 未处理！\n');
            else
                fprintf('✅ 数值范围看起来正常 (海水温盐范围)。\n');
            end
        else
            fprintf('⚠️ 所有数据都是 NaN？\n');
        end
        
    catch ME
        fprintf('❌ 读取或分析时出错: %s\n', ME.message);
    end
    
end

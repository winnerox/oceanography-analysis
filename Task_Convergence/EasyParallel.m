function EasyParallel_Pro()
    % --- 界面初始化 ---
    fig = uifigure('Name', 'MATLAB 并行专家版 (R2023b)', 'Position', [100 100 850 500]);
    
    state.workDir = 'D:\work'; 
    state.allScripts = string.empty; 
    
    % --- UI 布局 ---
    uilabel(fig, 'Position', [20 460 70 22], 'Text', '主目录:', 'FontWeight', 'bold');
    pathEdit = uieditfield(fig, 'text', 'Position', [90 460 550 22], 'Value', state.workDir, 'Editable', 'off');
    uibutton(fig, 'push', 'Text', '选择文件夹', 'Position', [650 460 100 22], 'ButtonPushedFcn', @(s,e) selectFolder());

    uilabel(fig, 'Position', [20 425 200 22], 'Text', '待选脚本:', 'FontWeight', 'bold');
    scriptList = uilistbox(fig, 'Position', [20 80 350 340], 'Multiselect', 'on');
    
    uilabel(fig, 'Position', [390 425 100 22], 'Text', '运行日志:', 'FontWeight', 'bold');
    logArea = uitextarea(fig, 'Position', [390 80 440 340], 'Editable', 'off', 'Value', ""); 

    runBtn = uibutton(fig, 'push', 'Text', '🚀 开始并行运行', 'Position', [20 20 350 45], ...
        'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', [0.4 0.8 1]);
    
    statusLamp = uilamp(fig, 'Position', [380 32 20 20], 'Color', 'green');
    progressLabel = uilabel(fig, 'Position', [410 32 300 22], 'Text', '就绪');

    refreshList();

    % ================= 核心逻辑 =================

    function selectFolder()
        selPath = uigetdir(state.workDir, '选择根目录');
        if selPath ~= 0, state.workDir = selPath; pathEdit.Value = selPath; refreshList(); end
    end

    function refreshList()
        if ~exist(state.workDir, 'dir'), return; end
        files = dir(fullfile(state.workDir, '**', '*.m'));
        displayNames = {}; fullPaths = [];
        for i = 1:length(files)
            if startsWith(files(i).name, '.') || contains(files(i).folder, 'private'), continue; end
            relFolder = strrep(files(i).folder, state.workDir, '');
            if isempty(relFolder), dName = files(i).name;
            else
                if startsWith(relFolder, filesep), relFolder = relFolder(2:end); end
                dName = fullfile(relFolder, files(i).name);
            end
            displayNames{end+1} = dName; fullPaths(end+1) = string(fullfile(files(i).folder, files(i).name));
        end
        scriptList.Items = displayNames; state.allScripts = fullPaths; 
    end

    runBtn.ButtonPushedFcn = @(src, event) startParallelTasks();

    function startParallelTasks()
        try
            selectedItems = scriptList.Value;
            if isempty(selectedItems) || isequal(selectedItems, {'无脚本文件'}), return; end
            [~, selectedIdx] = intersect(scriptList.Items, selectedItems, 'stable');
            
            logArea.Value = "▶️ 正在初始化任务..."; 
            statusLamp.Color = 'yellow';
            runBtn.Enable = 'off';
            
            pool = backgroundPool;
            total = length(selectedIdx);
            completed = 0;
            
            % 修复点：创建一个 futures 数组来追踪所有任务
            futures = parallel.Future.empty(0, total);

            for i = 1:total
                realIdx = selectedIdx(i);
                % 强制类型转换，防止 fileparts 报错
                scriptFullPath = char(state.allScripts(realIdx)); 
                [scriptDir, scriptName, ~] = fileparts(scriptFullPath);
                
                % 提交任务
                futures(i) = parfeval(pool, @run_wrapper, 0, scriptFullPath, scriptDir);
                
                % 设置回调
                afterEach(futures(i), @(~) taskDone(scriptName, "✅ 成功"), 0);
                afterEach(futures(i), @(err) taskDone(scriptName, "❌ 失败: " + err.message), 1);
            end
            
            % 修复点：监听整个 futures 数组
            afterAll(futures, @(~) finalize(), 0);

        catch ME
            % --- 终极捕获：如果循环内报错，直接显示在日志里 ---
            logArea.Value = [logArea.Value; "‼️ 界面线程崩溃!"; "原因: " + ME.message];
            statusLamp.Color = 'red';
            runBtn.Enable = 'on';
        end
        
        function taskDone(name, status)
            completed = completed + 1;
            msg = "[" + datestr(now, 'HH:MM:SS') + "] " + name + ": " + status;
            logArea.Value = [logArea.Value; string(msg)]; 
            scroll(logArea, 'bottom');
            progressLabel.Text = sprintf('正在处理 %d/%d...', completed, total);
        end

        function finalize()
            statusLamp.Color = 'green';
            runBtn.Enable = 'on';
            progressLabel.Text = '所有任务已完成！';
            logArea.Value = [logArea.Value; "🎉 运行结束。"];
        end
    end
end

function run_wrapper(fullPath, scriptDir)
    addpath(genpath(char(scriptDir)));
    cd(char(scriptDir));
    run(char(fullPath));
end
try
    m = matfile('D:/work/EN4_TSLA_Terms/EN4_Mean_State.mat');
    disp('Variables in EN4_Mean_State.mat:');
    disp(whos(m));
catch ME
    disp(ME.message);
end

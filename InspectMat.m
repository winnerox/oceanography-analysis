try
    m = matfile('D:/work/MAT_Data/EN4_TSLA_Terms_1to8_Average.mat');
    disp('Variables in EN4_TSLA_Terms_1to8_Average.mat:');
    disp(whos(m));
catch ME
    disp(ME.message);
end

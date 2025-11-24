function [s_class,coe,correct_flag] = checkparam(s_class,coe,Kotf,paramfiles)
% 检查并处理s_class参数的函数
% 输入: s_class - 包含kAmean, phaseA, modFac的结构体
% 输出: s_class - 可能更新后的结构体

% 加载已保存的数据
if exist(paramfiles, 'file')
    load(paramfiles);
else
    saved_kAmean = [];
    saved_phaseA = [];
    saved_modFac_sum = 0;
end

% 计算当前kAmean范数的平均值
current_norms = zeros(size(s_class.kAmean, 1), 1);
for i = 1:size(s_class.kAmean, 1)
    if i ==1||i ==5
        current_norms(i) = norm(s_class.kAmean(i, :))*2;
    elseif i ==2||i ==3
        current_norms(i) = norm(s_class.kAmean(i, :));
    else
        current_norms(i) = norm(s_class.kAmean(i, :))*sqrt(2);
    end
end
mean_norm = mean(current_norms);

% 检查是否有偏离平均值5%以上的范数
% max_deviation = max(abs(current_norms - mean_norm) / mean_norm);
max_deviation = max(abs(current_norms - mean_norm) );
if max_deviation > 8
    correct_flag=0;
    % 有偏离超过5%的情况，加载保存的数据
    if ~isempty(saved_kAmean)
        s_class.kAmean = saved_kAmean;
        s_class.phaseA = saved_phaseA;
        fprintf('已加载保存的kAmean和phaseA数据\n');
    else
        fprintf('检测到异常但没有保存的数据可加载\n');
    end
    s_class = PCMseparateF_3D(s_class);
    s_class.modFac = ModulationFactorF(s_class);
else
    correct_flag=1;
    % 准备保存当前数据
    current_modFac_sum = sum(s_class.modFac);

    % 判断是否需要保存
    if current_modFac_sum > saved_modFac_sum
        saved_kAmean = s_class.kAmean;
        saved_phaseA = s_class.phaseA;
        saved_modFac_sum = current_modFac_sum;
        save(paramfiles, 'saved_kAmean', 'saved_phaseA', 'saved_modFac_sum');
        fprintf('已保存新的kAmean、phaseA和modFac_sum数据\n');
    else
        fprintf('modFac求和未超过保存值，不保存数据\n');
    end
end
current_norms = zeros(size(s_class.kAmean, 1), 1);
for i = 1:size(s_class.kAmean, 1)
    if i ==1||i ==5
        current_norms(i) = norm(s_class.kAmean(i, :))*2;
    elseif i ==2||i ==3
        current_norms(i) = norm(s_class.kAmean(i, :));
    else
        current_norms(i) = norm(s_class.kAmean(i, :))*sqrt(2);
    end
end
mean_norm = mean(current_norms);
coe=mean_norm./Kotf;
s_class.kAmean=gpuArray(single(s_class.kAmean));
s_class.phaseA=gpuArray(single(s_class.phaseA));
end
% 获取所有变量信息
vars = whos;

% 筛选GPU变量
gpuVars = vars(strcmp({vars.class}, 'gpuArray'));

% 按显存占用排序（降序）
[~, idx] = sort([gpuVars.bytes], 'descend');
sortedGpuVars = gpuVars(idx);

% 打印排序结果
if ~isempty(sortedGpuVars)
    fprintf('\n%-20s %-15s %-15s %s\n', '变量名', '大小', '字节', '类型');
    for i = 1:length(sortedGpuVars)
        v = sortedGpuVars(i);
        sizeStr = sprintf('%dx', v.size);
        sizeStr = sizeStr(1:end-1);
        fprintf('%-20s %-15s %-15d %s\n', v.name, sizeStr, v.bytes, v.class);
    end
else
    fprintf('\n无GPU变量\n');
end
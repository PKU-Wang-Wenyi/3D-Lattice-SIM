% 检查内存变量
vars = whos;
[~, idx] = sort([vars.bytes], 'descend');
sortedVars = vars(idx);

disp('内存变量占用排序:');
fprintf('%-20s %-15s %-15s %s\n', '变量名', '大小', '字节', '类型');
for i = 1:length(sortedVars)
    v = sortedVars(i);
    sizeStr = strjoin(cellstr(num2str(v.size')), 'x');
    fprintf('%-20s %-15s %-15d %s\n', v.name, sizeStr, v.bytes, v.class);
end

% 检查显存变量
gpuVars = vars(strcmp({vars.class}, 'gpuArray'));
[~, gpuIdx] = sort([gpuVars.bytes], 'descend');
sortedGpuVars = gpuVars(gpuIdx);

disp('显存变量占用排序:');
if ~isempty(sortedGpuVars)
    fprintf('%-20s %-15s %-15s %s\n', '变量名', '大小', '字节', '类型');
    for i = 1:length(sortedGpuVars)
        v = sortedGpuVars(i);
        sizeStr = strjoin(cellstr(num2str(v.size')), 'x');
        fprintf('%-20s %-15s %-15d %s\n', v.name, sizeStr, v.bytes, v.class);
    end
else
    disp('无GPU变量');
end
function [rawdata, cut_ind] = rawdata_slice_cut_fun(rawdata, maxpos, rawdata_slice_cut, cut_ind, maxint)
% 截取使maxint和最大的连续片段
% 参数：
%   rawdata: 输入4D数组 (height x width x channels x frames)
%   maxpos: 原始最大值位置（索引）
%   rawdata_slice_cut: 需截取的总帧数（L）
%   cut_ind: 预设截取索引（若不为0则直接使用）
%   maxint: 各帧光强值向量（长度 = size(rawdata,4)）
% 返回：
%   rawdata: 截取后的4D数组
%   cut_ind: 实际截取的帧索引范围
m = size(rawdata, 4); % 总帧数
L = rawdata_slice_cut;
if cut_ind(1) ~= 0
    rawdata = rawdata(:, :, :, cut_ind);
    return;
elseif L <= 0 || L > m
    cut_ind = 1:m; % 无效L时截取全部
    return;
end
% ---- 核心修改：搜索maxint和最大的连续片段 ----
% 1. 确定搜索范围（以maxpos为中心，扩展±L）
searchStart = max(1, maxpos - L);
searchEnd = min(m, maxpos + L);
searchWindow = searchStart:searchEnd;

% 2. 计算搜索窗口内所有连续L长度子数组的maxint和
sums = zeros(1, numel(searchWindow) - L + 1);
for i = 1:(numel(searchWindow) - L + 1)
    idx = searchWindow(i:i + L - 1);
    sums(i) = sum(maxint(idx)); % 关键目标：最大化光强和
end

% 3. 找到和最大的子数组起始位置
[~, maxSumIdx] = max(sums);
startIdx = searchWindow(maxSumIdx);
endIdx = startIdx + L - 1;

% 4. 边界保护（理论上已通过searchWindow保证，此处冗余校验）
startIdx = max(1, startIdx);
endIdx = min(m, endIdx);

% ---- 执行截取 ----
rawdata = rawdata(:, :, :, startIdx:endIdx);
cut_ind = startIdx:endIdx;
end
% function [rawdata,cut_ind] = rawdata_slice_cut_fun(rawdata, maxpos, rawdata_slice_cut,cut_ind,maxint)
% % 截取数组中以最大值位置为中心的指定长度片段，处理边界情况
% % 参数：
% %   arr：输入数组
% %   maxpos：最大值所在位置（MATLAB索引从1开始）
% %   rawdata_slice_cut：需截取的总帧数（前后共多少帧）
% % 返回：
% %   截取后的数组片段
% if cut_ind(1)==0
%     m = size(rawdata,4);
%     L = rawdata_slice_cut;
%     % 处理无效的截取长度（如负数）
%     if L <= 0
%         return;
%     end
%
%     % 计算理想的起始和结束索引（MATLAB索引从1开始）
%     start_ideal = maxpos - floor((L - 1) / 2);
%     end_ideal = start_ideal + L - 1;
%
%     % 检查理想范围是否在数组边界内
%     if start_ideal >= 1 && end_ideal <= m
%         % 未超出边界，直接截取
%         rawdata = rawdata(:,:,:,start_ideal:end_ideal);
%         cut_ind=start_ideal:end_ideal;
%     else
%         % 超出边界，截取最后L个帧（若数组长度不足则截取全部）
%         start = max(1, m - L + 1);
%         rawdata = rawdata(:,:,:,start:min(m, start + L - 1));
%         cut_ind=start:min(m, start + L - 1);
%     end
% else
%     rawdata = rawdata(:,:,:,cut_ind);
% end
% end
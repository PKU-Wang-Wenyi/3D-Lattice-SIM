function f=imreadstack_ZT(info,nT)
numimage=info(1).nPhase*info(1).sizeZ;
start_image=(nT-1)*numimage+1;
end_image=(nT)*numimage;
target_frames = start_image:end_image;
frameCounts = info(1).frameCounts(:)';
cum_end = cumsum(frameCounts);
cum_start = [1, cum_end(1:end-1)+1];
file_idx = discretize(target_frames, [cum_start, inf]);
local_idx = target_frames - cum_start(file_idx) + 1;
[~, base_name, ext] = fileparts(info(1).name);
base_prefix = extractBefore(base_name, '.ome');
file_cache = containers.Map('KeyType','char','ValueType','any');
info_available = containers.Map('KeyType', 'char', 'ValueType', 'logical');
f=zeros(info(1).sizeX,info(1).sizeY,numimage,'uint16');
max_retry = 5;
retry_wait = 2;
for i = 1:numimage
    file_i = file_idx(i);
    if file_i == 1
        filename = fullfile(info(1).folder, [base_prefix '.ome' ext]);
    else
        filename = fullfile(info(1).folder, [base_prefix '_' num2str(file_i - 1) '.ome' ext]);
    end
    % if ~isKey(file_cache, filename)
    %     file_cache(filename) = imfinfo(filename);
    % end
    % file_info = file_cache(filename);
    if ~isKey(info_available, filename)
        try
            file_cache(filename) = imfinfo(filename);
            info_available(filename) = true;
        catch
            warning('无法获取 %s 的 imfinfo，使用非加速读取方式。', filename);
            info_available(filename) = false;
        end
    end
    success = false;
    attempt = 0;
    while ~success && attempt < max_retry
        try
            % f(:,:,i) = imread(filename, local_idx(i), 'Info', file_info);
            % success = true;
            if info_available(filename)
                f(:, :, i) = imread(filename, local_idx(i), 'Info', file_cache(filename));
            else
                f(:, :, i) = imread(filename, local_idx(i));
            end
            success = true;
        catch ME
            attempt = attempt + 1;
            if attempt >= max_retry
                error('读取图像 %s 的第 %d 张帧失败（重试 %d 次后仍未成功）。\n错误信息：%s', ...
                    filename, local_idx(i), attempt, ME.message);
            else
                warning('读取图像 %s 的第 %d 张帧失败，等待 %d 秒后重试（已尝试 %d 次）...', ...
                    filename, local_idx(i), retry_wait, attempt);
                pause(retry_wait);
            end
        end
    end
end
end
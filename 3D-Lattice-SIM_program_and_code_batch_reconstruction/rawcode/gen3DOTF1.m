function [OTF,param] = gen3DOTF_cache(obj, nz)
    % ===============================
    % 目标: 自动判断是否需要重新生成 OTF
    % 若参数一致 -> 直接加载已有 OTF
    % 若参数不一致 -> 重新计算并保存
    % ===============================

    % === 构造参数结构体 ===
    param.sizeX = obj.w;
    param.sizeY = obj.w;
    param.sizeZ = nz;
    param.pixelsize = obj.pixelsize; % in nm
    param.Zstep = obj.Zstep;         % in nm
    param.lambda = obj.lambda;       % in nm
    param.exlambda = obj.exlambda;   % in nm
    param.imgsize = [param.sizeX, param.sizeY, param.sizeZ];
    param.medium_refractive_index = 1.518;
    param.Sample_refractive_index = 1.518;
    param.coverglass_refractive_index = 1.518;
    param.obj_workingdistence = 0;
    param.nz = nz;
    param.w = obj.w;
    param.NA = obj.NA;

    % === 定义缓存文件路径 ===
    cacheFile = fullfile(pwd, 'OTF_cache.mat');

    % === 判断是否存在缓存文件 ===
    if isfile(cacheFile)
        try
            S = load(cacheFile, 'param', 'OTF');
            oldParam = S.param;
            % === 对比参数结构体 ===
            if isequaln(param, oldParam)
                fprintf('[gen3DOTF_cache] 参数一致，直接加载 OTF_cache.mat 中的 OTF。\n');
                OTF = S.OTF;
                return;
            else
                fprintf('[gen3DOTF_cache] 参数不同，重新生成 OTF。\n');
            end
        catch ME
            % 使用 ME.identifier（若没有则 fallback），并显示 message
            id = ME.identifier;
            if isempty(id)
                id = 'gen3DOTF_cache:LoadFailed';
            end
            warning(id, '[gen3DOTF_cache] 读取缓存失败: %s， 将重新生成 OTF。', ME.message);
            % 可选：打印简短报表以便调试（取消注释可查看堆栈）
            % fprintf(2, '%s\n', getReport(ME,'basic'));
        end
    else
        fprintf('[gen3DOTF_cache] 未检测到 OTF_cache.mat，开始生成。\n');
    end

    % === 调用原始 OTF 生成函数 ===
    [OTF, ~] = get_modelOTF_fast(param);
    OTF = gather(OTF);

    % === 保存缓存 ===
    try
        save(cacheFile, 'OTF', 'param', '-v7.3');
        fprintf('[gen3DOTF_cache] 已保存新的 OTF_cache.mat。\n');
    catch ME
        id = ME.identifier;
        if isempty(id)
            id = 'gen3DOTF_cache:SaveFailed';
        end
        % 报警：保存失败，但仍把 OTF 返回给调用者
        warning(id, '[gen3DOTF_cache] 保存 OTF 缓存失败: %s', ME.message);
        % 可选：打印扩展错误信息用于调试
        % fprintf(2, '%s\n', getReport(ME,'basic'));
    end
end

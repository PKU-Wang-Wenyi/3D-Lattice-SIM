function uniqueFiles = getUniqueTiffFiles(folderPath)
if nargin < 1
    folderPath = pwd;
end
filePattern = fullfile(folderPath,'/**/*ome.tif');
fileList = dir(filePattern);
totalFiles = length(fileList);
if totalFiles == 0
    uniqueFiles = [];
    disp('未找到任何tif文件');
    return;
end
% pagingPattern = '^(.*)_([0-9]+)\.ome\.tif$';
% mainPattern   = '^(.*)\.ome\.tif$';
pagingPattern = '^(.*)_(\d+)(\.ome\.tif)$';
mainPattern = '^(.*)(\.ome\.tif)$';
baseInfo = struct(...
    'coreName', cell(totalFiles, 1), ...
    'version', zeros(totalFiles, 1), ...
    'origIdx', 1:totalFiles, ...
    'folder', cell(totalFiles, 1), ...
    'SIZEX', zeros(totalFiles, 1), ...
    'SIZEY', zeros(totalFiles, 1), ...
    'Imglength', zeros(totalFiles, 1) ...
    );
for i = 1:totalFiles
    fullName = fileList(i).name;
    baseInfo(i).folder = fileList(i).folder;
    fullPath = fullfile(fileList(i).folder, fullName);

    try
        info = imfinfo(fullPath);
        baseInfo(i).SIZEX = info(1).Width;
        baseInfo(i).SIZEY = info(1).Height;
        baseInfo(i).Imglength = length(info);
    catch
        warning('读取 %s 的 imfinfo 失败，跳过尺寸记录。', fullPath);
    end

    tokens = regexp(fullName, pagingPattern, 'tokens');
    if ~isempty(tokens)
        coreNameWithoutVersion = tokens{1}{1};
        pageStr = tokens{1}{2};
        suffix = tokens{1}{3};
        coreName = [coreNameWithoutVersion suffix]; % 恢复完整coreName，带.ome.tif
        baseInfo(i).coreName = coreName;
        baseInfo(i).version = str2double(pageStr);
    else
        tokensMain = regexp(fullName, mainPattern, 'tokens');
        if ~isempty(tokensMain)
            coreName = [tokensMain{1}{1} tokensMain{1}{2}];
            baseInfo(i).coreName = coreName;
            baseInfo(i).version = 0;
        else
            warning('无法解析文件名: %s', fullName);
            baseInfo(i).coreName = fullName;
            baseInfo(i).version = 0;
        end
    end
end
uniqueCoreNames = unique({baseInfo.coreName});
numUnique = length(uniqueCoreNames);
uniqueFiles = struct(...
    'name', cell(numUnique, 1), ...
    'folder', repmat({folderPath}, numUnique, 1), ...
    'count', zeros(numUnique, 1), ...
    'versions', cell(numUnique, 1), ...
    'frameCounts', cell(numUnique, 1), ...
    'sizeX', zeros(1, 1), ...
    'sizeY', zeros(1, 1), ...
    'sizeZ', zeros(1, 1), ...
    'sizeT', zeros(1, 1), ...
    'Zstep', zeros(1, 1) ...
    );
for i = 1:numUnique
    currentCore = uniqueCoreNames{i};
    isCurrent = strcmp({baseInfo.coreName}, currentCore);
    relatedFiles = baseInfo(isCurrent);
    uniqueFiles(i).name = currentCore;
    uniqueFiles(i).count = length(relatedFiles);
    uniqueFiles(i).folder = relatedFiles(1).folder;
    uniqueFiles(i).sizeX = relatedFiles(1).SIZEX;
    uniqueFiles(i).sizeY = relatedFiles(1).SIZEY;
    uniqueFiles(i).versions = [relatedFiles.version];
    uniqueFiles(i).frameCounts = [relatedFiles.Imglength];
    %%
    paramsFile = fullfile(uniqueFiles(i).folder, 'params.txt');
    if exist(paramsFile, 'file')
        fileContent = fileread(paramsFile);
        ntPattern = 'nt\s*=\s*(\d+)\s*;';
        ntMatch = regexp(fileContent, ntPattern, 'tokens');
        if ~isempty(ntMatch)
            uniqueFiles(i).sizeT = str2double(ntMatch{1}{1});
        else
            warning(['无法从文件 ' paramsFile ' 中提取nt值']);
            uniqueFiles(i).sizeT = 1;
        end
        nzPattern = 'nz\s*=\s*(\d+)\s*;';
        nzMatch = regexp(fileContent, nzPattern, 'tokens');
        if ~isempty(nzMatch)
            uniqueFiles(i).sizeZ = str2double(nzMatch{1}{1});
        else
            warning(['无法从文件 ' paramsFile ' 中提取nz值']);
            uniqueFiles(i).sizeZ = 1;
        end
        dzPattern = 'dz\s*=\s*([\d.]+)\s*;';
        dzMatch = regexp(fileContent, dzPattern, 'tokens');
        if ~isempty(dzMatch)
            uniqueFiles(i).Zstep = str2double(dzMatch{1}{1})*1000;
        else
            warning(['无法从文件 ' paramsFile ' 中提取dz值']);
            uniqueFiles(i).Zstep = 125;
        end
    else
        try
            imgname=fullfile(uniqueFiles(i).folder,uniqueFiles(i).name);
            addpath(genpath('.\bfmatlab\'));
            reader = bfGetReader(imgname);
            omeMeta = reader.getMetadataStore();
            omeXML = char(omeMeta.dumpXML());
            expr = '<Plane[^>]*DeltaT="([\d\.]+)"[^>]*TheC="(\d+)"[^>]*TheT="(\d+)"[^>]*TheZ="(\d+)"';
            tokens = regexp(omeXML, expr, 'tokens');
            planeData = cellfun(@(x) str2double(x), tokens, 'UniformOutput', false);
            planeData = vertcat(planeData{:});
            [~, idx] = max(planeData(:,1));
            maxPlane = planeData(idx,:);
            fprintf('解析XML文件成功！最大DeltaT = %.2f ms, TheC = %d, TheT = %d, TheZ = %d\n', ...
                maxPlane(1), maxPlane(2), maxPlane(3), maxPlane(4));
            uniqueFiles(i).Zstep =0;
            uniqueFiles(i).sizeZ =maxPlane(4)+1;
            uniqueFiles(i).sizeT =maxPlane(2)+1;
            if uniqueFiles(i).sizeT==1&&uniqueFiles(i).sizeT*uniqueFiles(i).sizeZ*13>idx
                uniqueFiles(i).sizeZ=uniqueFiles(i).sizeZ-1;
            elseif uniqueFiles(i).sizeT>1&&uniqueFiles(i).sizeT*uniqueFiles(i).sizeZ*13>idx
                uniqueFiles(i).sizeT=uniqueFiles(i).sizeT-1;
            end
        catch
            disp([imgname, ' 数据处理不了！']);
            % continue
        end
    end
end
end
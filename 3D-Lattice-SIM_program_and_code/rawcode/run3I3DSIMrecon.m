function run3I3DSIMrecon(app)
basefolder=app.DataFolderEditField.Value;
tiffile = getUniqueTiffFiles(basefolder);
containsTubulin = contains({tiffile.name}, 'npc', 'IgnoreCase', true);
containsBeads = contains({tiffile.name}, 'beads', 'IgnoreCase', true);
containsOther1 = contains({tiffile.name}, 'neuron', 'IgnoreCase', true);
containsOther2= contains({tiffile.name}, 'act', 'IgnoreCase', true);
priorityFiles = containsTubulin | containsBeads |containsOther1|containsOther2;
% priorityFiles = containsTubulin | containsBeads |containsOther2;
tiffile = [tiffile(priorityFiles); tiffile(~priorityFiles)];
parallel.gpu.enableCUDAForwardCompatibility(true)
%% params
param.NA=1.49;
param.pixelsize=65;% in nm
param.RefractiveIndex=1.518;
param.nAngle=1;
param.nPhase=13;
param.Zstep=125;
param.slice_cut=0;
param.attStrength=app.attStrengthSpinner.Value;
param.attWidth=app.attWidthSpinner.Value;
param.w1=app.w1Spinner.Value;
param.w2=app.w2Spinner.Value;
param.RL=app.RLSpinner.Value;
param.Dark=app.DarkSpinner.Value;
%%
% max_retry = 1;
% retry_wait = 1;
for i=1:length(tiffile)
    % success = false;
    % attempt = 0;
    % tiffile(i).sizeT=2;
    % if tiffile(i).sizeT+tiffile(i).sizeZ>0
    fullfile(tiffile(i).folder,tiffile(i).name)
    param.imgname=fullfile(tiffile(i).folder,tiffile(i).name);
    param.rawinfo=tiffile(i);
    if contains(param.imgname,'Channel_')
        if contains(param.imgname,'Channel_561')
            param.exlambda=561;
            param.lambda=610;
        elseif contains(param.imgname,'Channel_488')
            param.exlambda=488;
            param.lambda=525;
        elseif contains(param.imgname,'Channel_640')
            param.exlambda=640;
            param.lambda=685;
        end
    else
        if contains(param.imgname,'561')
            param.exlambda=561;
            param.lambda=610;
        elseif contains(param.imgname,'488')
            param.exlambda=488;
            param.lambda=525;
        elseif contains(param.imgname,'640')
            param.exlambda=640;
            param.lambda=685;
        end
    end
    if tiffile(i).Zstep>0
        param.Zstep=tiffile(i).Zstep;
    end
    % try
    %     [SIM,WF]=recon3I_SIM_ZT(param);
    % catch
    %     disp(['没有处理！无法打开！',fullfile(tiffile(i).folder,tiffile(i).name)])
    %     continue
    % end
    [path,name]=fileparts(param.imgname);
    param.resultname=fullfile(path,[name,'TEMP_STR.tif']);
    % while ~success && attempt < max_retry
    % try
    [SIM, WF] = recon3I_SIM_ZT(param);
    % success = true;
    % catch ME
    %     attempt = attempt + 1;
    %     % if attempt >= max_retry
    %     %     disp(['失败：重试 ', num2str(attempt), ' 次仍无法处理！文件：', ...
    %     %         fullfile(tiffile(i).folder, tiffile(i).name)])
    %     %     disp(['错误信息：', ME.message])
    %     %     continue;
    %     % else
    %     %     disp(['处理失败，第 ', num2str(attempt), ' 次尝试后等待 ', ...
    %     %         num2str(retry_wait), ' 秒... 文件：', fullfile(tiffile(i).folder, tiffile(i).name)])
    %     %     pause(retry_wait);
    %     % end
    %     continue;
    % end
    % end
    % SIM=mat2gray(SIM).*65535;
    % imwritestack_16(SIM,fullfile(path,[name,'SIM.tif']))
    % imwritestack_16(WF,fullfile(path,[name,'WF.tif']))
end
end
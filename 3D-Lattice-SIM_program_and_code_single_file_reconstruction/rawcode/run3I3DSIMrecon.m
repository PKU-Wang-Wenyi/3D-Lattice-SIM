function run3I3DSIMrecon(app)
parallel.gpu.enableCUDAForwardCompatibility(true)
%% params
param.imgname =app.DataEditField.Value;
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
param.sizeT=app.sizeTSpinner.Value;
param.sizeZ=app.sizeZSpinner.Value;
%%
% max_retry = 1;
% retry_wait = 1;
param.Zstep=125;
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
[path,name]=fileparts(param.imgname);
param.resultname=fullfile(path,[name,'TEMP_STR.tif']);
[~, ~] = recon3I_SIM_ZT_single_data(param);
end
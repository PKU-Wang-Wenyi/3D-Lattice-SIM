%%%
clear;close all
parallel.gpu.enableCUDAForwardCompatibility(true)
if exist('reconinfo.txt','file')
    delete('reconinfo.txt')
end
diary('reconinfo.txt')
param.NA=1.49;
param.pixelsize=65;% in nm
param.imgname = '..\lattice_hex_488_EMTB_mBaojin-P_012\lattice_hex_488+561_COS7_bj-EM—BD—SEC_013_Channel_488_MMStack_Default.ome.tif';
%%
param.RefractiveIndex=1.518;
param.attStrength=0.90;
param.attWidth=1.5;
param.w1=1;
param.w2=0.9;
param.RL=0;
param.nAngle=1;
param.nPhase=13;
param.Zstep=125;
param.Dark=0;
param.slice_cut=0;
%%
param.sizeT=1;
param.sizeZ=6;
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
param.resultname=fullfile(path,'TEMP_STR.tif');
[SIM, WF] = recon3I_SIM_ZT_single_data(param);
diary off
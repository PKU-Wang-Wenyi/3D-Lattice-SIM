clear;close all
addpath("bfmatlab\")
parallel.gpu.enableCUDAForwardCompatibility(true)
basefolder='..\rawdata\';%rawdatafolder
% basefolder='P:\3D3ISIM\ZeissSquireLattice3DSIM\rawdata\';%rawdatafolder
rawfile=dir(fullfile(basefolder,'**/*.czi'));
rawfile = rawfile(~endsWith({rawfile.name}, 'SIM.czi'));
rawfile = rawfile(~endsWith({rawfile.name}, 'SIM².czi'));
param.NA=1.49;
param.pixelsize=65;
param.attStrength=0.9;
param.attWidth=2;
param.w1=1;
param.w2=0.1;
param.RL=0;
param.nAngle=1;
param.nPhase=13;
param.Dark=0;
for i=1:length(rawfile)
    rawfilename=fullfile(rawfile(i).folder,rawfile(i).name);
    [rawdata,param]=getczidata(rawfilename,param);
    savefolder=rawfilename(1:end-4);
    if ~exist(savefolder,"dir")
        mkdir(savefolder)
    end
    param.WFsavename=fullfile(savefolder,'WF.tif');
    param.SIMsavename=fullfile(savefolder,'SIM.tif');
    reconZeissSIM(rawdata,param);
end

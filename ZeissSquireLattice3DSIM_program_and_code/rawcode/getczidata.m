function [rawdata,param]=getczidata(imgname,param)
info=czifinfo(imgname);
% imgnum=length(info.sBlockList_P0);
textData=info.metadataXML;
sizeX = regexp(textData, '<SizeX>(\d+)</SizeX>', 'tokens');
sizeY = regexp(textData, '<SizeY>(\d+)</SizeY>', 'tokens');
sizeZ = regexp(textData, '<SizeZ>(\d+)</SizeZ>', 'tokens');
sizeH = regexp(textData, '<SizeH>(\d+)</SizeH>', 'tokens');
if isempty(sizeZ)
    sizeZ = 1;
else
    sizeZ = str2double(sizeZ{1}{1});
end
param.sizeZ=sizeZ;
sizeX = str2double(sizeX{1}{1});
sizeY = str2double(sizeY{1}{1});
sizeH = str2double(sizeH{1}{1});
exlambda = regexp(textData, '<ExcitationWavelength>(\d+\.\d+|\d+)</ExcitationWavelength>', 'tokens');
lambda = regexp(textData, '<EmissionWavelength>(\d+\.\d+|\d+)</EmissionWavelength>', 'tokens');
param.exlambda = str2double(exlambda{1}{1});
param.lambda = str2double(lambda{1}{1});
Zstep = regexp(textData, '<ScalingZ>(\d+\.\d+|\d+)e-007</ScalingZ>', 'tokens');
param.Zstep = str2double(Zstep{1}{1}).*100;
RefractiveIndex= regexp(textData, '<RefractiveIndex>(\d+\.\d+|\d+)</RefractiveIndex>', 'tokens');
param.RefractiveIndex=str2double(RefractiveIndex{1}{1});
rawdata=bfopen(imgname);
rawdata=rawdata{1};
rawdata = double(cell2mat(permute(rawdata(:,1),[3 2 1])));
rawdata=reshape(rawdata,[sizeX,sizeY,sizeZ,sizeH]);
rawdata=permute(rawdata,[1,2,4,3]);
param.w=max(sizeX,sizeY);
if sizeX~=sizeY
    temp=zeros(param.w,param.w,sizeH,sizeZ);
    temp(1:sizeX,1:sizeY,:,:)=rawdata;
    param.rawdata=temp;
    clear temp
end
end
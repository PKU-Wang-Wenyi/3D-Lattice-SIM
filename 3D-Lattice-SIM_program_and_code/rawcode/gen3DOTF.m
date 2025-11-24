function [OTF,param]=gen3DOTF(obj,nz)
param.sizeX=obj.w;
param.sizeY=obj.w;
param.sizeZ=nz;
param.pixelsize=obj.pixelsize;% in nm
param.Zstep=obj.Zstep;% in nm
param.lambda=obj.lambda;% in nm
param.exlambda=obj.exlambda;% in nm
param.imgsize=[param.sizeX,param.sizeY,param.sizeZ];
param.medium_refractive_index=1.518;%nikon
param.Sample_refractive_index=1.518;%water
param.coverglass_refractive_index=1.518;%schott-D263M-datasheet-en.pdf
% 
% param.medium_refractive_index=1.4;%nikon
% param.Sample_refractive_index=1.5;%water
% param.coverglass_refractive_index=1.518;%schott-D263M-datasheet-en.pdf
% param.obj_workingdistence=140;%in umto nm
param.obj_workingdistence=0;%in umto nm
param.nz=nz;
param.w=obj.w;
param.NA=obj.NA;

[OTF,~] = get_modelOTF_fast(param);
OTF=gather(OTF);
end
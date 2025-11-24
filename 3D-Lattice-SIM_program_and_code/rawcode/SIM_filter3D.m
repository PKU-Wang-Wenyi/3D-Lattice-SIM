function SIM=SIM_filter3D(s_class,obj,FFTSIM,qvector)
nz_RAW=obj.sizeZ;
if mod(nz_RAW,2)==0
    nz=nz_RAW+3;
else
    nz=nz_RAW+4;
end
FFTSIM=ifftn(ifftshift(FFTSIM));
FFTSIM_pad=zeros(size(FFTSIM,1),size(FFTSIM,2),nz,'like',FFTSIM);
FFTSIM_pad(:,:,1:nz_RAW)=FFTSIM;
clear FFTSIM
FFTSIM=fftshift(fftn(FFTSIM_pad));
clear FFTSIM_pad
[OTF3D,otfobj]=gen3DOTF(obj,nz);
qvector1=qvector./sqrt(3);
kxy = obj.w.*obj.pixelsize./sqrt(sum(qvector1.^2));
lambda = obj.RefractiveIndex/obj.exlambda;
kz = (lambda-sqrt(lambda^2-1/kxy.^2)).*nz.*obj.Zstep;
[OTFcombine,notch,OTF_double]=shiftOTF(OTF3D,s_class,obj,kz,qvector);
obj.kz=kz;
[filter] = filter3D(OTFcombine,s_class,otfobj,obj.w1,obj.w2,obj.Zstep,notch,OTF_double,obj,qvector);
FFTSIM=FFTSIM.*filter;
SIM=real(ifftn(ifftshift(FFTSIM)));
SIM=SIM(:,:,1:1:nz_RAW);
SIM(SIM<0)=0;
SIM=gather(SIM);
end
function reconZeissSIM(rawdata,obj)
obj.cutoff=1000/(0.5*obj.lambda/obj.NA);
obj.cyclesPerMicron=1000/(obj.w*obj.pixelsize);
obj.OtfProvider1=SimOtfProvider(obj,obj.NA,obj.lambda,1);
OTFo=single(obj.OtfProvider1.otf);
Kotf = OTFedgeF(OTFo);
PSFe = fspecial('gaussian',30,3.0);
WF=imresize(squeeze(mean(rawdata,3)),2,"nearest");
imwritestack_16(WF,obj.WFsavename);
rawdata=reshape(rawdata,[obj.w,obj.w,obj.nPhase*obj.sizeZ]);
rawdata=edgefilter(single(rawdata));
if obj.RL>0
    rawdata=uint16(rawdata);
    psf=single(abs(otf2psf(OTFo)));
    RL=obj.RL;
    parfor nz=1:size(rawdata,3)
        rawdata(:,:,nz)=deconvlucy(rawdata(:,:,nz),psf,RL);
    end
end
if obj.Dark>0
    rawdata=Dark(rawdata,obj.Dark);
    % rawdata=dark_sectioning(rawdata,obj.lambda,obj.NA,obj.pixelsize,obj.Dark);
end
rawdata=reshape(rawdata,[obj.w,obj.w,obj.nPhase,obj.sizeZ]);
s_class = single_orientation_class(OTFo,Kotf,rawdata,PSFe);
s_class= KapproxEstimationF(s_class);
[s_class.kAmean,coe] = kmeanF(s_class);
paramfiles=['saved_data_',num2str(obj.exlambda),'.mat'];
s_class.phaseA = phasesF(s_class);
s_class = PCMseparateF_3D(s_class);
s_class.modFac = ModulationFactorF(s_class);
[s_class,coe,~]=checkparam(s_class,coe,Kotf,paramfiles);
s_class.OBJpara = OBJpowerPara(s_class,OTFo,coe);
obj.OtfProvider2=SimOtfProvider(obj,obj.NA,obj.lambda,2.5);
Otf_filter=gpuArray.ones([2*s_class.w,2*s_class.w,size(s_class.Separated_all,3),1],'single');
Otf_filter(:,:,1,:) = applyOtf(Otf_filter(:,:,1,:),obj.OtfProvider2,1,0,0,1,0);
for i= 1:(size(s_class.Separated,3)-1)/2
    Otf_filter(:,:,(i-1)*2+2,:) = applyOtf(Otf_filter(:,:,(i-1)*2+2,:),obj.OtfProvider2,2,s_class.kAmean(i,2),s_class.kAmean(i,1),1,0);
    Otf_filter(:,:,(i-1)*2+3,:) = applyOtf(Otf_filter(:,:,(i-1)*2+3,:),obj.OtfProvider2,2,-s_class.kAmean(i,2),-s_class.kAmean(i,1),1,0);
end
co = 1;
s_class = WoFilterCenterF(s_class,co);
s_class = PCMfilteringF(s_class,co,Otf_filter);
Filter3D=0;
[SIM,Filter3D]=SIM_filter3D_nt(s_class,obj,s_class.fftDirectlyCombined,Kotf*coe,Filter3D);
imwritestack_16(SIM,obj.SIMsavename);
end
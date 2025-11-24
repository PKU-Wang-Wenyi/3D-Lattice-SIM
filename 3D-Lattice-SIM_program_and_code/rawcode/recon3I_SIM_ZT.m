function [SIM_all,WF_all]=recon3I_SIM_ZT(obj)
for numT=1:obj.rawinfo.sizeT
    obj.rawinfo.nPhase=obj.nPhase;
    rawdata=imreadstack_ZT(obj.rawinfo,numT);
    %%
    rawdata_slice_cut=obj.slice_cut;
    %%
    % rawdata=edgefilter(rawdata);
    [sizeX,sizeY,imgnum]=size(rawdata);
    obj.w=max(sizeX,sizeY);
    sizeZ=floor(imgnum/obj.nPhase);
    if sizeZ>80
        rawdata=rawdata(:,:,1:sizeZ*obj.nPhase);
        disp("too many Z will run out of memory!")
        rawdata=reshape(rawdata,[sizeX,sizeY,obj.nPhase,sizeZ]);
        sizeZ=ceil(sizeZ/2);
        obj.Zstep=obj.Zstep*2;
        rawdata=rawdata(:,:,:,1:2:end);
        rawdata=reshape(rawdata,[sizeX,sizeY,obj.nPhase*sizeZ]);
    else
        rawdata=rawdata(:,:,1:sizeZ*obj.nPhase);
    end
    rawdata=reshape(rawdata,[sizeX,sizeY,obj.nPhase,sizeZ]);
    if numT==1
        obj.sizeZ=sizeZ;
        obj.cutoff=1000/(0.5*obj.lambda/obj.NA);
        obj.cyclesPerMicron=1000/(obj.w*obj.pixelsize);
        obj.OtfProvider1=SimOtfProvider(obj,obj.NA,obj.lambda,1);
        OTFo=single(obj.OtfProvider1.otf);
        Kotf = OTFedgeF(OTFo);
        PSFe = fspecial('gaussian',30,3.0);
        tempWF=imresize(squeeze(mean(rawdata,3)),2,"nearest");
        % maxint=squeeze(max(max(tempWF,[],1),[],2));
        maxint=zeros(size(tempWF,3),1);
        for nztempWF=1:size(tempWF,3)
            % maxint(nztempWF)=func4(tempWF(:,:,nztempWF));
            maxint(nztempWF)=func6(tempWF(:,:,nztempWF));
        end
        % figure;plot(maxint);
        [~,maxpos]=max(maxint);
        cut_ind=[0,0];
    end
    WF=imresize(squeeze(mean(rawdata,3)),2,"nearest");
    % SIM_rawsize=zeros(size(WF));
    if rawdata_slice_cut>0
    [rawdata,cut_ind]=(rawdata_slice_cut_fun(rawdata,maxpos,rawdata_slice_cut,cut_ind,maxint));
    end
    obj.sizeZ=size(rawdata,4);
    rawdata=reshape(rawdata,[sizeX,sizeY,obj.nPhase*obj.sizeZ]);
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
    rawdata=reshape(rawdata,[sizeX,sizeY,obj.nPhase,obj.sizeZ]);
    if sizeX~=sizeY
        temp=zeros(obj.w,obj.w,obj.nPhase,obj.sizeZ,"like",rawdata);
        temp(1:sizeX,1:sizeY,:,:)=rawdata;
        rawdata=temp;
        clear temp
    end
    if numT==1
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
    else
        s_class = single_orientation_class_nt(s_class,rawdata);
        s_class = PCMseparateF_3D_nt(s_class);
        co = 1;
        s_class = WoFilterCenterF(s_class,co);
        s_class = PCMfilteringF(s_class,co,Otf_filter);
    end
    %%
    % if numT==1
    %     Otf_filter=ones(size(s_class.shifted(:,:,:,1)),'like',s_class.shifted);
    %     Otf_filter(:,:,1,:) = applyOtf(Otf_filter(:,:,1,:),obj.OtfProvider2,1,0,0,1,0);
    %     for i= 1:(size(s_class.Separated,3)-1)/2
    %         Otf_filter(:,:,(i-1)*2+2,:) = applyOtf(Otf_filter(:,:,(i-1)*2+2,:),obj.OtfProvider2,2,s_class.kAmean(i,2),s_class.kAmean(i,1),1,0);
    %         Otf_filter(:,:,(i-1)*2+3,:) = applyOtf(Otf_filter(:,:,(i-1)*2+3,:),obj.OtfProvider2,2,-s_class.kAmean(i,2),-s_class.kAmean(i,1),1,0);
    %     end
    % end
    % s_class.shifted=s_class.shifted.*Otf_filter;
    % s_class.shifted(:,:,1,:) = applyOtf(s_class.shifted(:,:,1,:),obj.OtfProvider2,1,0,0,1,0);
    % % fftDirectlyCombined = s_class.shifted(:,:,1,:);
    % for i= 1:(size(s_class.Separated,3)-1)/2
    %     s_class.shifted(:,:,(i-1)*2+2,:) = applyOtf(s_class.shifted(:,:,(i-1)*2+2,:),obj.OtfProvider2,2,s_class.kAmean(i,2),s_class.kAmean(i,1),1,0);
    %     s_class.shifted(:,:,(i-1)*2+3,:) = applyOtf(s_class.shifted(:,:,(i-1)*2+3,:),obj.OtfProvider2,2,-s_class.kAmean(i,2),-s_class.kAmean(i,1),1,0);
    %     % fftDirectlyCombined = fftDirectlyCombined +s_class.shifted(:,:,(i-1)*2+2,:)+s_class.shifted(:,:,(i-1)*2+3,:);
    % end
    % fftDirectlyCombined=(squeeze(sum(s_class.shifted,3)));
    % s_class.shifted=[];
    %%
    if numT==1
        Filter3D=0;
    end
    [SIM_all,Filter3D]=SIM_filter3D_nt(s_class,obj,s_class.fftDirectlyCombined,Kotf*coe,Filter3D);
    WF_all=WF;
    SIM_all=im2uint16(mat2gray(SIM_all));
    WF_all=im2uint16(mat2gray(WF_all));
    filename=['_T',num2str(numT)];
    imwritestack_16(SIM_all,strrep(obj.resultname,'TEMP_STR',['SIM',filename ]))
    imwritestack_16(WF_all,strrep(obj.resultname,'TEMP_STR',['WF',filename ]))
    % if numT==1
    %     SIM_all=SIM_rawsize;
    %     WF_all=WF;
    % else
    %     SIM_all=cat(3,SIM_all,SIM_rawsize);
    %     WF_all=cat(3,WF_all,WF);
    % end
    disp(num2str(numT./obj.rawinfo.sizeT))
end
end
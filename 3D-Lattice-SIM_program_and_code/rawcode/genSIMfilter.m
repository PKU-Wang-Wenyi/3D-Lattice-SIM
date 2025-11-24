function [Filter]=genSIMfilter(obj,cutoff)
if obj.sizeZ==1
    h=2*obj.w;
    w=h;
    otfTS=zeros(h,w,'single');
    otfTS=writeApoVector(otfTS,obj.OtfProvider,cutoff);      % Ideal OTF
    MaskF=gpuArray.zeros(h,w,'single');
    MaskF(otfTS~=0)=1;
    obj.imgSize=obj.w;obj.nrDirs=3;
    if obj.nPhase>7
        obj.nrBands=3;
        wFilter1=WienerFilterW1_3D(obj);
        wFilter2=WienerFilterW2_3D(obj);
    else
        obj.nrBands=2;
        wFilter1=WienerFilterW1_2D(obj);
        wFilter2=WienerFilterW2_2D(obj);
    end
    Wk1=otfTS./(wFilter1.wDenom+obj.w1^2);
    ApoFWHM=0.5*(cutoff-1);
    ApoFWHM=min(0.5,round(ApoFWHM*100)/100);
    apo= apodize_gauss([h,w], struct('rad',ApoFWHM));
    Wk2=apo./(wFilter2.wDenom+obj.w2^2);
    Filter=Wk1.*Wk2.*MaskF;
else
    Filter=SIM_filter3D_fusion(param,w1,w2,stepsize,sizeZ);
end
end
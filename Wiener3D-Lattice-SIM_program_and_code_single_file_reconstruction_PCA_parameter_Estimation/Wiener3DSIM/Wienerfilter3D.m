function [filter]=Wienerfilter3D (shift_OTF,SIM_param,obj,Wiener,stepsize,notch,OTF_double,param,qvector,ap_range)
Mask_final = zeros(size(shift_OTF),"single");
Mask_final(shift_OTF>1e-4)=1;
Apo=apodization(size(shift_OTF,1)/2, size(shift_OTF,3), qvector, ap_range, ...
    obj.Zstep, obj.exlambda);
WienerDenom = abs(Apo)./(abs(shift_OTF) + (Wiener)^2);
filter= WienerDenom .*Mask_final;
end

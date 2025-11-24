function [shift_OTF,notch,OTF]=shiftOTF(OTF,param,obj,kz,pitch)
OTF=placeFreq3D(OTF);
[notch]=getnotch(OTF,obj);
shift_OTF=OTF;
OTF_notch=OTF.*notch;
% figure;mesh(OTF_notch(:,:,floor(end/2)+1))
% figure;mesh(squeeze(OTF_notch(floor(end/2)+1,:,:)))
for i=1:length(param.kAmean)
    k=sqrt(param.kAmean(i,2).^2+param.kAmean(i,1).^2);
    if k<pitch*0.75
        shift_OTF=shift_OTF+spe_move(OTF_notch,[round([param.kAmean(i,2),param.kAmean(i,1)]),kz]);
        shift_OTF=shift_OTF+spe_move(OTF_notch,[round([-param.kAmean(i,2),-param.kAmean(i,1)]),kz]);
        shift_OTF=shift_OTF+spe_move(OTF_notch,[round([param.kAmean(i,2),param.kAmean(i,1)]),-kz]);
        shift_OTF=shift_OTF+spe_move(OTF_notch,[round([-param.kAmean(i,2),-param.kAmean(i,1)]),-kz]);
    else
        shift_OTF=shift_OTF+spe_move(OTF_notch,round([param.kAmean(i,2),param.kAmean(i,1)]));
        shift_OTF=shift_OTF+spe_move(OTF_notch,round([-param.kAmean(i,2),-param.kAmean(i,1)]));
    end
end
shift_OTF=real(shift_OTF);
end
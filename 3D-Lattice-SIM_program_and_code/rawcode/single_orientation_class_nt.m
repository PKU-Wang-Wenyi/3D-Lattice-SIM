function obj = single_orientation_class_nt(obj,Snoisy)
Snoisy=gpuArray(Snoisy);
% obj.fSnoisy_all=Snoisy;
for i=1:size(Snoisy,3)
    obj.Separated_all(:,:,i,:)=fftshift(fftn(Snoisy(:,:,i,:)));
end
% Zmedian=floor(size(Snoisy,4)/2)+1;
% obj.fSnoisy=obj.fSnoisy_all(:,:,:,Zmedian);
% clear Snoisy
end
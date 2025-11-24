function obj = single_orientation_class(OTFo,Kotf,Snoisy,PSFe)
obj.OTFo = OTFo;
obj.PSFe =PSFe;
obj.Kotf = Kotf;
obj.w = gpuArray(single(size(Snoisy,1))); % = OTFo size
obj.wo = obj.w/2;
x = linspace(0,obj.w-1,obj.w);
y = linspace(0,obj.w-1,obj.w);
[X,Y] = meshgrid(x,y);
obj.Cv = (X-obj.wo) + 1i*(Y-obj.wo);
obj.Ro = abs(obj.Cv);
Rg = obj.Ro.*(1/4); % 1/4 determined the range of notch-filter
obj.G = 1 - exp(-0.05*Rg.^1.2);
obj.fSnoisy_all=gpuArray.zeros(size(Snoisy),'single');
for i=1:size(Snoisy,3)
    obj.fSnoisy_all(:,:,i,:)=fftshift(fftn(Snoisy(:,:,i,:)));
end
Zmedian=floor(size(Snoisy,4)/2)+1;
obj.fSnoisy=obj.fSnoisy_all(:,:,:,Zmedian);
end
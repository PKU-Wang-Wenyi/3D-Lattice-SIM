function phase =  phasesF(obj)
npeak=(size(obj.Separated,3)-1)/2;
phase = zeros(npeak,1);
opt = 0;
for i = 1:npeak
    fSa1Tnoisy0 = obj.Separated(:,:,1);
    fSa1Tnoisy = obj.Separated(:,:,(i-1)*2+2);
    t = 2*obj.w;
    f_temp0 = gpuArray.zeros(t,t,'single');
    f_temp0(obj.wo+1:obj.w+obj.wo,obj.wo+1:obj.w+obj.wo) = fSa1Tnoisy0;

    f_temp =gpuArray.zeros(t,t,'single');
    f_temp(obj.wo+1:obj.w+obj.wo,obj.wo+1:obj.w+obj.wo) = fSa1Tnoisy;

    clear fSa1Tnoisy fSa1Tnoisy0;
    fSa1Tnoisy0 = f_temp0;
    fSa1Tnoisy = f_temp;
    clear f_temp f_temp0
    %----------------------------------
    if i == 1
        t = size(fSa1Tnoisy,1);
        u = linspace(0,t-1,t);
        v = linspace(0,t-1,t);
        [U,V] = meshgrid(u,v);
    end
    %                 phase(i) = obj.Ifreq2optF(obj.kAmean,fSa1Tnoisy,ifft2(fSa1Tnoisy),U,V,opt);
    phase(i) = Ifreq2optF(obj,obj.kAmean(i,:),fSa1Tnoisy0,ifft2(fSa1Tnoisy),U,V,opt);
end
end
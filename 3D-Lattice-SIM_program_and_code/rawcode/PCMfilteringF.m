function  [obj] = PCMfilteringF(obj,co,Otf_filter)
npeak=(size(obj.Separated,3)-1)/2;
t = 2*obj.w;
to = t/2;
u = linspace(0,t-1,t);
v = linspace(0,t-1,t);
[U,V] = meshgrid(u,v);
Aobj = obj.OBJpara(1);
Bobj = obj.OBJpara(2);
if ~isfield(obj, 'fftDirectlyCombined')
    obj.fftDirectlyCombined = gpuArray.zeros(t,t,1,size(obj.Separated_all,4),'single');
else
    obj.fftDirectlyCombined(:) = 0;
end
if ~isfield(obj, 'Temp')
    obj.Temp = gpuArray.zeros(t,t,1,size(obj.Separated_all,4),'single');
end
for i = 1:npeak
    fDp = obj.Separated_all(:,:,(i-1)*2+2,:);
    fDm = obj.Separated_all(:,:,(i-1)*2+3,:);
    kv = obj.kAmean(i,2) + 1i*obj.kAmean(i,1); % vector along illumination direction
    Rp = abs(obj.Cv-kv);
    Rm = abs(obj.Cv+kv);
    OBJp = exp(Aobj.*Rp + Bobj);
    OBJm = exp(Aobj.*Rm + Bobj);
    k3 = round(obj.kAmean(i,:));
    OBJp(obj.wo+1+k3(1),obj.wo+1+k3(2)) = 0.25*OBJp(obj.wo+2+k3(1),obj.wo+1+k3(2))...
        + 0.25*OBJp(obj.wo+1+k3(1),obj.wo+2+k3(2))...
        + 0.25*OBJp(obj.wo+0+k3(1),obj.wo+1+k3(2))...
        + 0.25*OBJp(obj.wo+1+k3(1),obj.wo+0+k3(2));
    OBJm(obj.wo+1-k3(1),obj.wo+1-k3(2)) = 0.25*OBJm(obj.wo+2-k3(1),obj.wo+1-k3(2))...
        + 0.25*OBJm(obj.wo+1-k3(1),obj.wo+2-k3(2))...
        + 0.25*OBJm(obj.wo+0-k3(1),obj.wo+1-k3(2))...
        + 0.25*OBJm(obj.wo+1-k3(1),obj.wo+0-k3(2));
    SFo = obj.modFac(i);
    [fDpf,~] = WoFilterSideLobeF(obj,fDp,co,OBJm,SFo);
    [fDmf,~] = WoFilterSideLobeF(obj,fDm,co,OBJp,SFo);
    if i==1
    obj.Temp(obj.wo+1:obj.w+obj.wo,obj.wo+1:obj.w+obj.wo,:) = obj.Separated_all(:,:,1,:);
    fDof = obj.Temp;
    end
    obj.Temp(obj.wo+1:obj.w+obj.wo,obj.wo+1:obj.w+obj.wo,:) = fDpf;
    fDpf = obj.Temp;
    obj.Temp(obj.wo+1:obj.w+obj.wo,obj.wo+1:obj.w+obj.wo,:) = fDmf;
    fDmf = obj.Temp;
    shiftmatrix=2*pi*(obj.kAmean(i,2)/t.*(U-to) + obj.kAmean(i,1)/t.*(V-to));
    fDpf = fftshift(fft2(ifft2(ifftshift(fDpf)).*exp( +1i.*shiftmatrix )));
    fDmf = fftshift(fft2(ifft2(ifftshift(fDmf)).*exp( -1i.*shiftmatrix )));
    Cv = (U-to) + 1i*(V-to);
    Ro = abs(Cv);
    Rp = abs(Cv-kv);
    k2 = sqrt(obj.kAmean(i,:)*obj.kAmean(i,:)');
    % frequency range over which corrective phase is determined
    Zmask = (Ro < 0.8*k2).*(Rp < 0.8*k2);
    % corrective phase
    Angle0 = angle( sum(sum( fDof.*conj(fDpf).*Zmask )) );
    % phase correction
    fDpf = exp(+1i*Angle0).*fDpf;
    fDmf = exp(-1i*Angle0).*fDmf;
    if i==1
        obj.fftDirectlyCombined=obj.fftDirectlyCombined+fDof.*Otf_filter(:,:,1,:);
    end
    obj.fftDirectlyCombined=obj.fftDirectlyCombined+fDpf.*Otf_filter(:,:,(i-1)*2+2,:);
    obj.fftDirectlyCombined=obj.fftDirectlyCombined+fDmf.*Otf_filter(:,:,(i-1)*2+3,:);
end
end
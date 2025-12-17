function  [obj] = PCMfilteringF1(obj,co)
npeak=(size(obj.Separated,3)-1)/2;
t = 2*obj.w;
to = t/2;
u = linspace(0,t-1,t);
v = linspace(0,t-1,t);
[U,V] = meshgrid(u,v);
Aobj = obj.OBJpara(1);
Bobj = obj.OBJpara(2);
% obj.shifted = gpuArray.zeros(t,t,size(obj.Separated_all,3),size(obj.Separated_all,4),'single');
% obj.shifted = zeros(t,t,size(obj.Separated_all,3),size(obj.Separated_all,4),'single');
% Np=gpuArray.zeros(1,1,size(obj.Separated_all,3)-1,size(obj.Separated_all,4),'single');
nPages = size(obj.Separated_all,3);            % 总页数
obj.shifted = cell(1, nPages);                   % cell 中存放每一页（gpuArray）
Temp = gpuArray.zeros(t,t,1,size(obj.Separated_all,4),'single');
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
    Temp(obj.wo+1:obj.w+obj.wo,obj.wo+1:obj.w+obj.wo,:) = obj.Separated_all(:,:,1,:);
    fDof = Temp;
    Temp(obj.wo+1:obj.w+obj.wo,obj.wo+1:obj.w+obj.wo,:) = fDpf;
    fDpf = Temp;
    Temp(obj.wo+1:obj.w+obj.wo,obj.wo+1:obj.w+obj.wo,:) = fDmf;
    fDmf = Temp;
    fDp1 = fftshift(fft2(ifft2(ifftshift(fDpf)).*exp( +1i.*2*pi*(obj.kAmean(i,2)/t.*(U-to) + obj.kAmean(i,1)/t.*(V-to)) )));
    fDm1 = fftshift(fft2(ifft2(ifftshift(fDmf)).*exp( -1i.*2*pi*(obj.kAmean(i,2)/t.*(U-to) + obj.kAmean(i,1)/t.*(V-to)) )));
    % figure;mesh(log(abs(fDp1(:,:,2)+1)))
    Cv = (U-to) + 1i*(V-to);
    Ro = abs(Cv);
    Rp = abs(Cv-kv);
    k2 = sqrt(obj.kAmean(i,:)*obj.kAmean(i,:)');
    % frequency range over which corrective phase is determined
    Zmask = (Ro < 0.8*k2).*(Rp < 0.8*k2);
    % corrective phase
    Angle0 = angle( sum(sum( fDof.*conj(fDp1).*Zmask )) );
    % phase correction
    fDp2 = exp(+1i*Angle0).*fDp1;
    fDm2 = exp(-1i*Angle0).*fDm1;
    % if i==1
    %     obj.shifted(:,:,1,:)=fDof ;
    % end
    % obj.shifted(:,:,(i-1)*2+2,:)=fDp2;
    % obj.shifted(:,:,(i-1)*2+3,:)=fDm2;
    if i == 1
        obj.shifted{1} = fDof;    % 中心页
    end
    % 目标索引（与原代码一致）
    idx1 = (i-1)*2 + 2;
    idx2 = (i-1)*2 + 3;
    obj.shifted{idx1} = fDp2;
    obj.shifted{idx2} = fDm2;
    % Np(:,:,(i-1)*2+1,:)=npDp;
    % Np(:,:,(i-1)*2+2,:)=npDm;
end
obj.shifted =cat(3, obj.shifted{:}); 
% clear obj.shifted
% obj.Np=Np;
% showft(obj.shifted(:,:,1,:))
% showft(obj.shifted(:,:,2,:))
% showft(obj.shifted(:,:,3,:))
% showft(obj.shifted(:,:,4,:))
% showft(obj.shifted(:,:,5,:))
% showft(obj.shifted(:,:,6,:))
% showft(obj.shifted(:,:,7,:))
% showft(obj.shifted(:,:,8,:))
% showft(obj.shifted(:,:,9,:))
% showft(obj.shifted(:,:,10,:))
% showft(obj.shifted(:,:,11,:))
% showft(obj.shifted(:,:,12,:))
% showft(obj.shifted(:,:,13,:))
end
function [obj] = KapproxEstimationF(obj)
if size(obj.fSnoisy,3)<13
    [Separated]=PCMseparate_tempF_2D(obj);
    Separated_temp=Separated.*obj.G.*(obj.Ro<obj.Kotf);
    Zmask = (obj.Ro>0.38*obj.Kotf ).*(obj.Ro<0.95*obj.Kotf );  % Spread spectrum coef <1.9
    kAo = zeros(3,2);
    kAo(1,:) = peakCCidx(obj,Separated_temp(:,:,1),Separated_temp(:,:,2),Zmask);
    kAo(2,:) = peakCCidx(obj,Separated_temp(:,:,1),Separated_temp(:,:,4),Zmask);
    kAo(3,:) = peakCCidx(obj,Separated_temp(:,:,1),Separated_temp(:,:,6),Zmask);
    obj.Separated=Separated;
    obj.kAo = kAo;
else
    % [Separated]=PCMseparate_tempF_3D(obj);
    [Separated]=PCMseparate_tempF_3D1(obj);
    Separated_temp=Separated.*obj.G.*(obj.Ro<obj.Kotf);
    % figure;mesh(log(abs(Separated_temp(:,:,2))+1))
    % figure;mesh(log(abs(Separated_temp(:,:,4))+1))
    % figure;mesh(log(abs(Separated_temp(:,:,6))+1))
    % figure;mesh(log(abs(Separated_temp(:,:,8))+1))
    % figure;mesh(log(abs(Separated_temp(:,:,10))+1))
    % figure;mesh(log(abs(Separated_temp(:,:,12))+1))
    
    Zmask1 = (obj.Ro>0.1*obj.Kotf ).*(obj.Ro<=0.6*obj.Kotf );
    Zmask2 = (obj.Ro>0.5*obj.Kotf ).*(obj.Ro<1.2*obj.Kotf );
    % showft(Separated_temp(:,:,8))
    % showft(Separated_temp(:,:,8).*Zmask1)
    % Zmask_all = (obj.Ro>0.05*obj.Kotf ).*(obj.Ro<1.2*obj.Kotf );
    kAo = zeros(6,2);
    kAo(1,:) = peakCCidx(obj,Separated_temp(:,:,1),Separated_temp(:,:,2),Zmask1);
    kAo(2,:) = peakCCidx(obj,Separated_temp(:,:,1),Separated_temp(:,:,4),Zmask2);
    kAo(3,:) = peakCCidx(obj,Separated_temp(:,:,1),Separated_temp(:,:,6),Zmask1);
    kAo(4,:) = peakCCidx(obj,Separated_temp(:,:,1),Separated_temp(:,:,8),Zmask1);
    kAo(5,:) = peakCCidx(obj,Separated_temp(:,:,1),Separated_temp(:,:,10),Zmask2);
    kAo(6,:) = peakCCidx(obj,Separated_temp(:,:,1),Separated_temp(:,:,12),Zmask2);
    obj.Separated=Separated;
    obj.kAo = kAo;
end
end
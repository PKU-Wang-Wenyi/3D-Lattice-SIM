function [obj] = WoFilterCenterF(obj,co)
if ~isfield(obj, 'Wiener_Filter')
    OTFpower = obj.OTFo.*conj(obj.OTFo);
    NoiseFreq = obj.Kotf + 20;
    Zo = obj.Ro>NoiseFreq;
    nNoise = obj.Separated(:,:,1).*Zo;
    NoisePower = sum(sum( nNoise.*conj(nNoise) ))./sum(sum(Zo));
    A = obj.OBJpara(1);
    B = obj.OBJpara(2);
    LogFe = A*obj.Ro + B ;
    OBJo = exp(LogFe);
    OBJpower = OBJo.^2 - 1.0*NoisePower;
    SFo = 1;
    obj.Wiener_Filter=(SFo.*conj(obj.OTFo)./NoisePower)./((SFo.^2).*OTFpower./NoisePower + co./OBJpower);
    % obj.npDo=NoisePower ;
end
obj.Separated_all(:,:,1,:) = obj.Separated_all(:,:,1,:).*obj.Wiener_Filter;
end
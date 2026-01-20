function [ separate ] = separateBandshifi_lattice( IrawFFT,phaOff,bands,fac )
    phaPerBand=bands;
    phases=zeros(1,phaPerBand);
    for p=1:phaPerBand
        phases(p)=(2*pi*(p-1))/phaPerBand+phaOff;                                                                              
    end
    separate=separateBands_final(IrawFFT,phases,bands,fac);
end

function [separate] = separateBands_final(IrawFFT,phases,bands,fac)
    if fac==0
        fac=zeros(1,bands);
        fac(1)=1;
        for I=2:bands
            fac(I)=0.5;
        end
    else
        for I=2:bands
            fac(I)=fac(I)*0.5;
        end
    end

    comp=zeros(1,bands);
    comp(1)=0;
    for I=2:(bands+1)/2
        comp((I-1)*2)=I-1;
        comp((I-1)*2+1)=-(I-1);
    end
    compfac=zeros(1,bands);
    compfac(1)=fac(1);
    for I=2:(bands+1)/2
        compfac((I-1)*2)=fac(I);
        compfac((I-1)*2+1)=fac(I);
    end

    W=exp(1i*phases'*comp);
    for I=1:bands
        W(I,:)=W(I,:).*compfac;
    end

    length=size(phases,2);
    siz=size(IrawFFT(:,:,1));
    S=reshape(IrawFFT,[prod(siz),length])*pinv(W)';
    Sk=reshape(S,[siz,bands]);
    separate=Sk;
end

function image_out = spe_move(image_in,shiftnum)
if length(shiftnum)==2
    [nx,ny,nz]=size(image_in);
    image_out =gpuArray.zeros([nx,ny,nz],'single');
    for i=1:nz
        temp=ifftshift(ifft2(image_in(:,:,i)));
        % figure;mesh(log(abs(temp)+1))
        t =gpuArray(single(nx));
        to = t/2;
        u = linspace(0,t-1,t);
        v = linspace(0,t-1,t);
        [U,V] = meshgrid(u,v);
        tempout=(temp).*exp( +1i.*2*pi*(shiftnum(1)/t.*(U-to) + shiftnum(2)/t.*(V-to)) );
        tempout=(fft2(fftshift(tempout)));
        image_out(:,:,i)= tempout;
    end
elseif length(shiftnum)==3
    [M, N, P] = size(image_in);
    % image_out =zeros([M, N, P]);
    M=gpuArray(single(M));
    N=gpuArray(single(N));
    P=gpuArray(single(P));
    fftVolume = ifftshift(ifftn(image_in));
    u= linspace(0,M-1,M);
    v= linspace(0,N-1,N);
    w= linspace(0,P-1,P);
    [X, Y, Z] = meshgrid(u, v, w);
    shiftX=shiftnum(1)/M.*(X-M/2);
    shiftY=shiftnum(2)/N.*(Y-N/2);
    shiftZ=shiftnum(3)/P.*(Z-P/2);
    shiftOperatorX = exp(+1i * 2 * pi * shiftX );
    shiftOperatorY = exp(+1i * 2 * pi * shiftY );
    shiftOperatorZ = exp(+1i * 2 * pi * shiftZ );
    shiftedFFT = fftVolume .* shiftOperatorX .* shiftOperatorY .* shiftOperatorZ;
    image_out = fftn(fftshift(shiftedFFT));
end
% image_out =double(shift(image_in,shiftnum));
% figure;mesh(log(abs(image_out)+1))
% figure;mesh(log(abs(image_out1)+1))
end
function [notch]=getnotch(OTF,obj)
[nx,ny,nz]=size(OTF);
nx=gpuArray(single(nx));
ny=gpuArray(single(ny));
nz=gpuArray(single(nz));
x=linspace(0,nx-1,nx);
y=linspace(0,ny-1,ny);
z=linspace(0,nz-1,nz);
[X,Y,Z]=meshgrid(x,y,z);
% Ro = sqrt((X-nx/2).^2*0.1+(Y-ny/2).^2*0.1+(Z-floor(nz/2)).^4);
Ro = sqrt((X-nx/2).^2*0.1+(Y-ny/2).^2*0.1+(Z-floor(nz/2)).^2*0.1);
% Rg = Ro.*(8);
Rg = 20*Ro.*(1-1*obj.attStrength);
notch = 1 - obj.attStrength*exp(-0.05*Rg);
% figure;mesh(notch(:,:,floor(nz/2)+1))
% figure;mesh(squeeze(notch(floor(nx/2)+1,:,:)))
end
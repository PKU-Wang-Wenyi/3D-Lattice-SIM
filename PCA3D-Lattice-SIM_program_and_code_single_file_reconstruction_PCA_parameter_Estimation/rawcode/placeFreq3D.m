function [out] = placeFreq3D( in )
siz=size(in);
w=siz(2);
h=siz(1);
out=gpuArray.zeros(2*w,2*h,size(in,3),size(in,4),'single');
out(h/2+1:h+h/2,w/2+1:w+w/2,:,:)=in;
end
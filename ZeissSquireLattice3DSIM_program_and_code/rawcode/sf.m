function sf(in)
in=squeeze(in);
if size(in,3)>1
    figure;mesh(log(abs(in(:,:,floor(end/2)+1))+1));
else
    figure;mesh(log(abs(in)+1));
end
end
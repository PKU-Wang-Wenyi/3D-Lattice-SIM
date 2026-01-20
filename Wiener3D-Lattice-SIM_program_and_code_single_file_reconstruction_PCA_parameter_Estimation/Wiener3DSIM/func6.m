function fz=func6(a)
%特南鲍姆梯度
% [Width,Height]=size(a);
Tx=[1 2 1;0 0 0;-1 -2 -1];
Ty=[-1 0 1;-2 0 2;-1 0 1];

Ax=conv2(a,Tx);
Ay=conv2(a,Ty);
fz = mean((Ax(:).^2) .* (Ay(:).^2));
% temp=0;
% for j=1:Width
%     for k=1:Height
%         A=(Ax(j,k)^2)*(Ay(j,k)^2);
%         temp=temp+A;
%     end
% end
% fz=temp/(Width*Height);

end
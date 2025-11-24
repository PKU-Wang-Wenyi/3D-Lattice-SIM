function fz=func4(a)
%拉普拉斯能量


E_lplc=[-1 -4 -1;-4 20 -4;-1 -4 -1];
A=conv2(a,E_lplc);
A=A.^2;
fz=sum(A(:));

end
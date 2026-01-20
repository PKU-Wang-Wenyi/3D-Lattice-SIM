function k2fa = peakCCidx(obj,fAo0,fAp0,Zmask)
            
            Zo = double( obj.Ro<obj.Kotf );
            
            %tic
            %Tbegin = toc;
            
            Zn = fftshift( ifft2( fft2(Zo).*conj(fft2(Zo)) ) );
            CC = fftshift( ifft2( fft2(fAo0.*Zo).*conj(fft2(fAp0.*Zo)) ) );
%             rho = abs(CC)./abs(Zn);
            rho = abs(CC)./(abs(Zn)+eps); 

                       
            temp = rho.*Zmask.*Zo;
            % figure;mesh(temp);
            [~, v] = max(temp(:));
            
            pmax = mod(v-1,obj.w) + 1; % row
            qmax = (v-pmax)/obj.w + 1; % column
            px0 = qmax - (obj.wo+1);
            py0 = pmax - (obj.wo+1);
            k2fa = [py0 px0]; % nearest pixel approximation

end
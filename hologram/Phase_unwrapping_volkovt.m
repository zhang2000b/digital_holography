
%Purpose: Performs 2D phase unwrapping; J. Magn. Reson. Imaging 27 (2008) 649
%Inputs: 'Pha_wr', wrapped phase
%Outouts: 'Pha_unwr', unwrapped phase;

function Pha_unwr=Phase_unwrapping_volkovt(Pha_wr)

        linhas=size(Pha_wr,1); 
        colunas=size(Pha_wr,2); 

        image=Pha_wr; 
        image=padarray(image,[abs(linhas-colunas),abs(linhas-colunas)],'replicate','post');%Faz padarray se a imagem nao for quadrada
        if mod(size(image,1),2)~=0
             image=padarray(image,[1,0],'replicate','post');%Faz padarray se as dimensoes nao forem par
             %在第一个维度上添加一行
        end
        if mod(size(image,2),2)~=0
             image=padarray(image,[0,1],'replicate','post');%Faz padarray se as dimensoes nao forem par
             %在第二个维度上添加一行
        end

        [sy,sx]=size(image);
        p=linspace(-sx/2,sx/2-1,sx); %x coordinate; 
        q=linspace(-sy/2,sy/2-1,sy); %y coordinate; 
        [p,q]=meshgrid(p,q); %2D coordinate; 
        N=sx;

        Gpsix=cos(image).*ifft2(fftshift((p*2*pi).*fftshift(fft2(sin(image)))))-sin(image).*ifft2(fftshift((p*2*pi).*fftshift(fft2(cos(image)))));
        Gpsiy=cos(image).*ifft2(fftshift((q*2*pi).*fftshift(fft2(sin(image)))))-sin(image).*ifft2(fftshift((q*2*pi).*fftshift(fft2(cos(image)))));

        Gfix=ifft2(fftshift((p*2*pi).*fftshift(fft2(image))));
        Gfiy=ifft2(fftshift((q*2*pi).*fftshift(fft2(image))));

        Gkx=(Gpsix-Gfix)/2/pi;

        Gky=(Gpsiy-Gfiy)/2/pi;

        p(sy/2+1,sx/2+1)=1;
        b=(fftshift(fft2(Gkx)).*p+fftshift(fft2(Gky)).*q)./(p.^2+q.^2);

        k=(ifft2(fftshift(b))); %Inverse Fourier transform.

        Pha_unwr=double(real(k)+image);
        Pha_unwr=Pha_unwr(1:size(Pha_wr,1),1:size(Pha_wr,2));

end


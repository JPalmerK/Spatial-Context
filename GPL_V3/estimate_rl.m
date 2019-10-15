function [spec_rl]=estimate_rl(cm,sp_orig,parm);
 

[k1]=find(cm);%indices of non-zero elements
new_slate=zeros(size(cm));
new_slate(k1)=sp_orig(k1);
pad=conv2(new_slate,ones(3,3)/9,'same');
k1x=find(pad);
noise_slate=sp_orig;
noise_slate(k1x)=0;%take out all values from call, left with only noise

[x1,y1]=find(noise_slate);%find non-zero elements, (row, col)
k0x=find(noise_slate);%indices of non-zero elements
z1=noise_slate(k0x);%only include non-zero vals
[x2,y2]=find(cm);
z2=griddata(x1,y1,z1,x2,y2,'nearest');%interpolates noise onto the template of the call
new_slate(k1)=new_slate(k1)-z2;%subtract noise from call
k=find(new_slate<0);
new_slate(k)=0;

%noise_slate(k1)=z2;  %imagesc(log([noise_slate,new_slate]))
spec_rl=2*new_slate.^2;%multiplication by 2 corrects for negative frequencies
spec_rl=spec_rl/parm.sample_freq/sum(hamming(parm.fftl).^2);%normalize by freq bin width
spec_rl=sum(spec_rl,1);%integrate across frequency
spec_rl=spec_rl*parm.sample_freq/parm.fftl;%scale by df (Hz)
spec_rl= sqrt(mean(spec_rl)); %RMS (uPa), to get dB, take 20*log10(RMS)

% spec_rl_sel=sum(spec_rl); %SEL, to get dB take 10*log10(SEL)

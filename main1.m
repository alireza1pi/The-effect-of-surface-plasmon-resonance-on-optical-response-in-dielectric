%figure1 code
close all
clc
clear

c=3e8;
syms L
rcore=10e-9;
rshell=15e-9;
nout=1.3;
nshell=1.8;
ncore=1.5;
ecore=ncore.^2;
eshell=nout.^2;
em=9.8;
kout=9;
f=(rcore/rshell)^(1/3);
alpha=4*pi*rshell.^2+((eshell-em)*(ecore+2*eshell)+f*(ecore-eshell)*(em+2*eshell))/(eshell+2*em)*(ecore+2*eshell)+f*(2*eshell-2*em)*(ecore-eshell);

%alpha_w=@(w) rcore.^3.*((eshell.*w-em)./(eshell.*w+2.*em));
sigma_ext=@(w) ((2*pi*w)/(c))*imag(alpha);


%{
s=[];
ww=[];
www=1.5;
for i=1:15
    
    
    a=sigma_ext(www);
    s(end+1)=a;
    ww(end+1)=www;
    www=www+.1;
end
plot(ww,s)
%}    
landa=15e-9;




s=[];
ww=[];
www=1.5;
for i=1:30
    w=www;
    
    s1=0;
    for L=1:5
        ss=(2*L+1)*real(aL(ncore,nout,nshell,kout,rcore,rshell,L)+bL(ncore,nout,nshell,kout,rcore,rshell,L));
        s1=s1+ss;
    end
    
    sigma_ext=(landa.^2/(2*pi))*s1;
    a=sigma_ext;
    s(end+1)=a;
    ww(end+1)=www;
    www=www+.05;
end






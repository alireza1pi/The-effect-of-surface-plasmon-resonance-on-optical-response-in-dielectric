function [bL,aL,AL,BL]=miefunc(ncore,nout,nshell,kout,rcore,rshell,L)

syms x;

phi_L=@(z) z*besselj(z,L);
X_L=@(z) z*bessely(z,L);
eps_L=@(z) z*besselh(z,L);

phi_L_p=@(z) (phi_L(z-1)-phi_L(z+1))/2;
X_L_p=@(z) (X_L(z-1)-X_L(z+1))/2;
eps_L_p=@(z) (eps_L(z-1)-eps_L(z+1))/2;



m1=ncore/nout;
m2=nshell/nout;
x=kout*rcore;
y=kout*rshell;


AL1=m2*phi_L(m2*x)*phi_L_p(m1*x)-m1*phi_L_p(m2*x)*phi_L(m1*x);
AL2=m2*X_L(m2*x)*phi_L_p(m1*x)-m1*X_L_p(m2*x)*phi_L(m1*x);
AL=AL1/AL2;

BL1=m2*phi_L(m1*x)*phi_L_p(m2*x)-m1*phi_L_p(m2*x)*phi_L_p(m1*x);
BL2=m2*X_L_p(m2*x)*phi_L(m1*x)-m1*X_L(m2*x)*phi_L_p(m1*x);
BL=BL1/BL2;



aL1=(phi_L(y)*(phi_L_p(m2*y)-AL*X_L_p(m2*y))-m2*phi_L_p(y)*(phi_L(m2*y)-AL*X_L(m2*y)));
aL2=(eps_L(y)*(phi_L_p(m2*y)-AL*X_L_p(m2*y))-m2*eps_L_p(y)*(phi_L(m2*y)-AL*X_L(m2*y)));
aL=aL1/aL2;


bL1=m2*phi_L(y)*(phi_L_p(m2*y)-BL*X_L_p(m2*y))-phi_L_p(y)*(phi_L(m2*y)-BL*X_L(m2*y));
bL2=m2*eps_L(y)*(phi_L_p(m2*y)-BL*X_L_p(m2*y))-eps_L_p(y)*(phi_L(m2*y)-BL*X_L(m2*y));
bL=bL1/bL2;




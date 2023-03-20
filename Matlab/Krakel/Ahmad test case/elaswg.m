z1=200;
z2=600;
rho1=1;
rho2=1.5;
f=25;
omega=2*pi*f;
cw=1500;
cp0=1700;
cs0=700;
rhow=1;
[k u v up vp sigma_zz sigma_zx phi psi z]=evs_efuncs(z1,z2,rho1,rho2,f,cw,cp0,cs0);
zs=180;
zr=30;
r=1:10:10000;
lmax=46;
l=1:lmax;
sigma_zz=sigma_zz(:,l);
k=k(l);
% Interpolate the modes at the source and receiver locations
Ps=-interp1(z,sigma_zz,zs);
Pr=-interp1(z,sigma_zz,zr);
kr=r.'*k;
Field=sqrt(2*pi)/rhow*(exp(i*kr)./sqrt(kr))*(Ps.*Pr).';
% We are plotting the potential, not the pressure
tloss=-20.*log10(abs(Field/(rhow*omega^2)));
plot(r/1000,tloss,'r');

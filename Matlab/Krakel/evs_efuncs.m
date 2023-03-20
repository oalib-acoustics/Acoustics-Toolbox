function [ev u v up vp tau_zz tau_zx phi psi z]=evs_efuncs(z1,z2,rho1,rho2,f,cw,cp0,cs0)
omega=2*pi*f;
kw=omega/cw;
atten_p=0; 
atten_s=0;
kp=omega/cp0*(1.+atten_p/54.58*i);
ks=omega/cs0*(1.+atten_s/54.58*i);
cp=omega/kp;
cs=omega/ks;
lambda1 = rho1*(cw^2);
lambda2 = rho2*(cp^2-2.*cs^2);
mu2 = rho2*cs^2;
dx=0.000001;
x=dx;
m=60000;
fact=0.5;
del=(fact-dx)/m;
f1=evfunc(x,kw,kp,ks,z1,z2,lambda1,lambda2,mu2);
x1=dx;
l=0;
kr=[];
for j=1:m
    dx=dx+del;
    x=dx;
    x2=x;
    f2=evfunc(x,kw,kp,ks,z1,z2,lambda1,lambda2,mu2);
    if f1*f2 <0
        alpha(1)=x1;
        alpha(2)=x2;
        t=fzero(@(x) evfunc(x,kw,kp,ks,z1,z2,lambda1,lambda2,mu2),alpha);
        l=l+1;
        kr(l)=t;
    end
    f1=f2;
    x1=x2;
 end
 lmax=l;
 ev=sort(kr(1:lmax),'descend');
 kr=ev;
 % Compute the eigenfunctions
 k1=sqrt(kw^2-kr.^2);
 k2=sqrt(kp^2-kr.^2);
 g2=sqrt(ks^2-kr.^2);
%	xr1 is the horizontal (x) displacement
%	xr2 is the vertical (z) displacement
%	xr3 is the shear stress
%	xr4 is the compressional sterss
%	xr5 is the compressional potential
%	xr6 is the shear potential
    nmax=1000;
    dz=z2/nmax;
    j=0:nmax;
    z=dz*j';
    zw=z(find(z<=z1));
    zb=z(find(z>z1));
    Mu=[zeros(1,length(zw)) mu2*ones(1,length(zb))].';
    Lambda=[lambda1*ones(1,length(zw)) lambda2*ones(1,length(zb))].';
    xsq=kr.^2;
    x1 = lambda1*kw^2;
	x2 = (lambda2+2.*mu2)*kp^2-2.*mu2*xsq;
	x3 = 2.*mu2*xsq;
	x4 = 2.*xsq-ks^2;
	y1 = -2./ks^2;
	y2 = x4/2.*y1;
	z21 = z2-z1;
	ck1 = cos(k1*z1);
	ck2 = cos(k2*z21);
	cg2 = cos(g2*z21);
	sk1 = sin(k1*z1);
	sk2 = sin(k2*z21);
	sg2 = sin(g2*z21);
	xk2 = ck2./(sk2.*k2);
	xg2 = sg2./(g2.*cg2);
	yy = k2.*sk2.*cg2.*(x2.*x4.*xk2.*xg2-2.*x3);
	d = ck1./(x2.*x4.*ck2.*sg2./g2-2.*k2.*x3.*sk2.*cg2);
	alpha1 = ones(size(zb))*(y2.*d.*x2.*x4.*sg2./g2);
	beta1 = ones(size(zb))*(-2.*x3.*y2.*d.*cg2);
	zeta1 = ones(size(zb))*(2.*x3.*y2.*d);
	alpha2 = ones(size(zb))*(-2.*k2.*x3.*y1.*d.*sk2);
	beta2 = ones(size(zb))*(-y1.*d.*x2.*x4.*ck2);
	zeta2 = ones(size(zb))*(y1.*d.*x2.*x4);
    xr5w = sin(zw*k1)./(ones(size(zw))*k1);
    xr5b = alpha1.*sin((zb-z2)*k2)./(ones(size(zb))*k2)+beta1.*cos((zb-z2)*k2)+zeta1.*cos((zb-z1)*k2);
    xr6w = zeros(size(xr5w));
    xr6b = alpha2.*cos((zb-z2)*g2)+beta2.*sin((zb-z2)*g2)./(ones(size(zb))*g2)+zeta2.*sin((zb-z1)*g2)./(ones(size(zb))*g2);
    xr5pw = cos(zw*k1);
    xr5pb= alpha1.*cos((zb-z2)*k2)-beta1.*(ones(size(zb))*k2).*sin((zb-z2)*k2)-zeta1.*(ones(size(zb))*k2).*sin((zb-z1)*k2);
    xr6pw = zeros(size(xr6w));
    xr6pb = -alpha2.*(ones(size(zb))*g2).*sin((zb-z2)*g2)+beta2.*cos((zb-z2)*g2)+zeta2.*cos((zb-z1)*g2);
    xr5ppw = (-ones(size(zw))*k1.^2).*xr5w;
    xr5ppb = (-ones(size(zb))*k2.^2).*xr5b;
    xr6ppw = zeros(size(xr6w));
    xr6ppb = (-ones(size(zb))*g2.^2).*xr6b;
    xr4w = -(ones(length(zw),length(k1))*x1).*xr5w;
    xr5 = [xr5w; xr5b];
    xr6 = [xr6w; xr6b];
    xr5p = [xr5pw; xr5pb];
    xr6p = [xr6pw; xr6pb];
    xr5pp = [xr5ppw; xr5ppb];
    xr6pp = [xr6ppw; xr6ppb];
    xr1 = xr5-xr6p;
    xr2 = xr5p - (ones(length(z),1)*xsq).*xr6;
    xr3 = (Mu*ones(1,length(kr))).*(2*xr5p - (ones(size(z))*x4).*xr6);
    xr4b = -(ones(size(zb))*x2).*xr5b - (ones(size(zb))*x3).*xr6pb;
    xr4 = [xr4w; xr4b];
    v = (ones(size(z))*kr).*xr1;
    u = xr2;
    tau_zz = xr4;
    tau_zx = (ones(size(z))*kr).*xr3;
    vp = (ones(size(z))*kr).*(xr5p - xr6pp);
    up = xr5pp - (ones(size(z))*xsq).*xr6p;
    phi = xr5;
    psi = xr6;
    
% Normalize the eigenfunctions

ker1 = (Mu*ones(1,length(kr))).*u.^2+((Lambda+2*Mu)*ones(1,length(kr))).*v.^2;
ker2 = (Mu*ones(1,length(kr))).*(2*u.*vp);
ker3 = (Lambda*ones(1,length(kr))).*(2*v.*up);
Kernel = 2*(ones(size(z))*kr).*ker1 + ker2 - ker3;

% Use the trapezoidal rule to integrate this

qvec = [dz/2 dz*ones(1,length(z)-2) dz/2];

for j=1:length(kr)
    S(j)=qvec*Kernel(:,j);
end

cof=sqrt(2.*kr./abs(S));

% The normalized eigenfunctions

 v = v*diag(cof);
 u = u*diag(cof);
 tau_zz = tau_zz*diag(cof);
 tau_zx = tau_zx*diag(cof);
 phi = phi*diag(cof);
 psi = psi*diag(cof);
 vp = vp*diag(cof);
 up = up *diag(cof);

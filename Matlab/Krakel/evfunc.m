    function y=evfunc(x,kw,kp,ks,z1,z2,lambda1,lambda2,mu2)
	xsq =  x^2;
    k1=sqrt(kw^2-xsq);
    k2=sqrt(kp^2-xsq);
    g2=sqrt(ks^2-xsq);
 	
	x1 = lambda1*kw^2;
	x2 = (lambda2+2.*mu2)*kp^2-2.*mu2*xsq;
	x3 = 2.*mu2*xsq;
	x4 = 2.*xsq-ks^2;
	y1 = -2./ks^2;
	z21 = z2-z1;
		

	ck1 = cos(k1*z1);
	ck2 = cos(k2*z21);
	cg2 = cos(g2*z21);
	sk1 = sin(k1*z1);
	sk2 = sin(k2*z21);
	sg2 = sin(g2*z21);
	
     
	y = 2.*x1*sk1/k1*(x2*x4*ck2*sg2/g2-2.*x3*k2*sk2*cg2)-...
        ck1*(4.*y1*x2*x3*x4*(1.-ck2*cg2)-y1*sk2*sg2/(k2*g2)*...
        (4.*x3^2*g2^2*k2^2+x4^2*x2^2));

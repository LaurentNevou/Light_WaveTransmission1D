%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% FabryPerot Cavity structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nL=1;
nR=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fabri-Perot cavity 

n1=1;
n2=3.6+0.0i;
lambda0=940e-9;      % Central wavelength
l1=1e-6;
l2=10*lambda0/(2*abs(n2));

alpha2=4*pi*imag(n2)./lambda;

layer=[
l1   n1
l2   n2
l1   n1
];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Formula computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r12 = (n2-n1)/(n1+n2);
r21 = (n1-n2)/(n1+n2);
RR = r12^2;

t12 = 2*n1/(n1+n2);
t21 = 2*n2/(n1+n2);
TT = t12*t21;

theta=pi;
phi=2*pi*n2*l2./lambda;

delta=2*(phi-theta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see the Chapther of Vincenzo Savona for the formula in the book:
% Confined Photon Systems
% Fundamentals and Applications Lectures from the Summerschool Held in Cargèse, Corsica, 3–15 August 1998
% Linear Optical Properties of Semiconductor Microcavities with Embedded Quantum Wells, Vincenzo Savona, pages 173-242
% 3) The Fabry-Perot resonator p184
% https://link.springer.com/chapter/10.1007/BFb0104383
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Emmanuel Rosencher, Optoelectronic
% Complement to Chapter 9
% 9D) Fabry-Perot cavities and Bragg reflectors, page 437
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tfp = t12*t21 ./ ( 1 +  r12*r21*exp(2i*(phi-theta)) ) ;
Tfp = (abs(tfp)).^2 .* exp(-alpha2*l2);

rfp = -t21/t21' * ( (r21'+r12*exp(2i*(phi-theta))) ./ (1+r12*r21*exp(2i*(phi-theta))));
Rfp = (abs(rfp)).^2;

Tf=Tfp;
Rf=Rfp;

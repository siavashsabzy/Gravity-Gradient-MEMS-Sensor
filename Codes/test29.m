function ds = test29(t,s)
ds = zeros(21,1);
Tc1=0;Td1=0;Tc2=0;Td2=0;Tc3=0;Td3=0;
I=[0.017,0,0;0,0.055,0;0,0,0.055];
Is = 6.9*10^-8;
deltaIs = 6.9*10^-8;
ks = 1.62*10^-4;
d = 0; 
mu=398576.0576;
RE=6400;
J2=1.08263*10^-3;
q=[s(4) s(5) s(6) s(7)];
A= quat2dcm(q); %Inertial to body Transformation
r=[s(8);s(9);s(10)];
rb=A*r;
rbx=hankel(rb,rb(end:-1:1));
rbx(eye(numel(rb))==1) = rb(1);
Tg=(3*mu/(sqrt(s(8)^2+s(9)^2+s(10)^2)^5))*rbx*I*rb;              %GravityGradientTorque due to Moment of inertia
ds(1) = Tc1 + Td1 + Tg(1,1) + ( - s(2) * s(3) * ( I(3,3) - I(2,2) ))/I(1,1);
ds(2) = Tc2 + Td2 + Tg(2,1) + ( - s(1) * s(3) * ( I(1,1) - I(3,3) ))/I(2,2);
ds(3) = Tc3 + Td3 + Tg(3,1) + ( - s(1) * s(2) * ( I(2,2) - I(1,1) ))/I(3,3);
ds(4) = 0.5*( s(3)*s(5) - s(2)*s(6) + s(1)*s(7));
ds(5) = 0.5*(-s(3)*s(4) + s(1)*s(6) + s(2)*s(7));
ds(6) = 0.5*( s(2)*s(4) - s(1)*s(5) + s(3)*s(7));
ds(7) = 0.5*(-s(1)*s(4) - s(2)*s(5) - s(3)*s(6));
ds(8) = s(11);
ds(9) = s(12);
ds(10)= s(13);
ds(11)= (-mu/(sqrt(s(8)^2+s(9)^2+s(10)^2)^3))*s(8) + (3*mu*J2*(RE^2)/(2*(sqrt(s(8)^2+s(9)^2+s(10)^2)^5)))...
    *(5*((s(10)^2/(sqrt(s(8)^2+s(9)^2+s(10)^2)^2))-1))*s(8);
ds(12)= (-mu/(sqrt(s(8)^2+s(9)^2+s(10)^2)^3))*s(9) + (3*mu*J2*(RE^2)/(2*(sqrt(s(8)^2+s(9)^2+s(10)^2)^5)))...
    *(5*((s(10)^2/(sqrt(s(8)^2+s(9)^2+s(10)^2)^2))-1))*s(9);
ds(13)= (-mu/(sqrt(s(8)^2+s(9)^2+s(10)^2)^3))*s(10) + (3*mu*J2*(RE^2)/(2*(sqrt(s(8)^2+s(9)^2+s(10)^2)^5)))...
    *(5*((s(10)^2/(sqrt(s(8)^2+s(9)^2+s(10)^2)^2))-1)*s(10) - 2*s(10));
R =sqrt(s(8)^2+s(9)^2+s(10));
[r, p, y] = quat2angle(q);
GGTr      = (3*mu)*deltaIs*sin(2*r)/(2*R^3);
GGTr45    = (3*mu)*deltaIs*sin(2*(r+45*pi/180))/(2*R^3);
GGTp      = (3*mu)*deltaIs*sin(2*p)/(2*R^3);
GGTp45    = (3*mu)*deltaIs*sin(2*(p+45*pi/180))/(2*R^3);
ds(14)= s(15);
ds(15)= (-ks*s(14) - d*s(15) + GGTr)/Is;
ds(16)= s(17);
ds(17)= (-ks*s(16) - d*s(17) + GGTr45)/Is;
ds(18)= s(19);
ds(19)= (-ks*s(18) - d*s(19) + GGTp)/Is;
ds(20)= s(21);
ds(21)= (-ks*s(20) - d*s(21) + GGTp45)/Is;
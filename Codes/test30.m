clear
clc
deltaIs = 6.9*10^-8;                %  |I_y - I_x|
ks      = 1.62*10^-4;               %  spring stifness
time0   = [2018 09 14 12 00 00];    %  Satellite Start Mission
mu      = 398576.0576;              %  Gravitational Coe.
RE      = 6371;                     %  in km
alt     = 200;                      %  altitude in km
a0      = (alt+RE);                 %  a      = semi-major axis, km
e0      = 0;                        %  e      = eccentricity
i0      = 98.6;                     %  i      = inclination, deg
Raan0   = -15;                      %  Raan   = big omega = right ascension of the ascending node, deg
aop0    = 0;                        %  aop    = little omega = argument of periapse, deg                           
ta0     = 0;                        %  ta     = true anomaly, deg
p0      = a0*(1-e0^2);              %  p0
Torbit  = 2*pi*sqrt((a0^3)/mu);     %  Torbit = Period of Orbit
[R0,V0] = COE2RV(e0,i0,Raan0,aop0,ta0,p0,mu); % position Transformation from classic2R,v
s01=[6*pi/180 6*pi/180 0 1 0 0 0 R0(1,1) R0(2,1) R0(3,1) V0(1,1) V0(2,1) V0(3,1) 0 0 0 0 0 0 0 0]; % Initial
t01=datenum(time0):0.01:datenum(time0)+Torbit/600; % initial Time
[t,s]=ode23tb(@test29,t01,s01);
q = [s(:,4) s(:,5) s(:,6) s(:,7)];
R = sqrt(s(:,8).^2+s(:,9).^2+s(:,10).^2);
thetar = s(:,14); thetar45 = s(:,16) ; thetap = s(:,18); thetap45 = s(:,20);
[roll, pitch, yaw] = quat2angle(q);
GGTroll    = (3*mu)*deltaIs*sin(2*roll)/(2*R.^3);
GGTpitch    = (3*mu)*deltaIs*sin(2*pitch)/(2*R.^3);
r   = -0.5.*asin((-ks.*thetar   .*2.*R.^3)/(3*mu*deltaIs));
r45 = 0.5.*asin((-ks.*thetar45  .*2.*R.^3)/(3*mu*deltaIs));
p   = -0.5.*asin((-ks.*thetap   .*2.*R.^3)/(3*mu*deltaIs));
p45 = 0.5.*asin((-ks.*thetap45  .*2.*R.^3)/(3*mu*deltaIs));
%% Sensor tuning
for i=1:numel(t01)
    if pitch(i,1) >= 45*pi/180
        pitchm(i,1) = p45(i)+45*pi/180;
    elseif pitch(i,1) <= -45*pi/180
        pitchm(i,1) = -p45(i)-45*pi/180;  
    else
        pitchm(i,1) = p(i);
    end
end
for i=1:numel(t01)
    if pitch(i,1) >= 90*pi/180
        pitchm(i,1) = -p(i)+90*pi/180;
    elseif pitch(i,1) <= -90*pi/180
        pitchm(i,1) = -p(i)-90*pi/180;  
    else
        pitchm(i,1) = pitchm(i,1);
    end
end
for i=1:numel(t01)
    if pitch(i,1) >= 135*pi/180
        pitchm(i,1) = -p45(i)+135*pi/180;
    elseif pitch(i,1) <= -135*pi/180
        pitchm(i,1) = +p45(i)-135*pi/180;  
    else
        pitchm(i,1) = pitchm(i,1);
    end
end
for i=1:numel(t01)
    if roll(i,1) >= 45*pi/180
        rollm(i,1) = r45(i)+45*pi/180;
    elseif roll(i,1) <= -45*pi/180
        rollm(i,1) = -r45(i)-45*pi/180;  
    else
        rollm(i,1) = r(i);
    end
end
for i=1:numel(t01)
    if roll(i,1) >= 90*pi/180
        rollm(i,1) = -r(i)+90*pi/180;
    elseif roll(i,1) <= -90*pi/180
        rollm(i,1) = -r(i)-90*pi/180;  
    else
        rollm(i,1) = rollm(i,1);
    end
end
for i=1:numel(t01)
    if roll(i,1) >= 135*pi/180
        rollm(i,1) = -r45(i)+135*pi/180;
    elseif roll(i,1) <= -135*pi/180
        rollm(i,1) = +r45(i)-135*pi/180;  
    else
        rollm(i,1) = rollm(i,1);
    end
end
f0 = sqrt(ks/deltaIs)/(2*pi);
%% figures
% figure
% hold on
% plot3(s(:,8),s(:,9),s(:,10));
% grid on
% xlabel('X(km)');
% ylabel('Y(km)');
% zlabel('Z(km)')
% hold on
figure
subplot(4,1,1);
plot(t01-datenum(time0),roll*180/pi,t01-datenum(time0),rollm*180/pi);
hold on
xlabel(['time(sec)' '   ' 'frequency:' num2str(f0) '   '  'w_r= '  num2str(s(1,1)*180/pi) ' deg/sec'])
ylabel('roll(deg)');
grid on
legend('roll(deg)','roll sensor')
subplot(4,1,2);
plot(t01-datenum(time0),GGTroll(:,1));
hold on
xlabel('time(sec)');
ylabel('GGT_r(n.m)');
grid on
subplot(4,1,3);
plot(t01-datenum(time0),pitch*180/pi,t01-datenum(time0),pitchm*180/pi);
hold on
xlabel(['time(sec)' '   ' 'frequency: ' num2str(f0) '       '  'w_p= '  num2str(s(1,2)*180/pi) ' deg/sec']);
ylabel('pitch(deg)');
grid on
legend('pitch(deg)','pitch sensor')
subplot(4,1,4);
plot(t01-datenum(time0),GGTpitch(:,1));
hold on
xlabel('time(sec)');
ylabel('GGT_p (n.m)');
grid on
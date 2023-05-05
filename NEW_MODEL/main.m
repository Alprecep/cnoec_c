clear all
clc
close all
%%
%file:///C:/Users/alpre/OneDrive/Desktop/Polimi/Year5/CNOEC/cnoec_c/1-s2.0-S0959152414002704-main.pdf
% Mill and feeder parameters
alfa_f = 0.055;
alfa_r = 0.465;
alfa_P = 1;
alfa_speed = 0.712;
alfa_ksif = 0.01;
sigma_Ps = 0.5;
sigma_Pv = 0.5;
Db = 7.85;
Ds = 3.2;
epsilon_sv = 0.6;
phi_b = 90;
phi_f = 29.6;
phi_r = 6.03;
ksi_Pmax = 1662;
v_Pmax = 0.34;
Vv = 84;
Xp = 0;

%Cyclone parameters
alfa_su = 0.87;
C1 = 0.6;
C2 = 0.7;
C3 = 4;
C4 = 4;
epsilon_c = 129;

vmill = 100;
Pmax = 1662;
XP = 0;
V_v = 84;

alfa_phi_f = 0.01;

% initial values
CFF = 374;
MIW = 4.64;
MFB = 5.69;
MFS = 65.2;
SFW = 140.5;

Xmw = 4.85;
Xms = 4.9;
Xmf = 1.09;
Xmr = 1.82;
Xmb = 8.51;
Xsw = 4.11;
Xss = 1.88;
Xsf = 0.42;

JT = 0.34;
SVOL = 5.99;
PSE = 0.67;


phi_Pmax =0.57;
%
dt = 0.01;
ylist = [];
xlist = [];

for i=1:1e3
    %Rheology factor
    phi = sqrt(max(0, 1-(1/epsilon_sv)-1)*Xms/Xmw);
    %The effect of the total charge on mill power (Zx)
    LOAD = Xmw+Xmr+Xms+Xmb;
    Zx = LOAD/(vmill*v_Pmax)-1;
    %The effect of the slurry rheology on mill power (Zr) is given by the empirically defined equation:
    Zr = phi/phi_Pmax-1;
    %The mill power draw for the Hulbert-model is defined as:
    Pmill = Pmax*(1-sigma_Pv*Zx*Zx-2*XP*sigma_Pv*sigma_Ps*Zx*Zr-sigma_Ps*Zr*Zr)*(alfa_speed)^(alfa_P);


    % Feeder modulewater (Vfwo), solids (Vfso), fines (Vffo), rocks (Vfro) and balls (Vfbo). The flow-rates are defined as:
    %The mill outlet flow-rates for water (Vmwo), solids (Vmso), fines (Vmfo), rocks (Vmro) and balls (Vmbo)

    Vccu = CFF*(Xss-Xsf)/(Xsw+Xss)*(1-C1*exp(-CFF/epsilon_c))*...
        (1-(Xss/(C2*(Xsw+Xss)))^(C3))*...
        (1-(Xsf/Xss)^(C4));
    Fu = 0.6-(0.6-Xss/(Xsw+Xss))*exp(-Vccu/(alfa_su*epsilon_c));

    Vcwu = (Xsw*(Vccu-Fu*Vccu))/(Fu*Xsw+Fu*Xsf-Xsf);
    Vcfu = Xsf*(Vccu-Fu*Vccu)/(Fu*Xsw+Fu*Xsf-Xsf);
    Vcsu = Vccu + (Xsf*(Vccu-Fu*Vccu))/(Fu*Xsw+Fu*Xsf-Xsf);
    Vsso = CFF*(Xss/SVOL);
    Vsfo = CFF*(Xsf/SVOL);
    Vcso = Vsso - Vcsu;
    Vcfo = Vsfo - Vcfu;

    Xmwdot = MIW-phi*V_v*Xmw*Xmw/(Xms+Xmw)+Vcwu;
    Xmsdot = MFS/Ds*(1-alfa_r)-phi*V_v*Xmw*Xmw/(Xms+Xmw) + Vcsu + phi*Pmill/(Ds*phi_r)*(Xmr/(Xmr+Xms));
    Xmfdot = MFS/Ds*(alfa_r)  -phi*V_v*Xmw*Xmf/(Xms+Xmw) + Vcfu +   1*Pmill/(Ds*phi_f)/(1+alfa_phi_f*((Xmw+Xmr+Xms+Xmb)/vmill - v_Pmax));
    Xmrdot = MFS/Ds*(alfa_r)  -phi*Pmill/(Ds*phi_r)*(Xmr/(Xmr+Xms));
    Xmbdot = MFB/Db           -phi*Pmill/(phi_b)*(Xmb/(Ds*(Xmr+Xms)+Db*Xmb));
    Xswdot = phi*Vv*Xmw*Xmw/(Xms+Xmw)-CFF*Xsw/(Xsw+Xss)+SFW;
    Xssdot = phi*Vv*Xmw*Xms/(Xms+Xmw)-CFF*Xss/(Xsw+Xss);
    Xsfdot = phi*Vv*Xmw*Xmf/(Xms+Xmw)-CFF*Xsf/(Xsw+Xss);

    Xmw = Xmw + Xmwdot*dt;
    Xms = Xms + Xmsdot*dt;
    Xmf = Xmf + Xmfdot*dt;
    Xmr = Xmr + Xmrdot*dt;
    Xmb = Xmb + Xmbdot*dt;
    Xsw = Xsw + Xswdot*dt;
    Xss = Xss + Xssdot*dt;
    Xsf = Xsf + Xsfdot*dt;

    xlist(i,:) = [Xmw;Xms;Xmf;Xmr;Xmb;Xsw;Xss;Xsf];
    ylist(i,:) = [JT;SVOL;PSE]';
    JT = (Xmw+Xms+Xmr+Xmb)/(vmill);
    %SVOL is the volume of slurry in the sump:
    SVOL = Xss+Xsw;
    PSE = Vcfo/(Vcso);
end


%%
plot(ylist,"rx:")
%%
plot(xlist,"k:")


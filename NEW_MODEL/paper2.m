clc
clear all
%% PARAMETERS
% Mill and feeder parameters
alfa_f = 0.05;
alfa_r = 0.47;
alfa_P = 1;
alfa_speed = 0.71;
alfa_psi_f = 0.01;

sigma_Ps = 0.5;
sigma_Pv = 0.5;

Db = 7.85;
Ds = 3.2;

epsilon_c =129;
epsilon_sv = 0.6;
phi_b = 90;
phi_f = 29.5;
phi_r = 6.00;

Pmax = 1662;
vmill = 100;
v_Pmax = 0.34;
Vv = 84;
Xp = 0;

%Cyclone parameters
alfa_su = 0.87;
C1 = 0.6;
C2 = 0.7;
C3 = 4;
C4 = 4;

% States
Xmb = 8.51;
Xmf = 0.6;
Xmr = 1.82;
Xms = 4.90;
Xmw = 4.85;
Xsf = 0.42;
Xss = 1.88;
Xsw = 4.11;

% Feeder modulewater (Vfwo), solids (Vfso), fines (Vffo), rocks (Vfro) and balls (Vfbo). The flow-rates are defined as:
MIW = 1;
MFS = 1;
MFB = 1;

Vfwo = MIW;
Vfso = MFS/Ds*(1-alfa_r);
Vffo = MFS/Ds*alfa_f;
Vfro = MFS/Ds*alfa_r;
Vfbo = MFB/Db;

% The population volume balance of the hold-up of water (Xmw),
%solids (Xms), fines (Xmf), rocks (Xmr) and balls (Xmb) in the mill are
%defined in terms of the inflow and outflow of each state:

Xmwdt = Vmwi - Vmwo ;
Xmsdt = Vmsi - Vmso + RC;
Xmfdt = Vmfi - Vmfo + FP;
Xmrdt = Vmri - RC;
Xmbdt = Vmbi - BC;

%The mill outlet flow-rates for water (Vmwo), solids (Vmso), fines (Vmfo), rocks (Vmro) and balls (Vmbo)
Vmwi = Vfwo+Vcwu;
Vmsi = Vfso+Vcsu;
Vmfi = Vffo+Vcfu;
Vmri = Vfro;
Vmbi = Vfbo;

%Rheology factor
ksi = sqrt(max(0, 1-(1/epsilon_sv)-1)*Xms/Xmw);

%The mill power draw for the Hulbert-model is defined as:
Pmill = Pmax*(1-sigma_Pv*Zx*Zx-2*XP*sigma_Pv*sigma_Ps*Zx*Zr-sigma_Ps*Zr*Zr)*(alfa_speed)^(alfa_P);

%The effect of the total charge on mill power (Zx)
Zx = LOAD/(v_mill*v_Pmax)-1;
LOAD = Xmw+Xmr+Xms+Xmb;

%The effect of the slurry rheology on mill power (Zr) is given by the empirically defined equation:
Zr = ksi/ksi_Pmax-1;

%RC r rock consumption is:
RC = Pmil*ksi/(Ds*ksi_r)*(Xmr/(Xmr+Xms));
%Ball consumption is defined as:
BC = Pmill*ksi/(phi_b)*(Xmb/(Ds*(Xmr+Xms)+Db*Xmb));
%The production of fines in the mill is defined as:
FP = (Pmill)/(Ds*phi_f*(1+alfa_phif)*(LOAD/vmill-v_Pmax));

%The discharge flow-rates of water (Vmwo), solids (Vmso), fines (Vmfo), rocks (Vmro) and balls (Vmbo) through the end-discharge
%screen are defined as:
Vmow = Vv*ksi*Xmw*(Xmw/(Xms+Xmw));
Vmso = Vv*ksi*Xmw*(Xms/(Xms+Xmw));
Vmfo = Vv*ksi*Xmw*(Xmwf(Xms+Xmw));
Vmro = 0;
Vmbo = 0;

%The population volume balance of the hold-up of water (Xsw), solids (Xss) and fines (Xsf) in the sump are defined as:
Xswdt = Vswi - Vswo + SFW;
Xssdt = Vssi - Vsso;
Xsfdt = xsfi - Vsfo;

%The sump discharge flow-rates of water (Vswo), solids (Vswo) and fines (Vsfo) are defined as:
Vswo = CFF * (Xsw/SVOL);
Vsso = CFF * (Xss/SVOL);
Vsfo = CFF * (Xsf/SVOL);

%SVOL is the volume of slurry in the sump:
SVOL = Xsw + Xss;


%where Fi is the fraction of solids in the total inflow volume and CFF is the cyclone feed flow-rate (
Fi = Vcsi / CFF;


%Thus, the sensitivity of the coarse material split to the total cyclone feed flow-rate is given by:
VccuVcci_1 = 1-C1*exp(-CFF/epsilon_c);
VccuVcci_2 = 1-(Fi/C2)^(C3);

%fraction fines in the feed solids (Pi). In terms of the inflow of fines
%(Vcfi) and solids (Vcsi), Pi is defined as:
Pi = Vcfi / Vcsi;

VccuVcci_3 = 1-Pi^C4;
VccuVcci = VccuVcci_1*VccuVcci_2*VccuVcci_3;

Fu = 0.6-(0.6-Fi)*exp(-Vccu/(alfa_su*epsilon_c));
Fu = (Vcfu+Vccu)/(Vcwu+Vcfu+Vccu);

Vcwu = Vcwi*(Vccu-Fu*Vccu)/(Fu*Vcwi+Fu*Vcfi-Vcfi);
Vcfu = Vcfi*(Vccu-Fu*Vccu)/(Fu*Vcwi+Fu*Vcfi-Vcfi);

Vcwi = Vswo;
Vcci = Vsso-Vsfo;
Vcsi = Vsso;
Vcfi = Vsfo;

PSE = Vcfo/(Vcco+Vcfo) ;
PSE = (Vcfi-Vcfu)/(Vcci-Vccu+Vcfo);


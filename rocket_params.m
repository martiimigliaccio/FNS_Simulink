%% rocket_params.m — Parametri geometrici e di massa

%% Geometria
L_rocket  = 1.00;        % lunghezza totale [m]
D         = 0.11;        % diametro [m]
R         = D/2;         % raggio [m]
L_nose    = 0.15;        % lunghezza ogiva ellittica [m]
A_ref     = pi*R^2;      % area di riferimento [m²]

%% Fin (Clipped Delta, 3 alette)
n_fins    = 3;
cr        = 1.5*D;       % root chord = 1.5D = 0.165 m
span_fin  = 1.1*D;       % span = 1.1D = 0.121 m
taper     = 0.5;         % taper ratio (ct/cr)
ct        = taper*cr;    % tip chord = 0.0825 m
t_fin     = 0.003;       % spessore [m]

%% Masse
m_total_0 = 1.945;       % massa iniziale totale [kg] (razzo + acqua piena)
% Tank propulsiva: 0.5 L acqua = 0.5 kg
% Tank payload:   0.5 L acqua = 0.5 kg
m_water_prop    = 0.5;   % [kg]
m_water_payload = 0.5;   % [kg]
m_dry           = m_total_0 - m_water_prop - m_water_payload; % massa a secco [kg]

%% CG a secco (da calcolare o misurare — stima geometrica)
% Stima: CG a secco ≈ 50% della lunghezza dalla punta (aggiornare con misura reale)
xcg_dry = 0.50;  % [m] dalla punta

%% Paracadute
D_chute   = 1.0;          % diametro [m]
Cd_chute  = 0.8;
A_chute   = pi*(D_chute/2)^2;
% Forza paracadute: F_chute = 0.5*rho*v^2 * Cd_chute * A_chute

%% Atmosfera (ISA semplificata)
rho0  = 1.225;  % [kg/m³] densità aria a livello del mare
g0    = 9.81;   % [m/s²]

load('neural_slosh_model.mat')
%% generate_slosh_data.m
% Genera dati di training per il Neural State Space dello sloshing
% Parametri calibrati sul volo reale (Axb picco ~550 m/s², Azb range ±2 m/s²)
% Durata: 10s — ~9000 campioni dinamicamente ricchi

clear; clc;

%% --------------------------------------------------------
%  1. PARAMETRI PENDOLO EQUIVALENTE (NASA SP-106, Config A)
%% --------------------------------------------------------
m_s  = 0.130;       % kg  massa oscillante
m_f  = 0.365;       % kg  massa fissa
fn   = 3.18;        % Hz  frequenza naturale
wn   = 2*pi*fn;     % rad/s
zeta = 0.03;        % smorzamento (tipico acqua in cilindro)
L_p  = 9.81/wn^2;  % m   lunghezza pendolo equivalente

fprintf('=== Generazione dati sloshing ===\n')
fprintf('L_p  = %.5f m\n', L_p)
fprintf('wn   = %.4f rad/s\n', wn)
fprintf('zeta = %.3f\n', zeta)

%% --------------------------------------------------------
%  2. GEOMETRIA TANK
%% --------------------------------------------------------
xcg_tank_pay = 0.35;  % m dalla punta del razzo

%% --------------------------------------------------------
%  3. GRIGLIA TEMPORALE
%% --------------------------------------------------------
dt    = 0.001;   % s
T_tot = 10;      % s  — bilanciamento velocità training / ricchezza dati
t     = (0:dt:T_tot)';
N     = length(t);

%% --------------------------------------------------------
%  4. PROFILO Axb CALIBRATO SUL VOLO REALE
%% --------------------------------------------------------
Axb = zeros(N,1);

% Fase propulsiva (0 - 0.083s): picco 550 m/s² decrescente
t_burn   = 0.083;
idx_burn = t <= t_burn;
Axb(idx_burn) = 550 * (1 - t(idx_burn)/t_burn);

% Coasting (0.083s - 1.1s): decelerazione drag
idx_coast = t > t_burn & t <= 1.1;
Axb(idx_coast) = -10 * exp(-(t(idx_coast) - t_burn)/0.4);

% Discesa (1.1s - fine): gravità dominante
idx_down = t > 1.1;
Axb(idx_down) = -9.81 * (1 - 0.3*exp(-(t(idx_down)-1.1)/5));

%% --------------------------------------------------------
%  5. PROFILO Azb CALIBRATO SUL VOLO REALE (±2 m/s²)
%% --------------------------------------------------------
rng(42)  % riproducibilità

Azb = 1.2*sin(2*pi*0.4*t)      ...   % oscillazione lenta di assetto
    + 0.6*sin(2*pi*1.1*t + 0.8) ...  % oscillazione media
    + 0.3*sin(2*pi*2.3*t + 1.5) ...  % oscillazione veloce
    + 0.15*randn(N,1);               % rumore

% Perturbazione iniziale durante il burn
Azb(idx_burn) = Azb(idx_burn) + 1.5*sin(2*pi*3*t(idx_burn));

% Smussamento ed clipping
Azb = smoothdata(Azb, 'gaussian', 30);
Azb = max(min(Azb, 2.0), -2.0);

%% --------------------------------------------------------
%  6. INTEGRAZIONE PENDOLO EQUIVALENTE
%  Equazione: δ̈ + 2ζωn·δ̇ + ωn²·δ = -Azb/L_p
%  Termine Axb·δ/L_p RIMOSSO: instabile con Axb >> g
%% --------------------------------------------------------
Axb_i = @(tt) interp1(t, Axb, tt, 'linear', 0);
Azb_i = @(tt) interp1(t, Azb, tt, 'linear', 0);

odefun = @(tt, y) [
    y(2);
    -2*zeta*wn*y(2) - wn^2*y(1) - Azb_i(tt)/L_p
];

fprintf('Integrazione ODE pendolo...\n')
opts = odeset('RelTol',1e-6, 'AbsTol',1e-8, 'MaxStep', dt);
[t_ode, y_ode] = ode45(odefun, [0, T_tot], [0; 0], opts);

% Ricampiona su griglia uniforme
delta      = interp1(t_ode, y_ode(:,1), t, 'linear');
delta_dot  = interp1(t_ode, y_ode(:,2), t, 'linear');
delta_ddot = gradient(delta_dot, dt);

%% --------------------------------------------------------
%  7. CALCOLO FORZE E MOMENTI DI SLOSHING
%% --------------------------------------------------------
Fz_sl = -m_s * (L_p * delta_ddot + Azb);         % forza laterale [N]
Fx_sl = -m_s * L_p * delta_ddot .* delta;         % forza assiale (2° ordine) [N]
My_sl = -m_s * Azb * xcg_tank_pay;               % momento di pitch [N·m]

fprintf('Integrazione completata.\n')
fprintf('delta:    max=%.4f rad = %.2f deg\n', max(abs(delta)), max(abs(delta))*180/pi)
fprintf('|Fz_sl|:  max=%.4f N\n',   max(abs(Fz_sl)))
fprintf('|Fx_sl|:  max=%.4f N\n',   max(abs(Fx_sl)))
fprintf('|My_sl|:  max=%.4f N·m\n', max(abs(My_sl)))

%% --------------------------------------------------------
%  8. NORMALIZZAZIONE Z-SCORE
%% --------------------------------------------------------
mu_in  = mean([Axb, Azb]);
std_in = std([Axb, Azb]);
std_in(std_in < 1e-6) = 1;

mu_out  = mean([Fz_sl, Fx_sl, My_sl]);
std_out = std([Fz_sl, Fx_sl, My_sl]);
std_out(std_out < 1e-6) = 1;

input_data  = ([Axb, Azb]            - mu_in)  ./ std_in;
output_data = ([Fz_sl, Fx_sl, My_sl] - mu_out) ./ std_out;
state_data  = [delta, delta_dot];

fprintf('\nNormalizzazione:\n')
fprintf('  Input  mu =[%.3f, %.4f]  std=[%.3f, %.4f]\n', ...
    mu_in(1),  mu_in(2),  std_in(1),  std_in(2))
fprintf('  Output mu =[%.5f, %.6f, %.5f]\n', mu_out(1), mu_out(2), mu_out(3))
fprintf('  Output std=[%.5f, %.6f, %.5f]\n', std_out(1), std_out(2), std_out(3))

%% --------------------------------------------------------
%  9. SALVA
%% --------------------------------------------------------
save('slosh_training_data.mat', ...
     't', 'dt', 'input_data', 'output_data', 'state_data', ...
     'm_s', 'm_f', 'wn', 'zeta', 'L_p', ...
     'mu_in', 'std_in', 'mu_out', 'std_out')

fprintf('\nDataset salvato: %d campioni, T_tot=%ds, dt=%.3fs\n', N, T_tot, dt)

%% --------------------------------------------------------
%  10. PLOT DI VERIFICA
%% --------------------------------------------------------
figure('Name','Dati training sloshing')
subplot(3,2,1); plot(t, Axb);          title('Axb [m/s²]');   grid on; xlabel('t [s]')
subplot(3,2,2); plot(t, Azb);          title('Azb [m/s²]');   grid on; xlabel('t [s]')
subplot(3,2,3); plot(t, delta*180/pi); title('δ [°]');         grid on; xlabel('t [s]')
subplot(3,2,4); plot(t, Fz_sl);        title('Fz_{sl} [N]');  grid on; xlabel('t [s]')
subplot(3,2,5); plot(t, Fx_sl);        title('Fx_{sl} [N]');  grid on; xlabel('t [s]')
subplot(3,2,6); plot(t, My_sl);        title('My_{sl} [N·m]');grid on; xlabel('t [s]')
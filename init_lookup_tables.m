%% init_lookup_tables.m
% Carica e prepara tutte le lookup table per il modello razzo

%% --- Cd (Drag Coefficient) ---
data_cd = readtable('cdfns.csv');
t_cd    = data_cd{:,1};
Cd_val  = data_cd{:,2};
Cd_val(isnan(Cd_val)) = 0.315;

%% --- Thrust dal profilo pulito ---
data_thr = readtable('thrust_clean.csv');
t_thr    = data_thr{:,1};   % colonna time_s
T_val    = data_thr{:,2};   % colonna thrust_N
T_val(isnan(T_val)) = 0;

LUT_Thrust = [t_thr, T_val];
fprintf('Thrust: %d punti, picco=%.1f N, burnout=%.4f s\n', ...
    length(t_thr), max(T_val), t_thr(find(T_val > 0, 1, 'last')))

%% --- Cm (Pitch Moment Coefficient) ---
data_cm = readtable('cmfns.csv');
t_cm    = data_cm{:,1};
Cm_val  = data_cm{:,2};
Cm_val(isnan(Cm_val)) = 0;

%% --- Cp (Center of Pressure position, in cm dalla punta) ---
data_cp = readtable('cpfns.csv');
t_cp    = data_cp{:,1};
Cp_val  = data_cp{:,2} / 100;
first_valid_cp = Cp_val(find(~isnan(Cp_val), 1));
Cp_val(isnan(Cp_val)) = first_valid_cp;

%% --- Fn (Normal Force Coefficient) ---
data_fn = readtable('fnfns.csv');
t_fn    = data_fn{:,1};
Fn_val  = data_fn{:,2};
Fn_val(isnan(Fn_val)) = 0;

%% --- Fl (Lateral Force Coefficient) ---
data_fl = readtable('flfns.csv');
t_fl    = data_fl{:,1};
Fl_val  = data_fl{:,2};
Fl_val(isnan(Fl_val)) = 0;

%% --- Packaging come vettori [t, val] per i Lookup Table blocks ---
LUT_Cd     = [t_cd,  Cd_val];
LUT_Thrust = [t_thr, T_val];
LUT_Cm     = [t_cm,  Cm_val];
LUT_Cp     = [t_cp,  Cp_val];
LUT_Fn     = [t_fn,  Fn_val];
LUT_Fl     = [t_fl,  Fl_val];

disp('Lookup tables caricate correttamente.');
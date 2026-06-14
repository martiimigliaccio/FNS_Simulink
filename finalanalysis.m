%% rocket_competition_analysis.m
% =========================================================================
%  WATER ROCKET — COMPLETE COMPETITION ANALYSIS TOOL
%  Airbus Sloshing Rocket Competition
% =========================================================================
%  INPUTS  (edit Section 1):
%    - Mass budget (wet & dry per component)
%    - Aerodynamic coefficients (LUT from workspace)
%    - Propulsion curve (LUT from workspace + pressure scaling)
%    - Launch pressure [bar absolute]
%
%  OUTPUTS (9 individual figures + Excel file):
%    Fig 1  — 2D Trajectory with orientation vectors
%    Fig 2  — Altitude vs Time
%    Fig 3  — Velocity vs Time
%    Fig 4  — Flight path angle and angle of attack
%    Fig 5  — CG and CP positions [mm from nose] vs time
%    Fig 6  — Static margin [calibres] vs time
%    Fig 7  — Total mass and thrust vs time
%    Fig 8  — Competition score breakdown
%    Fig 9  — Mass budget table
%    Excel  — rocket_competition_results.xlsx (4 sheets)
%
%  REQUIREMENTS: rocket_params + init_lookup_tables already run
% =========================================================================

clear; clc; close all;

if ~exist('LUT_Thrust','var')
    fprintf('[INFO] Loading workspace...\n')
    rocket_params;
    init_lookup_tables;
end

%% =========================================================================
%% SECTION 1 — INPUT PARAMETERS (edit here)
%% =========================================================================

%% 1.1 — LAUNCH PRESSURE
P_launch_bar = 10;   % absolute launch pressure [bar]
P_ref_bar    = 4.0;   % reference pressure of ID18 profile [bar]
P_atm_bar    = 1.0;   % atmospheric pressure [bar]

dP_launch  = P_launch_bar - P_atm_bar;
dP_ref     = P_ref_bar    - P_atm_bar;
T_scale    = dP_launch / dP_ref;
mdot_scale = sqrt(dP_launch / dP_ref);
Isp_scale  = T_scale / mdot_scale;

%% 1.2 — MASS BUDGET  {name, wet_g, dry_g}
mass_budget = {
    'Body structure',          550,  550;
    'Elliptical nose cone',     80,   80;
    'Fin set (3x)',             60,   60;
    'Propulsive tank (body)',  155,  155;
    'Propulsive water',        500,    0;
    'Payload tank (body)',     150,  150;
    'Payload water',           500,  500;
    'Parachute + housing',      80,   80;
    'Electronics / misc',       20,   20;
};

wet_vals      = cell2mat(mass_budget(:,2));
dry_vals      = cell2mat(mass_budget(:,3));
delta_vals    = wet_vals - dry_vals;
m_TOW_g       = sum(wet_vals);
m_dry_total_g = sum(dry_vals);
m_payload_g   = cell2mat(mass_budget(7,3));

fprintf('[CHECK] m_TOW = %.0f g  (expected 1945 g)\n', m_TOW_g)

%% 1.3 — PHYSICAL PARAMETERS
p.g0            = 9.81;
p.rho_air       = 1.225;
p.m_dry         = m_dry_total_g/1000 - 0.5;
p.m_wpay        = 0.500;
p.m_wprop0      = 0.500;
p.R_rocket      = 0.055;
p.A_ref         = pi * p.R_rocket^2;
p.D_ref         = 0.11;
p.L_rocket      = 1.00;
p.L_tank        = 0.30;
p.xcg_dry       = 0.50;
p.xcg_tank_prop = 0.85;
p.xcg_tank_pay  = 0.35;
p.Isp_eff       = 7.5 * Isp_scale;
p.Cd_chute      = 0.8;
p.A_chute       = pi * 0.5^2;
p.gamma0_deg    = 88;

%% 1.4 — LOOKUP TABLES (scale thrust for pressure)
T_scaled      = LUT_Thrust;
T_scaled(:,2) = LUT_Thrust(:,2) * T_scale;
p.LUT_T       = T_scaled;
p.LUT_Cd      = LUT_Cd;
p.LUT_Cm      = LUT_Cm;
p.LUT_Cp      = LUT_Cp;
p.LUT_Fn      = LUT_Fn;

fprintf('[INFO] Pressure scaling: T x%.3f, Isp x%.3f\n', T_scale, Isp_scale)
fprintf('[INFO] Peak thrust: %.1f N | Burnout: %.4f s\n', ...
    max(p.LUT_T(:,2)), p.LUT_T(find(p.LUT_T(:,2)>0,1,'last'),1))

%% =========================================================================
%% SECTION 2 — ODE SIMULATION
%% =========================================================================
gamma0_rad = p.gamma0_deg * pi/180;
y0         = [0; 0; 0; 0; gamma0_rad; 0; p.m_wprop0];
t_span     = [0, 60];
opts       = odeset('RelTol',1e-5,'AbsTol',1e-7,'MaxStep',5e-4,...
                    'Events',@(t,y) ground_event(t,y));

fprintf('[SIM] Running simulation...\n')
[t_sim, Y_sim] = ode45(@(t,y) rocket_ode_3dof(t,y,p), t_span, y0, opts);
fprintf('[SIM] Done: %d points, t_final=%.2f s\n', length(t_sim), t_sim(end))

%% Extract state variables
xe_sim    = Y_sim(:,1);  Vxe_sim = Y_sim(:,2);
ze_sim    = Y_sim(:,3);  Vze_sim = Y_sim(:,4);
th_sim    = Y_sim(:,5);  q_sim   = Y_sim(:,6);
mp_sim    = Y_sim(:,7);
V_sim     = sqrt(Vxe_sim.^2 + Vze_sim.^2);
gam_sim   = atan2(Vze_sim, Vxe_sim);
alpha_sim = th_sim - gam_sim;

%% Key metrics
[h_max, idx_ap] = max(ze_sim);
t_apogeo        = t_sim(idx_ap);
d_final         = xe_sim(end);
t_final         = t_sim(end);
V_max           = max(V_sim);
t_burnout_val   = p.LUT_T(find(p.LUT_T(:,2)>0,1,'last'),1);
idx_burn        = find(t_sim >= t_burnout_val, 1);
V_at_burnout    = V_sim(idx_burn);

%% Competition score
m_pay_kg = p.m_wpay;
m_TOW_kg = m_TOW_g / 1000;
score    = (sqrt(d_final^2 + h_max^2) + t_final) * (m_pay_kg / m_TOW_kg);

%% =========================================================================
%% SECTION 3 — STABILITY (CG, CP, Static Margin)
%% =========================================================================
xcg_sim = zeros(size(t_sim));
for i = 1:length(t_sim)
    xcg_sim(i) = cg_calc(mp_sim(i), p);
end

Cp_sim = zeros(size(t_sim));
for i = 1:length(t_sim)
    tc = min(t_sim(i), p.LUT_Cp(end,1));
    Cp_sim(i) = interp1(p.LUT_Cp(:,1), p.LUT_Cp(:,2), ...
        tc, 'linear', p.LUT_Cp(end,2));
end

SM_sim = (Cp_sim - xcg_sim) / p.D_ref;

%% =========================================================================
%% SECTION 4 — PRINT RESULTS
%% =========================================================================
if all(SM_sim > 0)
    stable_str = 'YES — SM > 0 throughout flight';
else
    stable_str = 'NO  — unstable phase detected, check CG/CP';
end

fprintf('\n')
fprintf('=============================================================\n')
fprintf('       WATER ROCKET — COMPETITION RESULTS SUMMARY            \n')
fprintf('=============================================================\n')
fprintf('  TRAJECTORY GEOMETRY\n')
fprintf('    h_max    = %7.2f m\n',  h_max)
fprintf('    d_final  = %7.2f m\n',  d_final)
fprintf('    t_flight = %7.2f s\n',  t_final)
fprintf('    t_apogee = %7.2f s\n',  t_apogeo)
fprintf('-------------------------------------------------------------\n')
fprintf('  PROPULSION  (P_launch = %.1f bar)\n', P_launch_bar)
fprintf('    T_max    = %7.1f N\n',  max(p.LUT_T(:,2)))
fprintf('    t_burn   = %7.4f s\n',  t_burnout_val)
fprintf('    V_burn   = %7.2f m/s\n', V_at_burnout)
fprintf('    V_max    = %7.2f m/s\n', V_max)
fprintf('-------------------------------------------------------------\n')
fprintf('  MASS BUDGET\n')
fprintf('    m_TOW    = %7.0f g\n',  m_TOW_g)
fprintf('    m_dry    = %7.0f g\n',  m_dry_total_g)
fprintf('    m_payload= %7.0f g  (payload water)\n', m_payload_g)
fprintf('    pay/TOW  = %7.4f\n',    m_pay_kg/m_TOW_kg)
fprintf('-------------------------------------------------------------\n')
fprintf('  STABILITY\n')
fprintf('    SM_min   = %7.3f cal\n', min(SM_sim))
fprintf('    SM_max   = %7.3f cal\n', max(SM_sim))
fprintf('    Stable?  = %s\n',        stable_str)
fprintf('-------------------------------------------------------------\n')
fprintf('  COMPETITION SCORE\n')
fprintf('    formula  = (sqrt(d^2+h^2) + t) x (m_pay / m_TOW)\n')
fprintf('    score    = (sqrt(%.2f^2+%.2f^2) + %.2f) x %.4f\n', ...
    d_final, h_max, t_final, m_pay_kg/m_TOW_kg)
fprintf('    SCORE    = %.4f\n', score)
fprintf('=============================================================\n\n')

%% =========================================================================
%% SECTION 5 — FIGURES
%% =========================================================================
c1  = [0.00 0.45 0.74];
c2  = [0.85 0.33 0.10];
c3  = [0.47 0.67 0.19];
fw  = 680; fh = 480;

%% helper: draw phase lines without polluting legend
    function draw_phase_lines(tb, tap)
        xline(tb,  'k:',  'HandleVisibility','off');
        xline(tap, 'k--', 'HandleVisibility','off');
        yl = ylim;
        text(tb,  yl(2), ' Burnout', 'FontSize',7, ...
            'VerticalAlignment','top','Color',[0.3 0.3 0.3])
        text(tap, yl(2), ' Apogee',  'FontSize',7, ...
            'VerticalAlignment','top','Color',[0.3 0.3 0.3])
    end

%% -- Figure 1: 2D Trajectory --
figure('Name','Fig 1 — 2D Trajectory','NumberTitle','off', ...
    'Position',[50 550 fw fh])
h_traj = plot(xe_sim, ze_sim, '-','Color',c1,'LineWidth',2, ...
    'DisplayName','Trajectory');
hold on
N = max(1, floor(length(t_sim)/18));
idx_arr = 1:N:length(t_sim);
h_ori = quiver(xe_sim(idx_arr), ze_sim(idx_arr), ...
    cos(gam_sim(idx_arr)), sin(gam_sim(idx_arr)), ...
    0.4, 'r','LineWidth',1,'DisplayName','Orientation');
h_lau = plot(xe_sim(1),      ze_sim(1),      'go', ...
    'MarkerSize',10,'MarkerFaceColor','g','DisplayName','Launch');
h_apo = plot(xe_sim(idx_ap), ze_sim(idx_ap), 'k^', ...
    'MarkerSize',10,'MarkerFaceColor','y','DisplayName','Apogee');
h_lan = plot(xe_sim(end),    ze_sim(end),    'rs', ...
    'MarkerSize',10,'MarkerFaceColor','r','DisplayName','Landing');
text(xe_sim(idx_ap)+0.05, ze_sim(idx_ap)+0.3, ...
    sprintf('h = %.2f m', h_max), 'FontSize',9)
text(xe_sim(end)+0.05, 0.3, ...
    sprintf('d = %.2f m\nt = %.2f s', d_final, t_final), 'FontSize',9)
xlabel('Horizontal distance x_e [m]')
ylabel('Altitude [m]')
title(sprintf('2D Trajectory  |  h_{max}=%.2f m  |  d=%.2f m  |  t=%.2f s', ...
    h_max, d_final, t_final))
legend([h_traj h_ori h_lau h_apo h_lan], 'Location','best')
grid on

%% -- Figure 2: Altitude vs Time --
figure('Name','Fig 2 — Altitude vs Time','NumberTitle','off', ...
    'Position',[100 500 fw fh])
plot(t_sim, ze_sim, '-','Color',c1,'LineWidth',2)
hold on
draw_phase_lines(t_burnout_val, t_apogeo)
xlabel('Time [s]'); ylabel('Altitude [m]')
title(sprintf('Altitude vs Time  |  h_{max}=%.2f m  at  t=%.2f s', ...
    h_max, t_apogeo))
grid on
text(0.05, 0.90, sprintf('h_{max} = %.2f m\nt_{apogee} = %.2f s', ...
    h_max, t_apogeo), 'Units','normalized','FontSize',9, ...
    'BackgroundColor','w','EdgeColor','k')

%% -- Figure 3: Velocity vs Time --
figure('Name','Fig 3 — Velocity vs Time','NumberTitle','off', ...
    'Position',[150 450 fw fh])
h1 = plot(t_sim, V_sim,   '-', 'Color',c2,'LineWidth',2, ...
    'DisplayName','V_{total}');
hold on
h2 = plot(t_sim, Vze_sim, '--','Color',c3,'LineWidth',1.5, ...
    'DisplayName','V_z (vertical)');
h3 = plot(t_sim, Vxe_sim, ':','Color',c1,'LineWidth',1.2, ...
    'DisplayName','V_x (horizontal)');
yline(0,'k--','HandleVisibility','off','LineWidth',0.8)
draw_phase_lines(t_burnout_val, t_apogeo)
xlabel('Time [s]'); ylabel('Velocity [m/s]')
title(sprintf('Velocity vs Time  |  V_{max}=%.2f m/s  |  V_{burnout}=%.2f m/s', ...
    V_max, V_at_burnout))
legend([h1 h2 h3], 'Location','best')
grid on

%% -- Figure 4: Gamma and Alpha --
figure('Name','Fig 4 — Gamma and Alpha','NumberTitle','off', ...
    'Position',[200 400 fw fh])
h1 = plot(t_sim, rad2deg(gam_sim),   '-', 'Color',c1,'LineWidth',2, ...
    'DisplayName','\gamma — Flight path angle');
hold on
h2 = plot(t_sim, rad2deg(alpha_sim), '--','Color',c2,'LineWidth',1.5, ...
    'DisplayName','\alpha — Angle of attack');
h3 = plot(t_sim, rad2deg(th_sim),    ':','Color',[0.5 0 0.5],'LineWidth',1.2, ...
    'DisplayName','\theta — Pitch angle');
yline(0,'k--','HandleVisibility','off','LineWidth',0.8)
draw_phase_lines(t_burnout_val, t_apogeo)
xlabel('Time [s]'); ylabel('Angle [deg]')
title('Flight Path Angle, Angle of Attack and Pitch Angle vs Time')
legend([h1 h2 h3], 'Location','best')
grid on

%% -- Figure 5: CG and CP vs Time --
figure('Name','Fig 5 — CG and CP Positions','NumberTitle','off', ...
    'Position',[250 350 fw fh])
yyaxis left
fill([t_sim; flipud(t_sim)], ...
    [xcg_sim*1000; flipud(Cp_sim*1000)], ...
    c3,'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off')
hold on
h1 = plot(t_sim, xcg_sim*1000, '-', 'Color',c2,'LineWidth',2, ...
    'DisplayName','CG');
h2 = plot(t_sim, Cp_sim*1000,  '--','Color',c1,'LineWidth',2, ...
    'DisplayName','CP');
draw_phase_lines(t_burnout_val, t_apogeo)
ylabel('Position from nose [mm]')
xlabel('Time [s]')
yyaxis right
h3 = plot(t_sim, SM_sim, ':','Color',[0.4 0 0.6],'LineWidth',1.5, ...
    'DisplayName','Static margin [cal]');
yline(0,'r--','HandleVisibility','off','LineWidth',1)
ylabel('Static margin [calibres]')
title('CG and CP Positions vs Time — Stability Overview')
legend([h1 h2 h3], 'Location','best')
grid on

%% -- Figure 6: Static Margin --
figure('Name','Fig 6 — Static Margin','NumberTitle','off', ...
    'Position',[300 300 fw fh])
area_pos = SM_sim; area_pos(SM_sim < 0) = 0;
area_neg = SM_sim; area_neg(SM_sim > 0) = 0;
area(t_sim, area_pos,'FaceColor',c3,'FaceAlpha',0.35,'EdgeColor','none', ...
    'HandleVisibility','off')
hold on
area(t_sim, area_neg,'FaceColor','r','FaceAlpha',0.45,'EdgeColor','none', ...
    'HandleVisibility','off')
plot(t_sim, SM_sim, '-','Color',[0.15 0.15 0.15],'LineWidth',2, ...
    'HandleVisibility','off')
yline(1.0,'--','Color',c1,'HandleVisibility','off')
yline(0,  'r--','HandleVisibility','off')
draw_phase_lines(t_burnout_val, t_apogeo)
xlabel('Time [s]'); ylabel('Static Margin [calibres]')
title(sprintf('Static Margin vs Time  |  min=%.3f cal  |  max=%.3f cal', ...
    min(SM_sim), max(SM_sim)))
ylim([min(SM_sim)-0.5, max(SM_sim)+0.5])
grid on
if all(SM_sim > 0)
    stab_txt = 'STABLE throughout flight';
    box_col  = [0.88 1.00 0.88];
else
    stab_txt = 'UNSTABLE phase detected!';
    box_col  = [1.00 0.80 0.80];
end
text(0.02, 0.93, sprintf('SM_{min} = %.3f cal   %s', min(SM_sim), stab_txt), ...
    'Units','normalized','FontSize',9,'FontWeight','bold', ...
    'BackgroundColor',box_col,'EdgeColor','k')

%% -- Figure 7: Mass and Thrust --
figure('Name','Fig 7 — Mass and Thrust','NumberTitle','off', ...
    'Position',[350 250 fw fh])
T_plot = arrayfun(@(tt) max(interp1(p.LUT_T(:,1), p.LUT_T(:,2), ...
    min(tt, p.LUT_T(end,1)),'linear',0), 0), t_sim);
yyaxis left
plot(t_sim, (p.m_dry + mp_sim + p.m_wpay)*1000, '-','Color',c2,'LineWidth',2)
ylabel('Total mass [g]')
xlabel('Time [s]')
yyaxis right
plot(t_sim, T_plot,'--','Color',c1,'LineWidth',2)
ylabel('Thrust [N]')
title(sprintf('Total Mass and Thrust vs Time  (P_{launch} = %.1f bar)', P_launch_bar))
xline(t_burnout_val,'k:','HandleVisibility','off')
legend('Total mass','Thrust','Location','NE')
grid on
xlim([0, min(0.25, t_final)])

%% -- Figure 8: Competition Score Breakdown --
figure('Name','Fig 8 — Competition Score','NumberTitle','off', ...
    'Position',[400 200 fw fh])
dist_term = sqrt(d_final^2 + h_max^2);
b = bar([1 2 4], [dist_term, t_final, score], 'FaceColor','flat','BarWidth',0.5);
b.CData = [c1; c2; c3];
xticks([1 2 4])
xticklabels({'sqrt(d^2+h^2)  [m]','t_{flight}  [s]','SCORE'})
ylabel('Value')
title(sprintf('Competition Score Breakdown  |  Score = %.4f', score))
grid on; box off
text(1, dist_term*1.04, sprintf('%.3f m', dist_term), ...
    'HorizontalAlignment','center','FontSize',10,'FontWeight','bold')
text(2, t_final*1.04,   sprintf('%.2f s',  t_final), ...
    'HorizontalAlignment','center','FontSize',10,'FontWeight','bold')
text(4, score*1.04,     sprintf('%.4f',    score), ...
    'HorizontalAlignment','center','FontSize',11,'FontWeight','bold','Color',c3)
text(0.05, 0.92, 'score = (sqrt(d^2+h^2) + t) x (m_{pay}/m_{TOW})', ...
    'Units','normalized','FontSize',8,'BackgroundColor','w','EdgeColor','k')
text(0.05, 0.82, sprintf('     = (%.3f + %.2f) x %.4f = %.4f', ...
    dist_term, t_final, m_pay_kg/m_TOW_kg, score), ...
    'Units','normalized','FontSize',8,'BackgroundColor','w')

%% -- Figure 9: Mass Budget Table --
figure('Name','Fig 9 — Mass Budget','NumberTitle','off', ...
    'Position',[450 150 760 420])
ax_t = axes; ax_t.Visible = 'off';
n_comp   = size(mass_budget,1);
% Build row_data with plain char or plain numeric — no nested cells
row_data = cell(n_comp+1, 4);
for ii = 1:n_comp
    row_data{ii,1} = mass_budget{ii,1};   % char
    row_data{ii,2} = wet_vals(ii);         % double
    row_data{ii,3} = dry_vals(ii);         % double
    row_data{ii,4} = delta_vals(ii);       % double
end
row_data{end,1} = 'TOTAL';
row_data{end,2} = sum(wet_vals);
row_data{end,3} = sum(dry_vals);
row_data{end,4} = sum(delta_vals);
t_obj = uitable(gcf,...
    'Data',       row_data,...
    'ColumnName', {'Component','Wet [g]','Dry [g]','Propellant [g]'},...
    'Position',   [20 20 720 360],...
    'ColumnWidth', {270, 90, 90, 110},...
    'FontSize',   11,...
    'RowName',    []);
addStyle(t_obj, uistyle('BackgroundColor',[0.82 0.93 0.78],'FontWeight','bold'), ...
    'row', n_comp+1)
title(ax_t, sprintf('Mass Budget  |  m_{TOW}=%.0f g  |  m_{dry}=%.0f g  |  m_{pay}/m_{TOW}=%.4f', ...
    sum(wet_vals), sum(dry_vals), m_pay_kg/m_TOW_kg), ...
    'FontSize',11,'FontWeight','bold')

%% =========================================================================
%% SECTION 6 — EXPORT EXCEL
%% =========================================================================
fname = 'rocket_competition_results.xlsx';

T_traj = table(t_sim, xe_sim, ze_sim, V_sim, Vze_sim, ...
    rad2deg(gam_sim), rad2deg(th_sim), mp_sim*1000, ...
    'VariableNames',{'time_s','xe_m','altitude_m','V_ms','Vz_ms', ...
    'gamma_deg','theta_deg','mprop_g'});
writetable(T_traj, fname, 'Sheet','Trajectory')

T_stab = table(t_sim, xcg_sim*1000, Cp_sim*1000, SM_sim, ...
    'VariableNames',{'time_s','CG_mm_from_nose','CP_mm_from_nose', ...
    'StaticMargin_cal'});
writetable(T_stab, fname, 'Sheet','Stability')

T_mass = table(mass_budget(:,1), num2cell(wet_vals), ...
    num2cell(dry_vals), num2cell(delta_vals), ...
    'VariableNames',{'Component','Wet_g','Dry_g','Propellant_g'});
writetable(T_mass, fname, 'Sheet','MassBudget')

summary_names  = {'h_max_m';'d_final_m';'t_flight_s';'t_apogee_s';'V_max_ms'; ...
    'V_burnout_ms';'T_max_N';'t_burnout_s';'m_TOW_g';'m_dry_g'; ...
    'm_payload_g';'SM_min_cal';'SM_max_cal';'P_launch_bar';'competition_score'};
summary_values = {h_max; d_final; t_final; t_apogeo; V_max; ...
    V_at_burnout; max(p.LUT_T(:,2)); t_burnout_val; m_TOW_g; m_dry_total_g; ...
    m_payload_g; min(SM_sim); max(SM_sim); P_launch_bar; score};
T_sum = table(summary_names, summary_values, ...
    'VariableNames',{'Parameter','Value'});
writetable(T_sum, fname, 'Sheet','Summary')

fprintf('[EXPORT] Saved: %s\n', fname)
fprintf('[EXPORT]   Sheet 1: Trajectory  (%d rows)\n', height(T_traj))
fprintf('[EXPORT]   Sheet 2: Stability   (%d rows)\n', height(T_stab))
fprintf('[EXPORT]   Sheet 3: MassBudget  (%d rows)\n', height(T_mass))
fprintf('[EXPORT]   Sheet 4: Summary     (%d rows)\n', height(T_sum))

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function dydt = rocket_ode_3dof(t, y, p)
    Vxe    = y(2);
    ze     = y(3);  Vze = y(4);
    theta  = y(5);  q   = y(6);
    m_prop = max(y(7), 0);

    V     = sqrt(Vxe^2 + Vze^2);
    V_eff = max(V, 0.5);

    m_tot = p.m_dry + m_prop + p.m_wpay;
    xcg   = cg_calc(m_prop, p);

    Iyy_dry  = p.m_dry  * (3*p.R_rocket^2 + p.L_rocket^2)/12;
    Iyy_prop = m_prop   * (3*p.R_rocket^2 + p.L_tank^2  )/12;
    Iyy_pay  = p.m_wpay * (3*p.R_rocket^2 + p.L_tank^2  )/12;
    Iyy = Iyy_dry + p.m_dry *(xcg-p.xcg_dry)^2          + ...
          Iyy_prop + m_prop  *(xcg-p.xcg_tank_prop)^2    + ...
          Iyy_pay  + p.m_wpay*(xcg-p.xcg_tank_pay)^2;

    tc  = min(t, p.LUT_T(end,1));
    Cd  = interp1(p.LUT_Cd(:,1), p.LUT_Cd(:,2), tc,'linear', p.LUT_Cd(end,2));
    Cm  = interp1(p.LUT_Cm(:,1), p.LUT_Cm(:,2), tc,'linear', 0);
    Cp  = interp1(p.LUT_Cp(:,1), p.LUT_Cp(:,2), tc,'linear', p.LUT_Cp(end,2));
    Fn  = interp1(p.LUT_Fn(:,1), p.LUT_Fn(:,2), tc,'linear', 0);
    T   = max(interp1(p.LUT_T(:,1), p.LUT_T(:,2), tc,'linear', 0), 0);

    q_dyn   = 0.5 * p.rho_air * V_eff^2;
    D_aero  = q_dyn * p.A_ref * Cd;
    FN_aero = q_dyn * p.A_ref * Fn;

    if m_prop > 1e-6 && T > 1
        m_prop_dot = -T / (p.Isp_eff * p.g0);
    else
        m_prop_dot = 0;
    end

    apogeo  = double(Vze < 0 && T < 1 && ze > 0.5);
    F_chute = apogeo * 0.5 * p.rho_air * V^2 * p.Cd_chute * p.A_chute;

    if V > 0.01
        Fx_drag = -(D_aero + F_chute) * Vxe/V;
        Fz_drag = -(D_aero + F_chute) * Vze/V;
    else
        Fx_drag = 0; Fz_drag = 0;
    end

    Fx_T = T * cos(theta);       Fz_T = T * sin(theta);
    Fx_N = -FN_aero*sin(theta);  Fz_N =  FN_aero*cos(theta);

    ax = (Fx_drag + Fx_T + Fx_N) / m_tot;
    az = (Fz_drag + Fz_T + Fz_N) / m_tot - p.g0;

    arm   = xcg - Cp;
    My    = q_dyn * p.A_ref * p.D_ref * Cm + FN_aero * arm;
    q_dot = My / Iyy;

    dydt = [Vxe; ax; Vze; az; q; q_dot; m_prop_dot];
end

function xcg = cg_calc(m_prop, p)
    m_tot = p.m_dry + m_prop + p.m_wpay;
    xcg   = (p.m_dry*p.xcg_dry + m_prop*p.xcg_tank_prop + ...
             p.m_wpay*p.xcg_tank_pay) / m_tot;
end

function [value, isterminal, direction] = ground_event(t, y)
    if t > 2.0
        value = y(3);
    else
        value = 1.0;
    end
    isterminal = 1;
    direction  = -1;
end
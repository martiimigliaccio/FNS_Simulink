%% train_neural_slosh.m — R2026a
% Addestra Neural State Space per sloshing
% Richiede: slosh_training_data.mat (generato da generate_slosh_data.m)

clear; clc;
load('slosh_training_data.mat')

fprintf('=== Neural State Space Sloshing — R2026a ===\n')
fprintf('Campioni totali: %d  |  dt = %.4f s\n', length(t), dt)

%% --------------------------------------------------------
%  1. SUDDIVISIONE DATASET  70/15/15
%% --------------------------------------------------------
N       = length(t);
N_train = round(0.70 * N);
N_val   = round(0.15 * N);

idx_train = 1 : N_train;
idx_val   = N_train+1 : N_train+N_val;
idx_test  = N_train+N_val+1 : N;

fprintf('Training:   %d campioni\n', length(idx_train))
fprintf('Validation: %d campioni\n', length(idx_val))
fprintf('Test:       %d campioni\n', length(idx_test))

%% --------------------------------------------------------
%  2. OGGETTI iddata
%% --------------------------------------------------------
data_train = iddata(output_data(idx_train,:), ...
                    input_data(idx_train,:), dt, ...
                    'OutputName', {'Fz_sl','Fx_sl','My_sl'}, ...
                    'InputName',  {'Axb','Azb'});

data_val   = iddata(output_data(idx_val,:), ...
                    input_data(idx_val,:), dt, ...
                    'OutputName', {'Fz_sl','Fx_sl','My_sl'}, ...
                    'InputName',  {'Axb','Azb'});

data_test  = iddata(output_data(idx_test,:), ...
                    input_data(idx_test,:), dt, ...
                    'OutputName', {'Fz_sl','Fx_sl','My_sl'}, ...
                    'InputName',  {'Axb','Azb'});

%% --------------------------------------------------------
%  3. CREA idNeuralStateSpace
%  nx=2: stato fisico = [delta, delta_dot]
%  nu=2: input  = [Axb, Azb]
%  ny=3: output = [Fz_sl, Fx_sl, My_sl]
%  HasFeedthrough=true: Fz dipende direttamente da Azb
%% --------------------------------------------------------
nx = 2;
nu = 2;
ny = 3;

sys_nss = idNeuralStateSpace(nx, ...
    'NumInputs',     nu,   ...
    'NumOutputs',    ny,   ...
    'HasFeedthrough', true, ...
    'Ts',            dt);

fprintf('\nModello base creato:\n')
fprintf('  nx=%d  nu=%d  ny=%d  Ts=%.4f\n', nx, nu, ny, dt)

%% --------------------------------------------------------
%  4. RETE DI STATO f(x,u) -> x_next   [32,32] tanh
%% --------------------------------------------------------
f_net = createMLPNetwork(sys_nss, 'state', ...
    'LayerSizes',         [32, 32], ...
    'Activations',        'tanh',   ...
    'WeightsInitializer', 'glorot');

sys_nss = setNetwork(sys_nss, 'state', f_net);
fprintf('Rete di stato assegnata [32,32] tanh.\n')
fprintf('  Parametri stato: %d\n', sum(cellfun(@numel, f_net.Learnables.Value)))

%% --------------------------------------------------------
%  5. RETE DI OUTPUT h(x,u) -> y   [32,16] tanh
%  createMLPNetwork per output funziona in R2026a
%% --------------------------------------------------------
h_net = createMLPNetwork(sys_nss, 'output', ...
    'LayerSizes',         [32, 16], ...
    'Activations',        'tanh',   ...
    'WeightsInitializer', 'glorot');

sys_nss = setNetwork(sys_nss, 'output', h_net);
fprintf('Rete di output assegnata [32,16] tanh.\n')
fprintf('  Parametri output: %d\n', sum(cellfun(@numel, h_net.Learnables.Value)))

%% Verifica che le reti abbiano parametri apprendibili


sys_nss.InputName  = {'Axb', 'Azb'};
sys_nss.OutputName = {'Fz_sl', 'Fx_sl', 'My_sl'};

%% --------------------------------------------------------
%  6. OPZIONI TRAINING  (Adam con learning rate schedule)
%% --------------------------------------------------------
opt = nssTrainingOptions('adam');
opt.MaxEpochs                  = 1000;
opt.LearnRate                  = 0.005;
opt.GradientDecayFactor        = 0.9;
opt.SquaredGradientDecayFactor = 0.999;
opt.WindowSize = 2048;
opt.LossFcn                    = 'MeanAbsoluteError';
opt.PlotLossFcn                = 1;
opt.LearnRateSchedule          = 'piecewise';
opt.LearnRateDropFactor        = 0.5;
opt.LearnRateDropPeriod        = 200;

fprintf('\nOpzioni training:\n')
fprintf('  MaxEpochs=%d | LearnRate=%.4f | MiniBatch=%d\n', ...
    opt.MaxEpochs, opt.LearnRate, opt.WindowSize)
fprintf('  LR drop: x%.1f ogni %d epoche\n', ...
    opt.LearnRateDropFactor, opt.LearnRateDropPeriod)

%% --------------------------------------------------------
%  7. TRAINING
%% --------------------------------------------------------
fprintf('\nAvvio training...\n')
sys_nss_trained = nlssest(data_train, sys_nss, opt);

%% --------------------------------------------------------
%  8. VALIDAZIONE
%% --------------------------------------------------------
fprintf('\n=== Risultati ===\n')
[~, fit_val]  = compare(sys_nss_trained, data_val);
[~, fit_test] = compare(sys_nss_trained, data_test);

fprintf('Fit VALIDATION:  Fz=%.1f%%  Fx=%.1f%%  My=%.1f%%\n', ...
    fit_val(1),  fit_val(2),  fit_val(3))
fprintf('Fit TEST:        Fz=%.1f%%  Fx=%.1f%%  My=%.1f%%\n', ...
    fit_test(1), fit_test(2), fit_test(3))

%% --------------------------------------------------------
%  9. PLOT CONFRONTO TEST SET
%% --------------------------------------------------------
t_test  = t(idx_test);
y_pred  = sim(sys_nss_trained, data_test);
idx_plt = 1:min(2000, length(t_test));

figure('Name','Neural State Space — Validazione finale')
titles = {'Fz_{sl} [N]', 'Fx_{sl} [N]', 'My_{sl} [N·m]'};
for k = 1:3
    subplot(3,1,k)
    plot(t_test(idx_plt), output_data(idx_test(idx_plt),k), ...
         'b',  'LineWidth', 1.5, 'DisplayName', 'Pendolo riferimento')
    hold on
    plot(t_test(idx_plt), y_pred.OutputData(idx_plt,k), ...
         'r--','LineWidth', 1.2, 'DisplayName', 'Neural SS')
    legend('Location','best')
    title(titles{k})
    grid on
end
xlabel('Tempo [s]')

%% --------------------------------------------------------
%  10. SALVA
%% --------------------------------------------------------
save('neural_slosh_model.mat', 'sys_nss_trained', ...
     'fit_val', 'fit_test', 'nx', 'dt', ...
     'mu_in', 'std_in', 'mu_out', 'std_out')

fprintf('\nModello salvato in neural_slosh_model.mat\n')
fprintf('Prossimo step: integrazione in Simulink.\n')
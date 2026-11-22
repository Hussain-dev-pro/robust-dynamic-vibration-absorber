%% DVA_Simulink_Build.m  -  builds the Simulink model from code
% Run DVA_Main.m first (it creates DVA_results.mat).
% Creates DVA_Simulink.slx: one unbalance force drives (a) the bare machine
% and (b) the machine + final robust absorber, both shown on one Scope.
clear; clc;
S = load('DVA_results.mat');                 % R, m1, k1, c1, m2, F0, w_op
m1 = S.m1; k1 = S.k1; c1 = S.c1; m2 = S.m2;
k2 = S.R.k2; c2 = S.R.c2; F0 = S.F0; w_op = S.w_op;

% State vector [x1 v1 x2 v2]  (x = displacement, v = velocity)
A_dva = [ 0             1              0       0;
         -(k1+k2)/m1  -(c1+c2)/m1     k2/m1   c2/m1;
          0             0              0       1;
          k2/m2         c2/m2         -k2/m2  -c2/m2 ];
B_dva = [0; 1/m1; 0; 0];
C_dva = [1 0 0 0];      D_dva = 0;           % output: machine x1
A_0   = [0 1; -k1/m1 -c1/m1];  B_0 = [0; 1/m1];  C_0 = [1 0];  D_0 = 0;

mdl = 'DVA_Simulink';
if bdIsLoaded(mdl), close_system(mdl, 0); end
if exist([mdl '.slx'], 'file'), delete([mdl '.slx']); end
new_system(mdl); open_system(mdl);

add_block('simulink/Sources/Sine Wave', [mdl '/Unbalance force F0 sin(w t)'], ...
    'Amplitude','F0', 'Frequency','w_op', 'Position',[40 120 90 160]);
add_block('simulink/Continuous/State-Space', [mdl '/Machine only'], ...
    'A','A_0','B','B_0','C','C_0','D','D_0', 'Position',[200 40 300 90]);
add_block('simulink/Continuous/State-Space', [mdl '/Machine + DVA'], ...
    'A','A_dva','B','B_dva','C','C_dva','D','D_dva', 'Position',[200 190 300 240]);
add_block('simulink/Signal Routing/Mux', [mdl '/Mux'], 'Inputs','2', 'Position',[380 110 385 170]);
add_block('simulink/Sinks/Scope', [mdl '/x1 compare'], 'Position',[450 120 490 160]);
add_block('simulink/Sinks/To Workspace', [mdl '/x1_none'], ...
    'VariableName','x1_none', 'SaveFormat','Timeseries', 'Position',[380 30 450 60]);
add_block('simulink/Sinks/To Workspace', [mdl '/x1_dva'], ...
    'VariableName','x1_dva',  'SaveFormat','Timeseries', 'Position',[380 220 450 250]);

src = 'Unbalance force F0 sin(w t)/1';
add_line(mdl, src, 'Machine only/1',  'autorouting','on');
add_line(mdl, src, 'Machine + DVA/1', 'autorouting','on');
add_line(mdl, 'Machine only/1',  'Mux/1',     'autorouting','on');
add_line(mdl, 'Machine + DVA/1', 'Mux/2',     'autorouting','on');
add_line(mdl, 'Mux/1',           'x1 compare/1','autorouting','on');
add_line(mdl, 'Machine only/1',  'x1_none/1', 'autorouting','on');
add_line(mdl, 'Machine + DVA/1', 'x1_dva/1',  'autorouting','on');

set_param(mdl, 'StopTime','5', 'Solver','ode45', 'MaxStep','1e-3', 'RelTol','1e-7');
save_system(mdl);

out = sim(mdl);
a = out.get('x1_none');  b = out.get('x1_dva');
figure('Color','w');
plot(a.Time, a.Data*1e3, 'k', b.Time, b.Data*1e3, 'r'); grid on;
xlabel('Time (s)'); ylabel('x_1 (mm)'); legend('No absorber','Robust DVA');
title('Simulink result at 955 rpm');
ss = a.Time > 4;
fprintf('Simulink steady amplitude: no DVA %.3f mm, with DVA %.3f mm\n', ...
    max(abs(a.Data(ss)))*1e3, max(abs(b.Data(ss)))*1e3);
fprintf('Compare with DVA_Main.m frequency-domain result -> they should agree.\n');

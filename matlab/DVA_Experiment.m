%% DVA_Experiment.m  -  table-top validation with a phone accelerometer
% PART A: size a small cantilever rig + absorber (same physics, scaled down)
% PART B: read phyphox CSV, FFT, find peaks, estimate damping, compare
%
% Rig: aluminium flat bar clamped to a table edge (primary system), phone
% taped at the tip (its mass is part of m1, its accelerometer is the sensor).
% Absorber: thin steel strip clamped at the tip with nuts/coins as m2.
% Test: tap the tip (impact test), record with phyphox
% "Acceleration (without g)", export CSV. Do it WITHOUT and WITH absorber.
clear; clc; close all;

%% PART A - rig design ------------------------------------------------------
E_al = 69e9;  rho_al = 2700;  b1 = 25e-3;  h1 = 3e-3;  L1 = 0.30;   % bar
m_tip = 0.20;                        % kg: phone (~0.19) + tape  -> WEIGH IT
I1  = b1*h1^3/12;   k1 = 3*E_al*I1/L1^3;
m1  = m_tip + 0.236*rho_al*b1*h1*L1; % Rayleigh effective tip mass
f1  = sqrt(k1/m1)/(2*pi);
mu  = 0.15;  m2 = mu*m1;             % absorber mass (nuts / coins)
f2  = f1/(1+mu);                     % Den Hartog tuning target
E_st = 200e9; rho_st = 7850; b2 = 12.7e-3; h2 = 0.8e-3;            % steel strip
I2  = b2*h2^3/12;
L2  = 0.05;                          % iterate for strip length incl. its mass
for it = 1:50
    meff = m2 + 0.236*rho_st*b2*h2*L2;
    L2   = (3*E_st*I2/(meff*(2*pi*f2)^2))^(1/3);
end
fprintf('PRIMARY : k1 = %.0f N/m, m1 = %.0f g, predicted f1 = %.2f Hz\n', k1, m1*1e3, f1);
fprintf('ABSORBER: m2 = %.0f g, target f2 = %.2f Hz -> free strip length L2 = %.0f mm\n', m2*1e3, f2, L2*1e3);
fprintf('Tip: make L2 a bit long, then slide the clamp to tune (measure f2 alone first).\n\n');

%% PART B - phyphox data ---------------------------------------------------
files = {'no_absorber.csv', 'with_absorber.csv'};   % your exports
axisCol = 4;          % phyphox columns: 1 time, 2 x, 3 y, 4 z (phone flat -> z)
figure('Color','w'); hold on; col = {'k','r'};
for j = 1:2
    if exist(files{j}, 'file')
        M = readmatrix(files{j});  t = M(:,1);  a = M(:,axisCol);
        src = 'measured';
    else                                   % demo signal so the script runs
        [t,a] = demoSignal(j, f1, f2, mu);  src = 'DEMO (no csv found)';
    end
    ok = ~isnan(t) & ~isnan(a);  t = t(ok);  a = a(ok);
    fs = 1/median(diff(t));                % phone sampling is not uniform:
    tu = (t(1):1/fs:t(end))';              % resample to a uniform grid
    au = interp1(t, a, tu, 'linear');
    au = au - mean(au);
    N  = numel(au);  win = 0.5 - 0.5*cos(2*pi*(0:N-1)'/(N-1));   % Hann
    A  = abs(fft(au.*win))/N;  A = A(1:floor(N/2));  f = (0:numel(A)-1)'*fs/N;
    sel = f > 1 & f < 4*f1;
    plot(f(sel), A(sel)/max(A(sel)), col{j}, 'LineWidth',1.2);
    pk = localPeaks(f(sel), A(sel), 2);
    fprintf('%s (%s): fs = %.0f Hz, peaks at %s Hz\n', files{j}, src, fs, mat2str(round(pk,2)));
    for p = pk(:)'
        fprintf('   peak %.2f Hz: damping ratio ~ %.3f (half-power)\n', p, halfPower(f(sel),A(sel),p));
    end
end
xlabel('Frequency (Hz)'); ylabel('Normalised |A|'); grid on;
legend('Without absorber','With absorber'); title('Measured spectra (tap test)');
fprintf('\nExpected: one peak near %.2f Hz without absorber; two peaks around it with absorber.\n', f1);
fprintf('Report: measured vs predicted f1 (%% error) and the split peaks vs the 2-DOF model.\n');

%% local functions
function pk = localPeaks(f, A, n)
% up to n peaks, each > 20% of the largest and at least 0.5 Hz apart
    idx = find(A(2:end-1) > A(1:end-2) & A(2:end-1) > A(3:end)) + 1;
    idx = idx(A(idx) > 0.2*max(A));
    [~, s] = sort(A(idx), 'descend');  idx = idx(s);  pk = [];
    for i = idx(:)'
        if all(abs(f(i) - pk) > 0.5), pk(end+1) = f(i); end %#ok<AGROW>
        if numel(pk) == n, break; end
    end
    pk = sort(pk);
end
function z = halfPower(f, A, fp)
    [~, i] = min(abs(f - fp));  h = A(i)/sqrt(2);
    il = i; while il > 1 && A(il) > h, il = il - 1; end
    ir = i; while ir < numel(A) && A(ir) > h, ir = ir + 1; end
    z = (f(ir) - f(il))/(2*fp);
end
function [t,a] = demoSignal(j, f1, f2, mu)
    fs = 400; t = (0:1/fs:20)'; a = 0.02*randn(size(t));
    if j == 1, fr = f1; else                     % split peaks of 2-DOF (undamped estimate)
        r = f2/f1; s = 1 + (1+mu)*r^2; q = sqrt(s^2 - 4*r^2);
        fr = f1*sqrt([(s-q)/2, (s+q)/2]);
    end
    for tap = 1                                  % ONE tap per record
        tt = t - tap; on = tt > 0;
        for fk = fr, a(on) = a(on) + exp(-2*pi*fk*0.02*tt(on)).*sin(2*pi*fk*tt(on)); end
    end
end

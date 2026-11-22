%% ========================================================================
%  DVA_Main.m  -  Robust Dynamic Vibration Absorber for a Pump at Resonance
%  ------------------------------------------------------------------------
%  Case study: a 200 kg motor-pump on a steel frame runs at 955 rpm.
%  The frame's vertical natural frequency is also 15.9 Hz (955 rpm), so the
%  machine sits on resonance: 27.8 mm/s rms -> ISO 10816-3 Zone D.
%  Goal: design a passive tuned mass damper (DVA) that brings vibration to
%  Zone B (< 4.5 mm/s) over the full +/-10% speed band, even when the frame
%  stiffness is only known to +/-15%.
%
%  Sections
%   1  Inputs (force derived from ISO 21940-11 balance grade G6.3)
%   2  Baseline problem
%   3  Mass-ratio selection (design criterion: ISO limit)
%   4  Four absorber designs: undamped, Den Hartog, minimax, ROBUST
%   5  Frequency response + results table
%   6  Robustness / tolerance study
%   7  Start-up simulation (ode45, unbalance force grows with speed^2)
%   8  Hardware sizing: mass block, helical springs, elastomer damper
%   9  Save results
%
%  Needs only base MATLAB (fminsearch, ode45). Tested logic mirrored in
%  Python; numbers in the report come from the same equations.
%  ========================================================================
clear; clc; close all;
figdir = 'figures'; if ~isfolder(figdir), mkdir(figdir); end

%% 1) INPUTS ---------------------------------------------------------------
m1    = 200;                       % kg   machine + effective frame mass
k1    = 2e6;                       % N/m  frame vertical stiffness
zeta1 = 0.02;                      % -    welded steel frame damping
c1    = 2*zeta1*sqrt(k1*m1);       % N.s/m (= 800)

n_op  = 955;                       % rpm  operating speed
w_op  = n_op*2*pi/60;              % rad/s

% Unbalance force from ISO 21940-11: e_per = G / omega
m_rot = 50;                        % kg   rotor mass
Gq    = 6.3;                       % mm/s balance grade G6.3 (pumps)
e_per = (Gq*1e-3)/w_op;            % m    permissible eccentricity
U     = m_rot*e_per;               % kg.m residual unbalance
F0    = U*w_op^2;                  % N    force amplitude at w_op (~31.5 N)

wn    = sqrt(k1/m1);               % rad/s
w     = linspace(0.5*wn, 1.5*wn, 3001);
band  = (w >= 0.9*w_op) & (w <= 1.1*w_op);    % operating speed band
kScale = linspace(0.85, 1.15, 13);            % frame stiffness uncertainty

ISO = [2.3 4.5 7.1];   % ISO 10816-3 group 2, flexible support: A/B, B/C, C/D (mm/s rms)
vrms = @(X,ww) X.*ww/sqrt(2)*1e3;              % displacement amp (m) -> mm/s rms

%% 2) BASELINE PROBLEM -----------------------------------------------------
X0 = abs(F0 ./ (k1 - m1*w.^2 + 1i*w*c1));
[~,iop] = min(abs(w - w_op));
fprintf('==== BASELINE ====\n');
fprintf('Unbalance U = %.0f g.mm  ->  F0 = %.1f N\n', U*1e6, F0);
fprintf('Natural freq = %.2f Hz (%.0f rpm), operating = %.2f Hz\n', wn/2/pi, wn*60/2/pi, w_op/2/pi);
fprintf('Amplitude at op = %.3f mm,  v_rms = %.1f mm/s  -> Zone %s\n', ...
        X0(iop)*1e3, vrms(X0(iop),w_op), isoZone(vrms(X0(iop),w_op),ISO));
fprintf('Force into foundation = %.0f N\n\n', abs(k1+1i*w_op*c1)*X0(iop));

%% 3) MASS-RATIO SELECTION -------------------------------------------------
mus = 0.02:0.01:0.30;  vBand = zeros(size(mus));
for j = 1:numel(mus)
    [k2j,c2j] = denHartog(mus(j)*m1, mus(j), wn);
    X = dva_frf(m1,k1,c1,mus(j)*m1,k2j,c2j,F0,w);
    vBand(j) = max(vrms(X(band),w(band)));
end
muMin = mus(find(vBand <= ISO(2), 1));
fprintf('==== MASS RATIO ====\n');
fprintf('Smallest mu meeting %.1f mm/s over the speed band: %.2f\n', ISO(2), muMin);
mu = 0.15;  m2 = mu*m1;        % chosen with margin (mu = 0.10 FAILS: 4.73 mm/s)
fprintf('Chosen mu = %.2f  ->  m2 = %.0f kg\n\n', mu, m2);

figure('Color','w'); plot(mus, vBand, 'o-', 'LineWidth',1.4); hold on;
yline(ISO(2),'r--','Zone B/C 4.5 mm/s'); yline(ISO(1),'g:','Zone A/B 2.3 mm/s');
xlabel('Mass ratio \mu = m_2/m_1'); ylabel('Worst v_{rms} in speed band (mm/s)');
title('Mass-ratio trade-off (Den Hartog tuning)'); grid on;
saveas(gcf, fullfile(figdir,'fig3_mass_ratio.png'));

%% 4) FOUR ABSORBER DESIGNS ------------------------------------------------
D = struct();
% (a) Undamped (classic)
D(1).name = 'Undamped';     D(1).k2 = m2*wn^2;  D(1).c2 = 0;
% (b) Den Hartog closed form
[k2dh, c2dh] = denHartog(m2, mu, wn);
D(2).name = 'Den Hartog';   D(2).k2 = k2dh;      D(2).c2 = c2dh;
% (c) Numerical minimax: minimise global peak of X1 (covers start-up too)
opts = optimset('TolX',1e-8,'TolFun',1e-12,'MaxIter',1e4,'MaxFunEvals',1e4,'Display','off');
Jg   = @(p) max(dva_frf(m1,k1,c1,m2,exp(p(1)),exp(p(2)),F0,w));
p    = fminsearch(Jg, log([D(2).k2 D(2).c2]), opts);
D(3).name = 'Minimax';      D(3).k2 = exp(p(1)); D(3).c2 = exp(p(2));
% (d) ROBUST: minimise worst v_rms in the speed band over k1 +/-15%
Jr   = @(p) worstBandVel(m1,k1*kScale,c1,m2,exp(p(1)),exp(p(2)),F0,w,band);
p    = fminsearch(Jr, log([D(2).k2 D(2).c2]), opts);
D(4).name = 'Robust (final)'; D(4).k2 = exp(p(1)); D(4).c2 = exp(p(2));

%% 5) FREQUENCY RESPONSE + RESULTS TABLE -----------------------------------
figure('Color','w'); hold on;
fill([0.9 1.1 1.1 0.9]*w_op/2/pi, [0 0 0.5 0.5], [1 .85 .6], 'EdgeColor','none','FaceAlpha',.4);
plot(w/2/pi, X0*1e3, 'k', 'LineWidth',1.6);
sty = {'b--','g-.','m:','r-'};
fprintf('==== DESIGNS (m2 = %.0f kg) ====\n', m2);
fprintf('%-15s %9s %8s %6s %6s %9s %9s %10s %8s\n', 'Design','k2 N/m','c2 Ns/m', ...
        'zeta2','f','Peak mm','v_op','v_band','Worst k1');
for d = 1:4
    [X1,Xrel] = dva_frf(m1,k1,c1,m2,D(d).k2,D(d).c2,F0,w);
    D(d).zeta = D(d).c2/(2*sqrt(D(d).k2*m2));
    D(d).f    = sqrt(D(d).k2/m2)/wn;
    D(d).peak = max(X1)*1e3;
    D(d).vop  = vrms(X1(iop),w_op);
    D(d).vband= max(vrms(X1(band),w(band)));
    D(d).vwc  = worstBandVel(m1,k1*kScale,c1,m2,D(d).k2,D(d).c2,F0,w,band);
    D(d).relMax = max(Xrel)*1e3;
    D(d).FT   = abs(k1+1i*w_op*c1)*X1(iop);
    plot(w/2/pi, X1*1e3, sty{d}, 'LineWidth',1.5);
    fprintf('%-15s %9.0f %8.0f %6.3f %6.3f %9.3f %8.2f %10.2f %8.2f\n', D(d).name, ...
        D(d).k2, D(d).c2, D(d).zeta, D(d).f, D(d).peak, D(d).vop, D(d).vband, D(d).vwc);
end
xlabel('Frequency (Hz)'); ylabel('Machine amplitude X_1 (mm)'); grid on; box on;
legend('Speed band \pm10%','No absorber', D.name, 'Location','northeast');
title('Machine response with and without absorber'); ylim([0 0.5]);
saveas(gcf, fullfile(figdir,'fig1_frf.png'));
R = D(4);   % final design
fprintf('\nFINAL: v_op %.1f -> %.1f mm/s (%.0f%% lower), foundation force %.0f -> %.0f N\n\n', ...
    vrms(X0(iop),w_op), R.vop, (1-R.vop/vrms(X0(iop),w_op))*100, abs(k1+1i*w_op*c1)*X0(iop), R.FT);

%% 6) ROBUSTNESS / TOLERANCE ----------------------------------------------
fprintf('==== TOLERANCE (final design, worst v_band over k1 +/-15%%) ====\n');
for s = [0.9 0.95 1.05 1.1]
    fprintf('k2 x %.2f -> %.2f mm/s\n', s, worstBandVel(m1,k1*kScale,c1,m2,R.k2*s,R.c2,F0,w,band));
end
for s = [0.5 0.75 1.5 2]
    fprintf('c2 x %.2f -> %.2f mm/s\n', s, worstBandVel(m1,k1*kScale,c1,m2,R.k2,R.c2*s,F0,w,band));
end
fprintf('\n');

%% 7) START-UP SIMULATION (ode45) -----------------------------------------
T = 8; tEnd = 12;                          % 0 -> 955 rpm in 8 s
odeo = odeset('RelTol',1e-7,'AbsTol',1e-11,'MaxStep',1e-3);
cases = {'No absorber',0,0; 'Undamped',D(1).k2,0; 'Robust (final)',R.k2,R.c2};
figure('Color','w'); hold on; col = {'k','b','r'};
fprintf('==== START-UP (0 -> %d rpm in %d s) ====\n', n_op, T);
for j = 1:3
    k2 = cases{j,2}; c2 = cases{j,3};
    [t,y] = ode45(@(t,y) runupRHS(t,y,m1,k1,c1,m2,k2,c2,U,w_op,T), [0 tEnd], zeros(4,1), odeo);
    if k2 == 0, rel = 0; else, rel = max(abs(y(:,1)-y(:,3)))*1e3; end
    fprintf('%-15s max x1 = %.3f mm, steady x1 = %.3f mm, max stroke = %.3f mm\n', ...
        cases{j,1}, max(abs(y(:,1)))*1e3, max(abs(y(t>tEnd-2,1)))*1e3, rel);
    plot(t, y(:,1)*1e3, col{j}, 'LineWidth',0.6);
end
xlabel('Time (s)'); ylabel('x_1 (mm)'); grid on; legend(cases(:,1));
title('Start-up: unbalance force grows with speed^2');
saveas(gcf, fullfile(figdir,'fig2_startup.png'));
fprintf('\n');

%% 8) HARDWARE SIZING ------------------------------------------------------
fprintf('==== HARDWARE (final design) ====\n');
rho = 7850;  Vb = m2/rho;                      % steel block volume
fprintf('Mass block: %.1f L steel, e.g. 200 x 200 x %.0f mm (+ 2 x 1 kg tuning plates)\n', Vb*1e3, Vb/0.04*1e3);

% 4 helical compression springs in parallel (EN 10270-1 / ASTM A229 oil-tempered)
G  = 79.3e9;  d = 8e-3;  Dm = 48e-3;  C = Dm/d;   nSpr = 4;
ke = R.k2/nSpr;
na = G*d^4/(8*Dm^3*ke);                        % active coils
nt = na + 2;  Ls = nt*d;                       % squared & ground ends
Fst = m2*9.81/nSpr;  yst = Fst/ke;
Fa  = ke*R.relMax*1e-3;                        % alternating force at max stroke
Ks  = 1 + 0.5/C;                               % direct shear factor (mean stress)
Kw  = (4*C-1)/(4*C-4) + 0.615/C;               % Wahl factor (alternating)
tm  = Ks*8*Fst*Dm/(pi*d^3)/1e6;  ta = Kw*8*Fa*Dm/(pi*d^3)/1e6;   % MPa
Sut = 1610/(d*1e3)^0.193;  Ssu = 0.67*Sut;     % A229 (Shigley Table 10-4)
Sse = 241/(1 - 379/Ssu);                       % Zimmerli unpeened -> Goodman
nF  = 1/(ta/Sse + tm/Ssu);
L0  = Ls + yst + R.relMax*1e-3 + 0.15*yst + 0.010;   % solid + travel + clash allowance
mS  = rho*pi*d^2/4*pi*Dm*na;  fSurge = 0.5*sqrt(ke/mS);
fprintf('Springs: 4 x (d = %.0f mm, D = %.0f mm, C = %.1f, %.1f active / %.1f total coils)\n', d*1e3, Dm*1e3, C, na, nt);
fprintf('  each k = %.1f kN/m, static deflection = %.2f mm, free length ~ %.0f mm, L0/D = %.2f (< 5.26: no buckling)\n', ke/1e3, yst*1e3, L0*1e3, L0/Dm);
fprintf('  tau_mean = %.1f MPa, tau_alt = %.2f MPa, Goodman fatigue FoS = %.0f\n', tm, ta, nF);
fprintf('  surge frequency = %.0f Hz (>> run speed %.1f Hz, 2x harmonic %.1f Hz): OK\n', fSurge, w_op/2/pi, w_op/pi);
eta = 2*R.zeta;
fprintf('Damper: c2 = %.0f N.s/m  ->  elastomer loss factor eta ~ %.2f at %.1f Hz\n', R.c2, eta, R.f*wn/2/pi);
fprintf('  (e.g. neoprene/butyl pads under the springs; confirm eta from datasheet at temperature)\n');
fprintf('  peak damper force ~ %.1f N, absorber stroke %.2f mm\n\n', R.c2*R.f*wn*R.relMax*1e-3, R.relMax);

%% 9) SAVE -----------------------------------------------------------------
save('DVA_results.mat','D','R','m1','k1','c1','m2','F0','U','w_op');
fprintf('Saved DVA_results.mat and figures/*.png\n');

%% ======================= LOCAL FUNCTIONS ================================
function [X1,Xrel] = dva_frf(m1,k1,c1,m2,k2,c2,F0,w)
% 2-DOF steady-state amplitudes for force F0 on m1 (complex impedance)
    Z11 = k1 + k2 - m1*w.^2 + 1i*w*(c1+c2);
    Z12 = -(k2 + 1i*w*c2);
    Z22 = k2 - m2*w.^2 + 1i*w*c2;
    Dt  = Z11.*Z22 - Z12.^2;
    X1  = abs(F0*Z22./Dt);
    Xrel= abs(F0*(-Z12 - Z22)./Dt);     % |X2 - X1| = absorber stroke
end

function [k2,c2] = denHartog(m2,mu,wn)
    f    = 1/(1+mu);
    zeta = sqrt(3*mu/(8*(1+mu)^3));
    wa   = f*wn;
    k2   = m2*wa^2;  c2 = 2*zeta*m2*wa;
end

function v = worstBandVel(m1,k1set,c1,m2,k2,c2,F0,w,band)
% worst v_rms (mm/s) inside the speed band over all frame stiffness cases
    v = 0;
    for k1 = k1set
        X = dva_frf(m1,k1,c1,m2,k2,c2,F0,w);
        v = max(v, max(X(band).*w(band))/sqrt(2)*1e3);
    end
end

function dy = runupRHS(t,y,m1,k1,c1,m2,k2,c2,U,w_op,T)
% linear speed ramp 0 -> w_op in T s; unbalance force U*w(t)^2*sin(theta)
    if t < T, wt = w_op*t/T; th = w_op*t^2/(2*T);
    else,     wt = w_op;     th = w_op*T/2 + w_op*(t-T); end
    F  = U*wt^2*sin(th);
    a1 = (F - c1*y(2) - k1*y(1) - c2*(y(2)-y(4)) - k2*(y(1)-y(3)))/m1;
    if k2 == 0, a2 = 0; else, a2 = (-c2*(y(4)-y(2)) - k2*(y(3)-y(1)))/m2; end
    dy = [y(2); a1; y(4); a2];
end

function z = isoZone(v,ISO)
    zones = 'ABCD';  z = zones(find([v <= ISO, true], 1));
end

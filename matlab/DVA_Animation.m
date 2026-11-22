%% DVA_Animation.m  -  animated proof: pump WITHOUT vs WITH the robust DVA
% Run DVA_Main.m first (creates DVA_results.mat).
% Scene 1 (6 s): both machines running at 955 rpm, side by side.
% Scene 2 (10 s): speed sweep 480 -> 1440 rpm; the cursor runs along the
%                 frequency-response curves while the machines move.
% Motion is magnified x150 and slowed so the eye can follow it. Amplitudes
% and phases are the exact steady-state 2-DOF results at each speed.
% Outputs (folder "animation"): DVA_animation.mp4 (or .avi), DVA_animation.gif,
% frame_955rpm.png, frame_sweep_resonance.png
clear; clc; close all;
S = load('DVA_results.mat');
m1 = S.m1; k1 = S.k1; c1 = S.c1; m2 = S.m2; F0 = S.F0; w_op = S.w_op;
k2 = S.R.k2; c2 = S.R.c2;
mag  = 150;          % motion magnification
fps  = 20;           % video frame rate
fvis = 1.5;          % visual oscillation rate (Hz) - slow motion
outdir = 'animation'; if ~isfolder(outdir), mkdir(outdir); end
ISO = [2.3 4.5 7.1];

% ---------- frequency-response curves for the lower panel ----------------
fH = linspace(8, 24, 800);
[A0,A1] = deal(zeros(size(fH)));
for j = 1:numel(fH)
    [X0,X1] = cplx(m1,k1,c1,m2,k2,c2,F0,2*pi*fH(j));
    A0(j) = abs(X0)*1e3;  A1(j) = abs(X1)*1e3;
end

% ---------- figure and static artwork ------------------------------------
fig = figure('Color','w','Position',[40 40 1280 720],'Name','DVA animation');
try, theme(fig,'light'); end %#ok<TRYNC>
k = [0 0 0];  steel = [0.27 0.47 0.72];  dva = [0.93 0.55 0.15];

axL = axes(fig,'Position',[0.03 0.34 0.30 0.60]);
axR = axes(fig,'Position',[0.35 0.34 0.30 0.60]);
for ax = [axL axR]
    hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');
    xlim(ax,[-210 210]); ylim(ax,[-30 640]); set(ax,'Color','w');
    plot(ax,[-170 170],[0 0],'Color',k,'LineWidth',2);
    for xh = -165:20:155, plot(ax,[xh xh+14],[0 -14],'Color',k); end
end
title(axL,'WITHOUT absorber','Color',[0.75 0.1 0.1],'FontSize',15);
title(axR,'WITH robust DVA (30 kg)','Color',[0.1 0.5 0.15],'FontSize',15);

% moving parts: springs, dampers, blocks, unbalance discs
G.L = mkParts(axL, steel, k, false, dva);
G.R = mkParts(axR, steel, k, true,  dva);
G.L.txt = text(axL,0,595,'','HorizontalAlignment','center','FontSize',12,'FontWeight','bold');
G.R.txt = text(axR,0,595,'','HorizontalAlignment','center','FontSize',12,'FontWeight','bold');

axF = axes(fig,'Position',[0.06 0.07 0.58 0.22]); hold(axF,'on'); box(axF,'on'); grid(axF,'on');
set(axF,'Color','w','XColor',k,'YColor',k,'GridColor',[0.6 0.6 0.6]);
fill(axF,[0.9 1.1 1.1 0.9]*w_op/2/pi,[0 0 0.45 0.45],[1 0.88 0.65],'EdgeColor','none');
plot(axF,fH,A0,'k','LineWidth',1.6); plot(axF,fH,A1,'r','LineWidth',2);
xlim(axF,[8 24]); ylim(axF,[0 0.45]);
xlabel(axF,'Frequency (Hz)','Color',k); ylabel(axF,'X_1 (mm)','Color',k);
cur = xline(axF,w_op/2/pi,'b-','LineWidth',1.5);
mk0 = plot(axF,NaN,NaN,'ko','MarkerFaceColor','k','MarkerSize',7);
mk1 = plot(axF,NaN,NaN,'ro','MarkerFaceColor','r','MarkerSize',7);
legend(axF,{'Speed band \pm10%','No absorber','Robust DVA'},'Location','northeast','TextColor',k,'Color','w');

axI = axes(fig,'Position',[0.67 0.07 0.31 0.87]); axis(axI,'off'); xlim(axI,[0 1]); ylim(axI,[0 1]);
text(axI,0,0.97,'Robust Dynamic Vibration Absorber','FontSize',15,'FontWeight','bold','Color',k);
info = {'Pump + frame: m_1 = 200 kg, k_1 = 2.0 MN/m, \zeta_1 = 2%', ...
        sprintf('Unbalance: G6.3, F_0 = %.1f N at 955 rpm', F0), ...
        sprintf('Absorber: m_2 = %.0f kg, k_2 = %.1f kN/m, c_2 = %.0f N\\cdot{}s/m', m2, k2/1e3, c2), ...
        'Limit: ISO 10816-3 Zone B/C = 4.5 mm/s'};
for j = 1:numel(info), text(axI,0,0.89-0.055*j,info{j},'FontSize',11,'Color',k); end
hScene = text(axI,0,0.58,'','FontSize',13,'FontWeight','bold','Color',[0.1 0.2 0.6]);
hFreq  = text(axI,0,0.51,'','FontSize',13,'Color',k);
text(axI,0,0.38,{'Motion magnified \times150 and slowed down.', ...
     'Amplitudes and phases are the exact', 'steady-state 2-DOF results at each speed.'}, ...
     'FontSize',10,'Color',[0.35 0.35 0.35],'VerticalAlignment','top');
text(axI,0,0.12,{'Generated in MATLAB by DVA_Animation.m', ['Run: ' char(datetime('now','Format','dd-MMM-yyyy HH:mm'))]}, ...
     'FontSize',9,'Color',[0.35 0.35 0.35],'Interpreter','none','VerticalAlignment','top');

% ---------- writers -------------------------------------------------------
try
    vw = VideoWriter(fullfile(outdir,'DVA_animation.mp4'),'MPEG-4');
catch
    vw = VideoWriter(fullfile(outdir,'DVA_animation.avi'),'Motion JPEG AVI');
end
vw.FrameRate = fps; vw.Quality = 95; open(vw);
gifFile = fullfile(outdir,'DVA_animation.gif');

N1 = 6*fps;  N2 = 10*fps;  N = N1 + N2;
fSweep = linspace(8, 24, N2);
[~, iPk] = min(abs(fSweep - w_op/2/pi));  iPk = iPk + N1;   % frame at resonance in sweep
sz = [];
for n = 1:N
    t  = (n-1)/fps;  th = 2*pi*fvis*t;
    if n <= N1
        f = w_op/2/pi;  set(hScene,'String','Scene 1: running at 955 rpm');
    else
        f = fSweep(n-N1); set(hScene,'String','Scene 2: speed sweep 480 \rightarrow 1440 rpm');
    end
    w = 2*pi*f;
    [X0,X1,X2] = cplx(m1,k1,c1,m2,k2,c2,F0,w);
    xL  = real(X0*exp(1i*th));  x1 = real(X1*exp(1i*th));  x2 = real(X2*exp(1i*th));
    upParts(G.L, mag*xL*1e3, 0, th, false);
    upParts(G.R, mag*x1*1e3, mag*x2*1e3, th, true);
    setReadout(G.L.txt, abs(X0), w, ISO);
    setReadout(G.R.txt, abs(X1), w, ISO);
    set(hFreq,'String',sprintf('Speed: %4.0f rpm  (%.1f Hz)', f*60, f));
    set(cur,'Value',f);
    set(mk0,'XData',f,'YData',abs(X0)*1e3);  set(mk1,'XData',f,'YData',abs(X1)*1e3);
    drawnow;
    fr = getframe(fig);  img = fr.cdata;
    if isempty(sz), sz = size(img); else, img = fitTo(img, sz); end
    writeVideo(vw, img);
    if mod(n,2) == 1                                   % GIF at 10 fps, half size
        [A,map] = rgb2ind(img(1:2:end,1:2:end,:), 128, 'nodither');
        if n == 1, imwrite(A,map,gifFile,'gif','LoopCount',Inf,'DelayTime',2/fps);
        else,      imwrite(A,map,gifFile,'gif','WriteMode','append','DelayTime',2/fps); end
    end
    if n == round(N1/2), imwrite(img, fullfile(outdir,'frame_955rpm.png')); end
    if n == iPk,         imwrite(img, fullfile(outdir,'frame_sweep_resonance.png')); end
end
close(vw);
fprintf('Animation saved: %s, %s (%d frames)\n', vw.Filename, 'DVA_animation.gif', N);
fprintf('Key frames: frame_955rpm.png, frame_sweep_resonance.png\n');

%% ===================== local functions ==================================
function [X0,X1,X2] = cplx(m1,k1,c1,m2,k2,c2,F0,w)
% complex steady-state amplitudes for force F0*exp(i w t) on m1
    X0  = F0/(k1 - m1*w^2 + 1i*w*c1);
    Z11 = k1 + k2 - m1*w^2 + 1i*w*(c1+c2);
    Z12 = -(k2 + 1i*w*c2);
    Z22 = k2 - m2*w^2 + 1i*w*c2;
    D   = Z11*Z22 - Z12^2;
    X1  = F0*Z22/D;   X2 = -F0*Z12/D;
end

function P = mkParts(ax, steel, k, withDVA, dva)
    P.s1 = plot(ax,NaN,NaN,'Color',k,'LineWidth',1.8);       % primary spring
    P.d1 = plot(ax,NaN,NaN,'Color',k,'LineWidth',1.8);       % primary damper
    P.m  = patch(ax,NaN,NaN,steel,'EdgeColor',k,'LineWidth',1.5);
    P.disc = plot(ax,NaN,NaN,'w','LineWidth',2);
    P.ecc  = plot(ax,NaN,NaN,'o','MarkerSize',10,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','w');
    P.lab  = text(ax,0,0,'Pump 200 kg','Color','w','FontWeight','bold','HorizontalAlignment','center','FontSize',10);
    if withDVA
        P.s2 = plot(ax,NaN,NaN,'Color',k,'LineWidth',1.8);
        P.d2 = plot(ax,NaN,NaN,'Color',k,'LineWidth',1.8);
        P.a  = patch(ax,NaN,NaN,dva,'EdgeColor',k,'LineWidth',1.5);
        P.alab = text(ax,0,0,'DVA','Color','k','FontWeight','bold','HorizontalAlignment','center','FontSize',10);
    end
end

function upParts(P, d1, d2, th, withDVA)
% d1, d2 = magnified displacements (plot units) of machine and absorber
    y0b = 200; Hm = 120; Wm = 230; y0a = y0b + Hm + 110; Ha = 60; Wa = 130;
    yb = y0b + d1;
    [x,y] = springXY(-55, 0, yb, 7, 22);   set(P.s1,'XData',x,'YData',y);
    [x,y] = dashpotXY(55, 0, yb, y0b, 20); set(P.d1,'XData',x,'YData',y);
    set(P.m,'XData',[-Wm Wm Wm -Wm]/2,'YData',[yb yb yb+Hm yb+Hm]);
    yc = yb + Hm/2 - 8;  r = 30;  a = linspace(0,2*pi,60);
    set(P.disc,'XData',r*cos(a),'YData',yc + r*sin(a));
    set(P.ecc,'XData',r*sin(th),'YData',yc + r*cos(th));      % vertical force ~ cos(th)
    set(P.lab,'Position',[0 yb+Hm-14 0]);
    if withDVA
        ya = y0a + d2;  yt = yb + Hm;  L0 = y0a - (y0b + Hm);
        [x,y] = springXY(-30, yt, ya, 5, 16);   set(P.s2,'XData',x,'YData',y);
        [x,y] = dashpotXY(30, yt, ya, L0, 12);  set(P.d2,'XData',x,'YData',y);
        set(P.a,'XData',[-Wa Wa Wa -Wa]/2,'YData',[ya ya ya+Ha ya+Ha]);
        set(P.alab,'Position',[0 ya+Ha/2 0]);
    end
end

function [x,y] = springXY(xc, y0, y1, nc, wd)
    L = y1 - y0;  lead = 0.12*L;
    yz = linspace(y0+lead, y1-lead, 2*nc+1);
    xz = xc + wd*[0, repmat([1 -1],1,nc)];  xz(end) = xc;
    x = [xc xz xc];  y = [y0 yz y1];
end

function [x,y] = dashpotXY(xc, y0, y1, L0, w)
    yc = y0 + 0.10*L0;  yt = y0 + 0.85*L0;  yp = y1 - 0.55*L0;
    x = [xc xc NaN xc-w xc-w xc+w xc+w NaN xc xc NaN xc-0.8*w xc+0.8*w];
    y = [y0 yc NaN yt yc yc yt NaN y1 yp NaN yp yp];
end

function setReadout(h, X, w, ISO)
    v = abs(X)*w/sqrt(2)*1e3;  z = 'ABCD';  iz = find([v <= ISO, true], 1);
    col = [0.1 0.55 0.15; 0.45 0.6 0.1; 0.9 0.5 0.0; 0.8 0.1 0.1];
    set(h,'String',{sprintf('X_1 = %.3f mm', abs(X)*1e3), sprintf('v_{rms} = %.1f mm/s  \\rightarrow  Zone %s', v, z(iz))}, ...
          'Color',col(iz,:));
end

function img = fitTo(img, sz)
% keep every video frame the same size
    img = img(1:min(end,sz(1)), 1:min(end,sz(2)), :);
    if size(img,1) < sz(1) || size(img,2) < sz(2)
        img(sz(1),sz(2),3) = 0;
    end
end

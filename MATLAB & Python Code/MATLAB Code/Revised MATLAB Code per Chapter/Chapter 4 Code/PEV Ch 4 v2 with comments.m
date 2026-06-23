%% ================================================================
% PEV Chapter 4 Code v2 with Comments.m
%
% Chapter 4: Finite-Horizon Model-Based Pursuit-Evasion (only)  
%
% This script computes the Riccati saddle gains, evaluates three evader policies on
% the nonlinear actuator-limited kinematic-bicycle rollout, and exports only
% Chapter 4 figures and numerical summaries.
%% ================================================================

clear;
clc;
close all;

rng(11,'twister'); % A fixed seed keeps the Chapter 4 nominal benchmark reproducible when the figures are regenerated.

%% ================================================================
% 0. THESIS FIGURE AND LEGEND DEFAULTS
% White figure/axes backgrounds and black text match the Chapter 6 export.
%% ================================================================
set(groot,'defaultFigureColor','w');
set(groot,'defaultFigureInvertHardcopy','off');
set(groot,'defaultAxesColor','w');
set(groot,'defaultAxesXColor','k');
set(groot,'defaultAxesYColor','k');
set(groot,'defaultTextColor','k');
set(groot,'defaultLegendTextColor','k');
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
set(groot,'defaultAxesFontSize',12);
set(groot,'defaultTextFontSize',12);

%% ================================================================
% 1. TIME, MODEL, AND STATE DIMENSIONS
%% ================================================================

Ts = 0.02;
N  = 400;

Lp = 0.256;
Le = 0.256;

nx  = 8;
nup = 2;
nue = 2;
nz  = nx + nup + nue;

% State:
% This is the eight-state two-vehicle construction introduced in Chapter 3 and used consistently through Chapter 6.
% X = [xP yP psiP vP xE yE psiE vE]'

%% ================================================================
% 2. INITIAL CONDITION AND NOMINAL OPERATING POINT
%% ================================================================


% The same initial engagement is retained across the nominal Chapters 4--5 comparison so the results stay comparable.

x0 = [ ...
    -2.00; ...   % xP
     0.00; ...   % yP
     0.00; ...   % psiP
     0.75; ...   % vP
     0.00; ...   % xE
     0.60; ...   % yE
     0.00; ...   % psiE
     0.60];      % vE

vbarP = 0.75;
vbarE = 0.60;

xbar = [ ...
    0;
    0;
    0;
    vbarP;
    0;
    0;
    0;
    vbarE];

captureRadius = 0.35;


% These limits are enforced during the nonlinear rollout, so the Chapter 4 results reflect bounded vehicle commands rather than an ideal unconstrained game.

%% ================================================================
% 3. PHYSICAL LIMITS
%% ================================================================

amaxP = 1.25;
amaxE = 0.75;

deltaMaxP = deg2rad(27);
deltaMaxE = deg2rad(22);

vminP = 0.0;
vmaxP = 1.45;

vminE = 0.0;
vmaxE = 1.05;

evaderPolicyScale = 0.55;

% Chapter 4 synthesizes the local saddle policy from this linear model; the nonlinear bicycle equations are used later for the actual rollout.

%% ================================================================
% 4. LINEARIZED KINEMATIC BICYCLE MODEL
%% ================================================================


% The discrete matrices below are the finite-horizon game model used by the Riccati recursion in Chapter 4.

AcP = [ ...
    0  0  0      1;
    0  0  vbarP  0;
    0  0  0      0;
    0  0  0      0];

BcP = [ ...
    0  0;
    0  0;
    0  vbarP/Lp;
    1  0];

AcE = [ ...
    0  0  0      1;
    0  0  vbarE  0;
    0  0  0      0;
    0  0  0      0];

BcE = [ ...
    0  0;
    0  0;
    0  vbarE/Le;
    1  0];

%% ================================================================
% 5. ZERO-ORDER-HOLD DISCRETIZATION
%% ================================================================

MP = expm([AcP BcP; zeros(2,6)]*Ts);
AP = MP(1:4,1:4);
BP = MP(1:4,5:6);

ME = expm([AcE BcE; zeros(2,6)]*Ts);
AE = ME(1:4,1:4);
BE = ME(1:4,5:6);

A = blkdiag(AP,AE);

Bp = [ ...
    BP;
    zeros(4,2)];

Be = [ ...
    zeros(4,2);
    BE];

fprintf('\n================================================\n');
fprintf('MODEL INFORMATION\n');
fprintf('================================================\n');
fprintf('A dimensions      = %d x %d\n',size(A,1),size(A,2));
fprintf('Bp dimensions     = %d x %d\n',size(Bp,1),size(Bp,2));
fprintf('Be dimensions     = %d x %d\n',size(Be,1),size(Be,2));
fprintf('Sampling time     = %.3f s\n',Ts);
fprintf('Horizon time      = %.3f s\n',N*Ts);
fprintf('Initial distance  = %.4f m\n',norm(x0(1:2)-x0(5:6)));
fprintf('Capture radius    = %.4f m\n\n',captureRadius);

%% ================================================================
% 6. COST STRUCTURE
%% ================================================================

% The cost is written in relative vehicle coordinates so the game rewards closing position, heading, and speed differences rather than tracking an arbitrary world-frame point.

Crel = [eye(4), -eye(4)];

QrBase = diag([35,35,5,1]);

Qbase = Crel'*QrBase*Crel + 1e-6*eye(nx);

RpBase = diag([0.35,0.18]);
ReBase = diag([8.0,8.0]);

%% ================================================================
% 7. SEARCH FOR A VALID SADDLE GAME
%% ================================================================


% This small search is a numerical well-posedness check, not controller tuning after seeing results. A candidate is retained only when the stage-wise saddle conditions are satisfied.

qScaleCandidates   = [1.25, 1.00, 0.80, 0.60];
terminalCandidates = [12, 10, 8, 6, 4];
rpScaleCandidates  = [1.00, 1.25, 1.50, 2.00, 2.50];
reScaleCandidates  = [1.00, 1.25, 1.50, 1.75, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00, 6.00, 8.00];

validGameFound = false;

selectedQScale = NaN;
selectedTerminal = NaN;
selectedRpScale = NaN;
selectedReScale = NaN;

for qScale = qScaleCandidates

    Qtrial = qScale*Qbase;

    for terminalMultiplier = terminalCandidates

        QfTrial = terminalMultiplier*Qtrial;

        for rpScale = rpScaleCandidates

            RpTrial = rpScale*RpBase;

            for reScale = reScaleCandidates

                ReTrial = reScale*ReBase;

                Ptrial  = zeros(nx,nx,N+1);
                KpTrial = zeros(nup,nx,N);
                KeTrial = zeros(nue,nx,N);

                Ptrial(:,:,N+1) = QfTrial;

                candidateValid = true;

                minPcurvCandidate = inf;
                maxEschurCandidate = -inf;
                maxCondCandidate = 0;

                % The backward recursion creates a stage-dependent gain sequence because the amount of time remaining changes over the finite-horizon.
                for k = N:-1:1

                    Pnext = Ptrial(:,:,k+1);

                    Gpp = RpTrial + Bp'*Pnext*Bp;
                    Gpe =           Bp'*Pnext*Be;
                    Gep =           Be'*Pnext*Bp;
                    Gee = -ReTrial + Be'*Pnext*Be;

                    GppSym = (Gpp+Gpp')/2;
                    minEigP = min(eig(GppSym));

                    if minEigP <= 1e-9
                        candidateValid = false;
                        break;
                    end

                    SchurE = Gee - Gep*(Gpp\Gpe);
                    SchurESym = (SchurE+SchurE')/2;
                    maxEigE = max(eig(SchurESym));

                    if maxEigE >= -1e-9
                        candidateValid = false;
                        break;
                    end

                    G = [Gpp Gpe; Gep Gee];

                    if rcond(G) < 1e-11
                        candidateValid = false;
                        break;
                    end

                    F = [ ...
                        Bp'*Pnext*A;
                        Be'*Pnext*A];

                    Kjoint = G\F;

                    KpTrial(:,:,k) = Kjoint(1:nup,:);
                    KeTrial(:,:,k) = Kjoint(nup+1:end,:);

                    Pk = Qtrial + A'*Pnext*A - F'*(G\F);
                    Ptrial(:,:,k) = (Pk+Pk')/2;

                    minPcurvCandidate = min(minPcurvCandidate,minEigP);
                    maxEschurCandidate = max(maxEschurCandidate,maxEigE);
                    maxCondCandidate = max(maxCondCandidate,cond(G));

                end

                if candidateValid

                    validGameFound = true;

                    Q = Qtrial;
                    Qf = QfTrial;
                    Rp = RpTrial;
                    Re = ReTrial;

                    Pexact = Ptrial;
                    KpExact = KpTrial;
                    KeExact = KeTrial;

                    selectedQScale = qScale;
                    selectedTerminal = terminalMultiplier;
                    selectedRpScale = rpScale;
                    selectedReScale = reScale;

                    bestMinPcurv = minPcurvCandidate;
                    bestMaxEschur = maxEschurCandidate;
                    bestMaxCond = maxCondCandidate;

                    break;

                end

            end

            if validGameFound
                break;
            end

        end

        if validGameFound
            break;
        end

    end

    if validGameFound
        break;
    end

end

if ~validGameFound
    error('No valid saddle game found. Increase Rp/Re or reduce Q/Qf.');
end

fprintf('================================================\n');
fprintf('VALID SADDLE GAME FOUND\n');
fprintf('================================================\n');
fprintf('Q scale                   = %.3f\n',selectedQScale);
fprintf('Terminal multiplier       = %.3f\n',selectedTerminal);
fprintf('Rp scale                  = %.3f\n',selectedRpScale);
fprintf('Re scale                  = %.3f\n',selectedReScale);
fprintf('Minimum Pursuer Curvature = %.6e\n',bestMinPcurv);
fprintf('Maximum Evader Schur Eig. = %.6e\n',bestMaxEschur);
fprintf('Maximum Saddle Condition  = %.6e\n\n',bestMaxCond);

% All three cases use the same pursuer gain schedule. Only the evader response changes, which isolates how adversarial behavior changes the nominal engagement.

%% ================================================================
% 8. MODEL-BASED NOMINAL EVADER-POLICY SWEEP
%% ================================================================

evaderModes = ["straight","limited_saddle","saddle"];
displayLabels = ["Straight","Limited Saddle","Full Saddle"];
numModes = numel(evaderModes);

capturedAll = false(numModes,1);
captureTimeAll = nan(numModes,1);
minDistanceAll = nan(numModes,1);
finalDistanceAll = nan(numModes,1);
costAll = nan(numModes,1);
pursuerEffortAll = nan(numModes,1);
evaderEffortAll = nan(numModes,1);
pursuerSatAll = zeros(numModes,1);
evaderSatAll = zeros(numModes,1);
stepsUsedAll = N*ones(numModes,1);

Xall  = zeros(nx,N+1,numModes);
UPall = zeros(nup,N,numModes);
UEall = zeros(nue,N,numModes);

bike = @(z,u,L) [ ...
    z(4)*cos(z(3)); ...
    z(4)*sin(z(3)); ...
    z(4)/L*tan(u(2)); ...
    u(1)];

for modeIndex = 1:numModes

    X = zeros(nx,N+1);
    UP = zeros(nup,N);
    UE = zeros(nue,N);
    X(:,1) = x0;

    captured = false;
    captureStep = N;
    cumulativeCost = 0;
    pursuerEffort = 0;
    evaderEffort = 0;
    pursuerSat = 0;
    evaderSat = 0;

    for k = 1:N

        currentState = X(:,k);
        errorState = currentState - xbar;

        upRaw = -KpExact(:,:,k)*errorState;

        switch evaderModes(modeIndex)
            case "straight"
                ueRaw = [0;0];
            case "limited_saddle"
                ueRaw = -evaderPolicyScale*KeExact(:,:,k)*errorState;
            case "saddle"
                ueRaw = -KeExact(:,:,k)*errorState;
            otherwise
                error('Unknown evader mode.');
        end

        % The Riccati policy first produces an unconstrained mathematical command. It is then clipped to the feasible acceleration and steering limits before propagation.
        up = upRaw;
        ue = ueRaw;

        up(1) = min(max(up(1),-amaxP),amaxP);
        up(2) = min(max(up(2),-deltaMaxP),deltaMaxP);
        ue(1) = min(max(ue(1),-amaxE),amaxE);
        ue(2) = min(max(ue(2),-deltaMaxE),deltaMaxE);

        if any(abs(up-upRaw) > 1e-10)
            pursuerSat = pursuerSat + 1;
        end
        if any(abs(ue-ueRaw) > 1e-10)
            evaderSat = evaderSat + 1;
        end

        UP(:,k) = up;
        UE(:,k) = ue;

        zP = currentState(1:4);
        zE = currentState(5:8);

        % The policy was synthesized from the local linear model, but performance is evaluated here on the nonlinear kinematic bicycle dynamics using RK4 integration.
        % Chapter 4 evaluates the Riccati controller on the nonlinear kinematic bicycle model rather than only on the linearized synthesis model. 
        % RK4 integration is used here to propagate each vehicle over one sample interval using four slope evaluations, which provides a more accurate approximation of curved, steering-dependent vehicle motion than a single forward-Euler update.

        k1P = bike(zP,up,Lp);
        k2P = bike(zP+0.5*Ts*k1P,up,Lp);
        k3P = bike(zP+0.5*Ts*k2P,up,Lp);
        k4P = bike(zP+Ts*k3P,up,Lp);
        zPnext = zP + Ts/6*(k1P+2*k2P+2*k3P+k4P);

        k1E = bike(zE,ue,Le);
        k2E = bike(zE+0.5*Ts*k1E,ue,Le);
        k3E = bike(zE+0.5*Ts*k2E,ue,Le);
        k4E = bike(zE+Ts*k3E,ue,Le);
        zEnext = zE + Ts/6*(k1E+2*k2E+2*k3E+k4E);

        zPnext(4) = min(max(zPnext(4),vminP),vmaxP);
        zEnext(4) = min(max(zEnext(4),vminE),vmaxE);
        zPnext(3) = atan2(sin(zPnext(3)),cos(zPnext(3)));
        zEnext(3) = atan2(sin(zEnext(3)),cos(zEnext(3)));

        nextState = [zPnext;zEnext];
        X(:,k+1) = nextState;

        cumulativeCost = cumulativeCost + errorState'*Q*errorState + up'*Rp*up - ue'*Re*ue;
        pursuerEffort = pursuerEffort + Ts*(up'*up);
        evaderEffort = evaderEffort + Ts*(ue'*ue);

        % Capture is a sampled event: the rollout ends at the first time the planar separation enters the prescribed capture radius.
        separation = norm(zPnext(1:2)-zEnext(1:2));
        if separation <= captureRadius
            captured = true;
            captureStep = k;
            break;
        end
    end

    used = 1:(captureStep+1);
    distance = sqrt((X(1,used)-X(5,used)).^2 + (X(2,used)-X(6,used)).^2);

    capturedAll(modeIndex) = captured;
    if captured
        captureTimeAll(modeIndex) = captureStep*Ts;
    else
        captureTimeAll(modeIndex) = N*Ts;
    end
    minDistanceAll(modeIndex) = min(distance);
    finalDistanceAll(modeIndex) = distance(end);
    costAll(modeIndex) = cumulativeCost;
    pursuerEffortAll(modeIndex) = pursuerEffort;
    evaderEffortAll(modeIndex) = evaderEffort;
    pursuerSatAll(modeIndex) = pursuerSat;
    evaderSatAll(modeIndex) = evaderSat;
    stepsUsedAll(modeIndex) = captureStep;

    Xall(:,:,modeIndex) = X;
    UPall(:,:,modeIndex) = UP;
    UEall(:,:,modeIndex) = UE;
end

fprintf('================================================\n');
fprintf('CHAPTER 4: MODEL-BASED EVADER-POLICY RESULTS\n');
fprintf('================================================\n');
for modeIndex = 1:numModes
    fprintf('%s: captured=%d, t_c=%.2f s, d_min=%.4f m, margin=%.1f mm, pursuer saturation=%.1f%%\n', ...
        displayLabels(modeIndex), capturedAll(modeIndex), captureTimeAll(modeIndex), ...
        minDistanceAll(modeIndex), 1000*(captureRadius-minDistanceAll(modeIndex)), ...
        100*pursuerSatAll(modeIndex)/max(stepsUsedAll(modeIndex),1));
end




% The exported trajectories and summary plots are the source figures for the nominal-results discussion in Chapter 4.

%% ================================================================
% 9. CHAPTER 4 FIGURE EXPORT
%% ================================================================

scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    scriptDir = pwd;
else
    scriptDir = fileparts(scriptPath);
end
outputFolder = fullfile(scriptDir,'PEV_Ch4_Model_Based_Output_v11');
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

figSize = [100 100 1050 700];
exportResolution = 300;

% Figures 4.1--4.3: trajectory figures.
% The legend identifies the model-based game objects and the sampled
% capture event represented in each trajectory figure.


for modeIndex = 1:numModes
    finalStep = stepsUsedAll(modeIndex);
    idx = 1:(finalStep+1);

    fig = figure('Color','w','InvertHardcopy','off','Position',figSize);

    hP = plot(Xall(1,idx,modeIndex),Xall(2,idx,modeIndex), ...
        '-','LineWidth',1.9);
    hold on;
    hE = plot(Xall(5,idx,modeIndex),Xall(6,idx,modeIndex), ...
        '--','LineWidth',1.9);

    hP0 = plot(x0(1),x0(2),'o','MarkerSize',9,'LineWidth',1.6);
    hE0 = plot(x0(5),x0(6),'s','MarkerSize',9,'LineWidth',1.6);

    capturePX = Xall(1,idx(end),modeIndex);
    capturePY = Xall(2,idx(end),modeIndex);
    captureEX = Xall(5,idx(end),modeIndex);
    captureEY = Xall(6,idx(end),modeIndex);

    hPc = plot(capturePX,capturePY,'x','MarkerSize',11,'LineWidth',2.0);
    hEc = plot(captureEX,captureEY,'x','MarkerSize',11,'LineWidth',2.0);

    theta = linspace(0,2*pi,300);
    hCap = plot(captureEX + captureRadius*cos(theta), ...
                captureEY + captureRadius*sin(theta), ...
                ':','LineWidth',1.5);

    xlabel('x position [m]','Color','k');
    ylabel('y position [m]','Color','k');
    title(sprintf('Nonlinear Kinematic-Bicycle Trajectories: %s Evader', ...
        displayLabels(modeIndex)),'Color','k');

    axis equal;
    grid on;
    box on;
    set(gca,'Color','w','XColor','k','YColor','k', ...
        'FontName','Times New Roman','FontSize',12,'LineWidth',1.0);

    lgd = legend([hP hE hP0 hE0 hPc hEc hCap], ...
        {'Pursuer trajectory', ...
         'Evader trajectory', ...
         'Pursuer initial position', ...
         'Evader initial position', ...
         'Pursuer position at first sampled capture', ...
         'Evader position at first sampled capture', ...
         'Capture boundary'}, ...
        'Location','best');
    set(lgd,'Color','w','TextColor','k','EdgeColor','k', ...
        'FontName','Times New Roman','FontSize',10.5);

    set(fig,'PaperPositionMode','auto');
    print(fig,fullfile(outputFolder,sprintf('fig4_%d_trajectory_%s.png', ...
        modeIndex,lower(strrep(char(displayLabels(modeIndex)),' ', '_')))), ...
        '-dpng',sprintf('-r%d',exportResolution));
    close(fig);
end

% Saturation is reported as a realization diagnostic. It shows how often the bounded rollout differs from the raw local-game request, which motivates the execution-side study in Chapter 6.
orderedCategories = categorical(displayLabels,displayLabels,'Ordinal',true);
captureMarginMM = 1000*(captureRadius-minDistanceAll);
pursuerSaturationPercent = 100*pursuerSatAll./max(stepsUsedAll,1);

% Figures 4.4--4.6 contain one model-based series. The single legend entry
% makes that design explicit without implying a learned-policy comparison.
fig = figure('Color','w','InvertHardcopy','off','Position',figSize);
hBar = bar(orderedCategories,captureTimeAll);
ylabel('Capture time [s]','Color','k');
title('Capture Time Across Evader Policies','Color','k');
grid on;
box on;
set(gca,'Color','w','XColor','k','YColor','k', ...
    'FontName','Times New Roman','FontSize',12,'LineWidth',1.0);
lgd = legend(hBar,'Model-based saddle-policy rollout','Location','northwest');
set(lgd,'Color','w','TextColor','k','EdgeColor','k', ...
    'FontName','Times New Roman','FontSize',10.5);
set(fig,'PaperPositionMode','auto');
print(fig,fullfile(outputFolder,'fig4_4_capture_time.png'), ...
    '-dpng',sprintf('-r%d',exportResolution));
close(fig);

fig = figure('Color','w','InvertHardcopy','off','Position',figSize);
hBar = bar(orderedCategories,captureMarginMM);
ylabel('Sampled capture margin [mm]','Color','k');
title('Sampled Capture Margin Across Evader Policies','Color','k');
grid on;
box on;
set(gca,'Color','w','XColor','k','YColor','k', ...
    'FontName','Times New Roman','FontSize',12,'LineWidth',1.0);
lgd = legend(hBar,'Model-based saddle-policy rollout','Location','northwest');
set(lgd,'Color','w','TextColor','k','EdgeColor','k', ...
    'FontName','Times New Roman','FontSize',10.5);
set(fig,'PaperPositionMode','auto');
print(fig,fullfile(outputFolder,'fig4_5_capture_margin.png'), ...
    '-dpng',sprintf('-r%d',exportResolution));
close(fig);

fig = figure('Color','w','InvertHardcopy','off','Position',figSize);
hBar = bar(orderedCategories,pursuerSaturationPercent);
ylabel('Pursuer saturation [% of pre-capture steps]','Color','k');
title('Pursuer Input Saturation Across Evader Policies','Color','k');
grid on;
box on;
set(gca,'Color','w','XColor','k','YColor','k', ...
    'FontName','Times New Roman','FontSize',12,'LineWidth',1.0);
lgd = legend(hBar,'Model-based saddle-policy rollout','Location','northwest');
set(lgd,'Color','w','TextColor','k','EdgeColor','k', ...
    'FontName','Times New Roman','FontSize',10.5);
set(fig,'PaperPositionMode','auto');
print(fig,fullfile(outputFolder,'fig4_6_pursuer_saturation.png'), ...
    '-dpng',sprintf('-r%d',exportResolution));
close(fig);

% The CSV preserves the numerical values behind the Chapter 4 figures and tables for straightforward thesis traceability.
summaryTable = table(displayLabels',capturedAll,captureTimeAll,minDistanceAll,captureMarginMM, ...
    pursuerSaturationPercent,pursuerEffortAll,evaderEffortAll,costAll, ...
    'VariableNames',{'EvaderPolicy','Captured','CaptureTime_s','MinimumSeparation_m', ...
    'SampledCaptureMargin_mm','PursuerSaturation_percent','PursuerEffort','EvaderEffort','GameCost'});
writetable(summaryTable,fullfile(outputFolder,'chapter4_model_based_summary.csv'));

fid = fopen(fullfile(outputFolder,'run_manifest.txt'),'w');
if fid ~= -1
    fprintf(fid,'Chapter 4 model-based-only run manifest\n');
    fprintf(fid,'This script reports model-based Chapter 4 results only.\n');
    fprintf(fid,'Seed: rng(11,''twister'')\n');
    fprintf(fid,'Ts=%.6f s, N=%d, capture radius=%.6f m\n',Ts,N,captureRadius);
    fprintf(fid,'Evader policies: straight, limited_saddle (alpha=%.2f), saddle.\n',evaderPolicyScale);
    fclose(fid);
end

fprintf('\nChapter 4 v11 model-based figures and CSV summary saved to:\n%s\n',outputFolder);

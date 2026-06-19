%% ================================================================
%  COMPLETE PEV NOMINAL SWEEP + THESIS FIGURE EXPORT
%
%  Functions within the MATLAB Command Window.
%
%  This program does the following:
%    1. Builds the two-vehicle kinematic-bicycle LQ game.
%    2. Searches for a valid finite-horizon saddle-point game.
%    3. Recovers the Riccati gains using fitted quadratic Q-learning.
%    4. Sweeps evader modes: straight, limited_saddle, saddle.
%    5. Prints numerical results.
%    6. Generates thesis-ready figures.
%    7. Saves all figures as PNG files.
%    8. Saves CSV summaries, a LaTeX thesis table, and a run manifest.
%% ================================================================

clear;
clc;
close all;

rng(11,'twister');

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
% X = [xP yP psiP vP xE yE psiE vE]'

%% ================================================================
% 2. INITIAL CONDITION AND NOMINAL OPERATING POINT
%% ================================================================

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

%% ================================================================
% 4. LINEARIZED KINEMATIC BICYCLE MODEL
%% ================================================================

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

Crel = [eye(4), -eye(4)];

QrBase = diag([35,35,5,1]);

Qbase = Crel'*QrBase*Crel + 1e-6*eye(nx);

RpBase = diag([0.35,0.18]);
ReBase = diag([8.0,8.0]);

%% ================================================================
% 7. SEARCH FOR A VALID SADDLE GAME
%% ================================================================

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

%% ================================================================
% 8. FITTED QUADRATIC Q-LEARNING RECOVERY
%% ================================================================

nFeatures = nz*(nz+1)/2;

samplesPerStage = max(850,7*nFeatures);
lambda = 1e-10;

Plearn  = zeros(nx,nx,N+1);
KpLearn = zeros(nup,nx,N);
KeLearn = zeros(nue,nx,N);

bellmanResidual = zeros(N,1);
featureRank = zeros(N,1);
featureCondition = zeros(N,1);

Plearn(:,:,N+1) = Qf;

for k = N:-1:1

    Pnext = Plearn(:,:,k+1);

    Phi = zeros(samplesPerStage,nFeatures);
    Y   = zeros(samplesPerStage,1);

    for sample = 1:samplesPerStage

        x = [ ...
            1.20*randn;
            1.20*randn;
            0.30*randn;
            0.30*randn;
            1.20*randn;
            1.20*randn;
            0.30*randn;
            0.30*randn];

        up = [ ...
            0.45*randn;
            0.18*randn];

        ue = [ ...
            0.40*randn;
            0.16*randn];

        xnext = A*x + Bp*up + Be*ue;

        target = ...
            x'*Q*x + ...
            up'*Rp*up - ...
            ue'*Re*ue + ...
            xnext'*Pnext*xnext;

        zeta = [x;up;ue];

        phi = zeros(nFeatures,1);
        idx = 1;

        for row = 1:nz
            for col = row:nz

                if row == col
                    phi(idx) = zeta(row)*zeta(col);
                else
                    phi(idx) = 2*zeta(row)*zeta(col);
                end

                idx = idx + 1;

            end
        end

        Phi(sample,:) = phi';
        Y(sample) = target;

    end

    featureRank(k) = rank(Phi);
    featureCondition(k) = cond(Phi'*Phi + lambda*eye(nFeatures));

    if featureRank(k) < nFeatures
        error('Feature matrix rank deficient at stage %d.',k);
    end

    theta = (Phi'*Phi + lambda*eye(nFeatures))\(Phi'*Y);

    H = zeros(nz,nz);
    idx = 1;

    for row = 1:nz
        for col = row:nz

            H(row,col) = theta(idx);
            H(col,row) = theta(idx);

            idx = idx + 1;

        end
    end

    H = (H+H')/2;

    Hxx = H(1:nx,1:nx);

    pIdx = nx+1:nx+nup;
    eIdx = nx+nup+1:nz;

    Hpp = H(pIdx,pIdx);
    Hpe = H(pIdx,eIdx);
    Hep = H(eIdx,pIdx);
    Hee = H(eIdx,eIdx);

    Hpx = H(pIdx,1:nx);
    Hex = H(eIdx,1:nx);

    Hxp = H(1:nx,pIdx);
    Hxe = H(1:nx,eIdx);

    Huu = [Hpp Hpe; Hep Hee];
    Hux = [Hpx; Hex];
    Hxu = [Hxp Hxe];

    HppSym = (Hpp+Hpp')/2;
    SchurE = Hee - Hep*(Hpp\Hpe);
    SchurESym = (SchurE+SchurE')/2;

    if min(eig(HppSym)) <= 1e-8
        error('Learned pursuer curvature invalid at stage %d.',k);
    end

    if max(eig(SchurESym)) >= -1e-8
        error('Learned evader Schur condition invalid at stage %d.',k);
    end

    KjointLearn = Huu\Hux;

    KpLearn(:,:,k) = KjointLearn(1:nup,:);
    KeLearn(:,:,k) = KjointLearn(nup+1:end,:);

    PkLearn = Hxx - Hxu*(Huu\Hux);
    Plearn(:,:,k) = (PkLearn+PkLearn')/2;

    prediction = Phi*theta;

    bellmanResidual(k) = norm(prediction-Y)/max(norm(Y),eps);

end

pursuerGainError = norm(KpLearn(:)-KpExact(:))/max(norm(KpExact(:)),eps);
evaderGainError  = norm(KeLearn(:)-KeExact(:))/max(norm(KeExact(:)),eps);

fprintf('================================================\n');
fprintf('FITTED-Q RECOVERY RESULTS\n');
fprintf('================================================\n');
fprintf('Augmented Q dimension       = %d\n',nz);
fprintf('Unique Q coefficients       = %d\n',nFeatures);
fprintf('Samples per stage           = %d\n',samplesPerStage);
fprintf('Relative pursuer gain error = %.6e\n',pursuerGainError);
fprintf('Relative evader gain error  = %.6e\n',evaderGainError);
fprintf('Maximum Bellman residual    = %.6e\n',max(bellmanResidual));
fprintf('Minimum feature rank        = %d of %d\n',min(featureRank),nFeatures);
fprintf('Maximum feature condition   = %.6e\n\n',max(featureCondition));

%% ================================================================
% 9. NOMINAL EVADER-MODE SWEEP
%% ================================================================

evaderModes = ["straight","limited_saddle","saddle"];
numModes = numel(evaderModes);

controllerNames = ["Exact LQ","Learned Q"];
numControllers = 2;

capturedAll = false(numModes,numControllers);
captureTimeAll = nan(numModes,numControllers);
minDistanceAll = nan(numModes,numControllers);
finalDistanceAll = nan(numModes,numControllers);
costAll = nan(numModes,numControllers);
pursuerEffortAll = nan(numModes,numControllers);
evaderEffortAll = nan(numModes,numControllers);
pursuerSatAll = zeros(numModes,numControllers);
evaderSatAll = zeros(numModes,numControllers);
stepsUsedAll = N*ones(numModes,numControllers);

Xall  = zeros(nx,N+1,numModes,numControllers);
UPall = zeros(nup,N,numModes,numControllers);
UEall = zeros(nue,N,numModes,numControllers);

timeVector = 0:Ts:N*Ts;

bike = @(z,u,L) [ ...
    z(4)*cos(z(3));
    z(4)*sin(z(3));
    z(4)/L*tan(u(2));
    u(1)];

for modeIndex = 1:numModes

    evaderModeCurrent = evaderModes(modeIndex);

    for controller = 1:numControllers

        X = zeros(nx,N+1);
        UP = zeros(nup,N);
        UE = zeros(nue,N);

        X(:,1) = x0;

        captured = false;
        captureStep = NaN;
        stepsUsed = N;

        cumulativeCost = 0;
        pursuerEffort = 0;
        evaderEffort = 0;

        pursuerSat = 0;
        evaderSat = 0;

        for k = 1:N

            currentState = X(:,k);
            errorState = currentState - xbar;

            if controller == 1
                KpCurrent = KpExact(:,:,k);
                KeCurrent = KeExact(:,:,k);
            else
                KpCurrent = KpLearn(:,:,k);
                KeCurrent = KeLearn(:,:,k);
            end

            upRaw = -KpCurrent*errorState;

            switch evaderModeCurrent

                case "straight"
                    ueRaw = [0;0];

                case "limited_saddle"
                    ueRaw = -evaderPolicyScale*KeCurrent*errorState;

                case "saddle"
                    ueRaw = -KeCurrent*errorState;

                otherwise
                    error('Unknown evader mode.');
            end

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

            cumulativeCost = cumulativeCost + ...
                errorState'*Q*errorState + ...
                up'*Rp*up - ...
                ue'*Re*ue;

            pursuerEffort = pursuerEffort + Ts*(up'*up);
            evaderEffort  = evaderEffort  + Ts*(ue'*ue);

            separation = norm(zPnext(1:2)-zEnext(1:2));

            if separation <= captureRadius

                captured = true;
                captureStep = k;
                stepsUsed = k;
                break;

            end

        end

        distance = sqrt((X(1,1:stepsUsed+1)-X(5,1:stepsUsed+1)).^2 + ...
                        (X(2,1:stepsUsed+1)-X(6,1:stepsUsed+1)).^2);

        capturedAll(modeIndex,controller) = captured;

        if captured
            captureTimeAll(modeIndex,controller) = captureStep*Ts;
        end

        minDistanceAll(modeIndex,controller) = min(distance);
        finalDistanceAll(modeIndex,controller) = distance(end);
        costAll(modeIndex,controller) = cumulativeCost;
        pursuerEffortAll(modeIndex,controller) = pursuerEffort;
        evaderEffortAll(modeIndex,controller) = evaderEffort;
        pursuerSatAll(modeIndex,controller) = pursuerSat;
        evaderSatAll(modeIndex,controller) = evaderSat;
        stepsUsedAll(modeIndex,controller) = stepsUsed;

        Xall(:,:,modeIndex,controller) = X;
        UPall(:,:,modeIndex,controller) = UP;
        UEall(:,:,modeIndex,controller) = UE;

    end

end

%% ================================================================
% 10. PRINT RESULTS
%% ================================================================

fprintf('================================================\n');
fprintf('DISTINGUISHABLE NOMINAL EVADER-MODE SWEEP\n');
fprintf('================================================\n');

for modeIndex = 1:numModes

    fprintf('\nEvader mode: %s\n',evaderModes(modeIndex));
    fprintf('------------------------------------------------\n');

    for controller = 1:numControllers

        fprintf('%s controller\n',controllerNames(controller));
        fprintf('  Captured                 = %d\n',capturedAll(modeIndex,controller));

        if capturedAll(modeIndex,controller)
            fprintf('  Capture time             = %.4f s\n',captureTimeAll(modeIndex,controller));
        else
            fprintf('  Capture time             = not achieved\n');
        end

        fprintf('  Minimum separation       = %.6f m\n',minDistanceAll(modeIndex,controller));
        fprintf('  Final separation         = %.6f m\n',finalDistanceAll(modeIndex,controller));
        fprintf('  Cumulative game cost     = %.6e\n',costAll(modeIndex,controller));
        fprintf('  Pursuer control effort   = %.6e\n',pursuerEffortAll(modeIndex,controller));
        fprintf('  Evader control effort    = %.6e\n',evaderEffortAll(modeIndex,controller));
        fprintf('  Pursuer saturation count = %d\n',pursuerSatAll(modeIndex,controller));
        fprintf('  Evader saturation count  = %d\n',evaderSatAll(modeIndex,controller));

    end

end

%% ================================================================
% 11. THESIS-READY FIGURE EXPORT
%% ================================================================

% This section only controls exported figures and thesis summary files.
% It does not change the simulation, controller, fitted-Q recovery, or
% numerical results computed above.

outputFolder = fullfile(pwd,'PEV_Thesis_Figures');

if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

fprintf('\nSaving thesis figures to:\n%s\n\n',outputFolder);

desiredModeNames = ["Straight","limited_Saddle","Saddle"];
displayLabels = ["Straight","Limited Saddle","Full Saddle"];

modeOrder = zeros(1,numel(desiredModeNames));

for modeIndex = 1:numel(desiredModeNames)

    locatedIndex = find(evaderModes == desiredModeNames(modeIndex),1);

    if isempty(locatedIndex)
        error('Evader mode "%s" was not found.',desiredModeNames(modeIndex));
    end

    modeOrder(modeIndex) = locatedIndex;

end

orderedCategories = categorical( ...
    displayLabels, ...
    displayLabels, ...
    'Ordinal',true);

minDistancePlot = minDistanceAll(modeOrder,:);
captureTimePlot = captureTimeAll(modeOrder,:);
pursuerSatPlot  = pursuerSatAll(modeOrder,:);
stepsUsedPlot   = stepsUsedAll(modeOrder,:);

captureMarginMM = 1000*(captureRadius-minDistancePlot);
pursuerSaturationPercent = 100*pursuerSatPlot ./ max(stepsUsedPlot,1);

figureWidth  = 1100;
figureHeight = 720;
exportResolution = 300;
lineWidth = 1.8;

%% ------------------------------------------------
% 11A. Capture margin figure
%% ------------------------------------------------

fig = figure('Color','white','Position',[100 100 figureWidth figureHeight]);
bar(orderedCategories,captureMarginMM);

ylabel('Capture margin [mm]');
title('Capture Margin Across Evader Policies');

legend('Exact LQ','Learned Q','Location','best');

grid on;
set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1.0,'Box','on');
set(fig,'PaperPositionMode','auto');

print(fig, ...
    fullfile(outputFolder,'01_capture_margin_across_evader_modes.png'), ...
    '-dpng',sprintf('-r%d',exportResolution));

close(fig);

%% ------------------------------------------------
% 11B. Capture time figure
%% ------------------------------------------------

fig = figure('Color','white','Position',[100 100 figureWidth figureHeight]);
bar(orderedCategories,captureTimePlot);

ylabel('Capture time [s]');
title('Capture Time Across Evader Policies');

legend('Exact LQ','Learned Q','Location','best');

grid on;
set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1.0,'Box','on');
set(fig,'PaperPositionMode','auto');

print(fig, ...
    fullfile(outputFolder,'02_capture_time_across_evader_modes.png'), ...
    '-dpng',sprintf('-r%d',exportResolution));

close(fig);

%% ------------------------------------------------
% 11C. Saturation percentage figure
%% ------------------------------------------------

fig = figure('Color','white','Position',[100 100 figureWidth figureHeight]);
bar(orderedCategories,pursuerSaturationPercent);

ylabel('Saturated time steps [%]');
title('Pursuer Input Saturation Across Evader Policies');

legend('Exact LQ','Learned Q','Location','best');

grid on;
set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1.0,'Box','on');
set(fig,'PaperPositionMode','auto');

print(fig, ...
    fullfile(outputFolder,'03_pursuer_saturation_percentage.png'), ...
    '-dpng',sprintf('-r%d',exportResolution));

close(fig);

%% ------------------------------------------------
% 11D. Trajectory and separation figures
%% ------------------------------------------------

for orderedModeIndex = 1:numel(desiredModeNames)

    originalModeIndex = modeOrder(orderedModeIndex);
    modeTitle = displayLabels(orderedModeIndex);
    modeFileName = char(desiredModeNames(orderedModeIndex));

    exactControllerIndex = 1;
    learnedControllerIndex = 2;

    exactFinalStep = stepsUsedAll(originalModeIndex,exactControllerIndex);
    learnedFinalStep = stepsUsedAll(originalModeIndex,learnedControllerIndex);

    exactIndices = 1:(exactFinalStep+1);
    learnedIndices = 1:(learnedFinalStep+1);

    %% Trajectory

    fig = figure('Color','white','Position',[100 100 figureWidth figureHeight]);

    plot(Xall(1,exactIndices,originalModeIndex,exactControllerIndex), ...
         Xall(2,exactIndices,originalModeIndex,exactControllerIndex), ...
         '-', ...
         'LineWidth',lineWidth);

    hold on;

    plot(Xall(5,exactIndices,originalModeIndex,exactControllerIndex), ...
         Xall(6,exactIndices,originalModeIndex,exactControllerIndex), ...
         '--', ...
         'LineWidth',lineWidth);

    markerSpacing = max(1,floor(numel(learnedIndices)/14));
    markerIndices = 1:markerSpacing:numel(learnedIndices);

    plot(Xall(1,learnedIndices(markerIndices),originalModeIndex,learnedControllerIndex), ...
         Xall(2,learnedIndices(markerIndices),originalModeIndex,learnedControllerIndex), ...
         'o', ...
         'LineStyle','none', ...
         'MarkerSize',5, ...
         'LineWidth',1.1);

    plot(Xall(5,learnedIndices(markerIndices),originalModeIndex,learnedControllerIndex), ...
         Xall(6,learnedIndices(markerIndices),originalModeIndex,learnedControllerIndex), ...
         's', ...
         'LineStyle','none', ...
         'MarkerSize',5, ...
         'LineWidth',1.1);

    plot(x0(1),x0(2),'o','MarkerSize',9,'LineWidth',1.6);
    plot(x0(5),x0(6),'s','MarkerSize',9,'LineWidth',1.6);

    capturePX = Xall(1,exactIndices(end),originalModeIndex,exactControllerIndex);
    capturePY = Xall(2,exactIndices(end),originalModeIndex,exactControllerIndex);
    captureEX = Xall(5,exactIndices(end),originalModeIndex,exactControllerIndex);
    captureEY = Xall(6,exactIndices(end),originalModeIndex,exactControllerIndex);

    plot(capturePX,capturePY,'x','MarkerSize',11,'LineWidth',2.0);
    plot(captureEX,captureEY,'x','MarkerSize',11,'LineWidth',2.0);

    circleAngle = linspace(0,2*pi,300);

    plot(captureEX + captureRadius*cos(circleAngle), ...
         captureEY + captureRadius*sin(circleAngle), ...
         ':', ...
         'LineWidth',1.3);

    xlabel('x position [m]');
    ylabel('y position [m]');

    title(sprintf('Pursuit-Evasion Trajectories: %s Evader',modeTitle));

    legend( ...
        'Exact pursuer', ...
        'Exact evader', ...
        'Learned pursuer samples', ...
        'Learned evader samples', ...
        'Pursuer initial position', ...
        'Evader initial position', ...
        'Pursuer capture position', ...
        'Evader capture position', ...
        'Capture boundary', ...
        'Location','best');

    grid on;
    axis equal;
    set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1.0,'Box','on');
    set(fig,'PaperPositionMode','auto');

    print(fig, ...
        fullfile(outputFolder,sprintf('%02d_trajectory_%s_evader.png', ...
        2*orderedModeIndex+2,modeFileName)), ...
        '-dpng',sprintf('-r%d',exportResolution));

    close(fig);

    %% Separation

    exactTime = 0:Ts:(exactFinalStep*Ts);
    learnedTime = 0:Ts:(learnedFinalStep*Ts);

    exactDistance = sqrt( ...
        (Xall(1,exactIndices,originalModeIndex,exactControllerIndex) - ...
         Xall(5,exactIndices,originalModeIndex,exactControllerIndex)).^2 + ...
        (Xall(2,exactIndices,originalModeIndex,exactControllerIndex) - ...
         Xall(6,exactIndices,originalModeIndex,exactControllerIndex)).^2);

    learnedDistance = sqrt( ...
        (Xall(1,learnedIndices,originalModeIndex,learnedControllerIndex) - ...
         Xall(5,learnedIndices,originalModeIndex,learnedControllerIndex)).^2 + ...
        (Xall(2,learnedIndices,originalModeIndex,learnedControllerIndex) - ...
         Xall(6,learnedIndices,originalModeIndex,learnedControllerIndex)).^2);

    fig = figure('Color','white','Position',[100 100 figureWidth figureHeight]);

    plot(exactTime,exactDistance,'-','LineWidth',lineWidth);

    hold on;

    markerSpacing = max(1,floor(numel(learnedTime)/15));
    markerIndices = 1:markerSpacing:numel(learnedTime);

    plot(learnedTime(markerIndices),learnedDistance(markerIndices), ...
         'o', ...
         'LineStyle','none', ...
         'MarkerSize',5, ...
         'LineWidth',1.1);

    yline(captureRadius,':','Capture radius','LineWidth',1.5);

    xlabel('Time [s]');
    ylabel('Separation distance [m]');

    title(sprintf('Pursuer-Evader Separation: %s Evader',modeTitle));

    legend('Exact LQ','Learned Q samples','Capture radius','Location','best');

    grid on;
    set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1.0,'Box','on');
    set(fig,'PaperPositionMode','auto');

    print(fig, ...
        fullfile(outputFolder,sprintf('%02d_separation_%s_evader.png', ...
        2*orderedModeIndex+3,modeFileName)), ...
        '-dpng',sprintf('-r%d',exportResolution));

    close(fig);

end

%% ------------------------------------------------
% 11E. Bellman residual
%% ------------------------------------------------

fig = figure('Color','white','Position',[100 100 figureWidth figureHeight]);

semilogy(1:N,bellmanResidual,'LineWidth',lineWidth);

xlabel('Finite-horizon stage');
ylabel('Relative Bellman residual');

title('Stagewise Fitted-Q Bellman Residual');

grid on;
set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1.0,'Box','on');
set(fig,'PaperPositionMode','auto');

print(fig, ...
    fullfile(outputFolder,'10_fitted_q_bellman_residual.png'), ...
    '-dpng',sprintf('-r%d',exportResolution));

close(fig);

%% ------------------------------------------------
% 11F. Gain-recovery error
%% ------------------------------------------------

gainErrorP = zeros(N,1);
gainErrorE = zeros(N,1);

for k = 1:N

    gainErrorP(k) = ...
        norm(KpLearn(:,:,k)-KpExact(:,:,k),'fro') / ...
        max(norm(KpExact(:,:,k),'fro'),eps);

    gainErrorE(k) = ...
        norm(KeLearn(:,:,k)-KeExact(:,:,k),'fro') / ...
        max(norm(KeExact(:,:,k),'fro'),eps);

end

fig = figure('Color','white','Position',[100 100 figureWidth figureHeight]);

semilogy(1:N,gainErrorP,'-','LineWidth',1.5);

hold on;

semilogy(1:N,gainErrorE,'--','LineWidth',1.5);

xlabel('Finite-horizon stage');
ylabel('Relative gain error');

title('Fitted-Q Recovery Error for the Kinematic-Bicycle LQ Game');

legend('Pursuer gain error','Evader gain error','Location','best');

grid on;
set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1.0,'Box','on');
set(fig,'PaperPositionMode','auto');

print(fig, ...
    fullfile(outputFolder,'11_fitted_q_gain_recovery_error.png'), ...
    '-dpng',sprintf('-r%d',exportResolution));

close(fig);

%% ------------------------------------------------
% 11G. Save numerical summaries
%% ------------------------------------------------

summaryTable = table( ...
    displayLabels', ...
    captureTimePlot(:,1), ...
    captureTimePlot(:,2), ...
    minDistancePlot(:,1), ...
    minDistancePlot(:,2), ...
    captureMarginMM(:,1), ...
    captureMarginMM(:,2), ...
    pursuerSaturationPercent(:,1), ...
    pursuerSaturationPercent(:,2), ...
    'VariableNames',{ ...
        'EvaderMode', ...
        'ExactCaptureTime_s', ...
        'LearnedCaptureTime_s', ...
        'ExactMinimumSeparation_m', ...
        'LearnedMinimumSeparation_m', ...
        'ExactCaptureMargin_mm', ...
        'LearnedCaptureMargin_mm', ...
        'ExactPursuerSaturation_percent', ...
        'LearnedPursuerSaturation_percent'});

writetable(summaryTable, ...
    fullfile(outputFolder,'nominal_evader_mode_summary_wide.csv'));

summaryTableCompact = table( ...
    displayLabels', ...
    captureTimePlot(:,1), ...
    captureTimePlot(:,2), ...
    minDistancePlot(:,1), ...
    minDistancePlot(:,2), ...
    captureMarginMM(:,1), ...
    captureMarginMM(:,2), ...
    pursuerSaturationPercent(:,1), ...
    pursuerSaturationPercent(:,2), ...
    'VariableNames',{ ...
        'Mode', ...
        't_LQ_s', ...
        't_Q_s', ...
        'dmin_LQ_m', ...
        'dmin_Q_m', ...
        'margin_LQ_mm', ...
        'margin_Q_mm', ...
        'sat_LQ_pct', ...
        'sat_Q_pct'});

writetable(summaryTableCompact, ...
    fullfile(outputFolder,'nominal_evader_mode_summary_compact.csv'));

longMode = strings(numModes*numControllers,1);
longController = strings(numModes*numControllers,1);
longCaptured = false(numModes*numControllers,1);
longCaptureTime = nan(numModes*numControllers,1);
longMinSeparation = nan(numModes*numControllers,1);
longCaptureMargin = nan(numModes*numControllers,1);
longPursuerSatPercent = nan(numModes*numControllers,1);
longPursuerEffort = nan(numModes*numControllers,1);
longEvaderEffort = nan(numModes*numControllers,1);
longCost = nan(numModes*numControllers,1);

rowCounter = 0;

for orderedModeIndex = 1:numModes

    originalModeIndex = modeOrder(orderedModeIndex);

    for controller = 1:numControllers

        rowCounter = rowCounter + 1;

        longMode(rowCounter) = displayLabels(orderedModeIndex);
        longController(rowCounter) = controllerNames(controller);
        longCaptured(rowCounter) = capturedAll(originalModeIndex,controller);
        longCaptureTime(rowCounter) = captureTimeAll(originalModeIndex,controller);
        longMinSeparation(rowCounter) = minDistanceAll(originalModeIndex,controller);
        longCaptureMargin(rowCounter) = 1000*(captureRadius - minDistanceAll(originalModeIndex,controller));
        longPursuerSatPercent(rowCounter) = ...
            100*pursuerSatAll(originalModeIndex,controller) / ...
            max(stepsUsedAll(originalModeIndex,controller),1);
        longPursuerEffort(rowCounter) = pursuerEffortAll(originalModeIndex,controller);
        longEvaderEffort(rowCounter) = evaderEffortAll(originalModeIndex,controller);
        longCost(rowCounter) = costAll(originalModeIndex,controller);

    end

end

summaryTableLong = table( ...
    longMode, ...
    longController, ...
    longCaptured, ...
    longCaptureTime, ...
    longMinSeparation, ...
    longCaptureMargin, ...
    longPursuerSatPercent, ...
    longPursuerEffort, ...
    longEvaderEffort, ...
    longCost, ...
    'VariableNames',{ ...
        'EvaderMode', ...
        'Controller', ...
        'Captured', ...
        'CaptureTime_s', ...
        'MinimumSeparation_m', ...
        'CaptureMargin_mm', ...
        'PursuerSaturation_percent', ...
        'PursuerEffort', ...
        'EvaderEffort', ...
        'GameCost'});

writetable(summaryTableLong, ...
    fullfile(outputFolder,'nominal_evader_mode_summary_long.csv'));

%% ------------------------------------------------
% 11H. Write a thesis-readable LaTeX table
%% ------------------------------------------------

latexFile = fullfile(outputFolder,'nominal_evader_mode_summary_vertical_table.tex');
fid = fopen(latexFile,'w');

if fid == -1
    warning('Could not create LaTeX summary table file.');
else

    fprintf(fid,'%% Auto-generated by MATLAB script.\n');
    fprintf(fid,'%% Paste into the thesis and adjust caption/label if needed.\n\n');

    fprintf(fid,'\\begin{table}[H]\n');
    fprintf(fid,'\\centering\n');
    fprintf(fid,'\\caption{Nominal evader-mode summary for the Riccati and fitted-Q policies.}\n');
    fprintf(fid,'\\label{tab:nominal_evader_mode_summary}\n');
    fprintf(fid,'\\renewcommand{\\arraystretch}{1.18}\n');
    fprintf(fid,'\\setlength{\\tabcolsep}{4pt}\n');
    fprintf(fid,'\\small\n');
    fprintf(fid,'\\begin{tabular}{p{0.22\\textwidth} p{0.16\\textwidth} c c c c}\n');
    fprintf(fid,'\\toprule\n');
    fprintf(fid,'\\textbf{Evader mode} & \\textbf{Controller} & \\textbf{Captured} & \\textbf{$t_c$ [s]} & \\textbf{$d_{\\min}$ [m]} & \\textbf{Sat. [\\%%]} \\\\\n');
    fprintf(fid,'\\midrule\n');

    for r = 1:height(summaryTableLong)

        if summaryTableLong.Captured(r)
            capturedText = 'Yes';
        else
            capturedText = 'No';
        end

        fprintf(fid,'%s & %s & %s & %.2f & %.4f & %.2f \\\\\n', ...
            char(summaryTableLong.EvaderMode(r)), ...
            char(summaryTableLong.Controller(r)), ...
            capturedText, ...
            summaryTableLong.CaptureTime_s(r), ...
            summaryTableLong.MinimumSeparation_m(r), ...
            summaryTableLong.PursuerSaturation_percent(r));

        if mod(r,2) == 0 && r < height(summaryTableLong)
            fprintf(fid,'\\addlinespace[0.25em]\n');
        end

    end

    fprintf(fid,'\\bottomrule\n');
    fprintf(fid,'\\end{tabular}\n');
    fprintf(fid,'\\end{table}\n');

    fclose(fid);

end

%% ------------------------------------------------
% 11I. Save a reproducibility manifest
%% ------------------------------------------------

manifestFile = fullfile(outputFolder,'run_manifest.txt');
fid = fopen(manifestFile,'w');

if fid == -1
    warning('Could not create run manifest.');
else

    fprintf(fid,'PEV thesis figure generation manifest\n');
    fprintf(fid,'=====================================\n\n');

    fprintf(fid,'Random generator: rng(11,''twister'')\n');
    fprintf(fid,'Sampling time Ts: %.12g s\n',Ts);
    fprintf(fid,'Horizon stages N: %d\n',N);
    fprintf(fid,'Horizon time: %.12g s\n',N*Ts);
    fprintf(fid,'Wheelbase Lp: %.12g m\n',Lp);
    fprintf(fid,'Wheelbase Le: %.12g m\n',Le);
    fprintf(fid,'Capture radius: %.12g m\n',captureRadius);
    fprintf(fid,'Initial condition:\n');

    for i = 1:numel(x0)
        fprintf(fid,'  x0(%d) = %.16g\n',i,x0(i));
    end

    fprintf(fid,'\nSelected game parameters:\n');
    fprintf(fid,'  Q scale: %.12g\n',selectedQScale);
    fprintf(fid,'  Terminal multiplier: %.12g\n',selectedTerminal);
    fprintf(fid,'  Rp scale: %.12g\n',selectedRpScale);
    fprintf(fid,'  Re scale: %.12g\n',selectedReScale);
    fprintf(fid,'  Minimum pursuer curvature: %.16e\n',bestMinPcurv);
    fprintf(fid,'  Maximum evader Schur eigenvalue: %.16e\n',bestMaxEschur);
    fprintf(fid,'  Maximum saddle condition: %.16e\n',bestMaxCond);

    fprintf(fid,'\nFitted-Q diagnostics:\n');
    fprintf(fid,'  Unique Q coefficients: %d\n',nFeatures);
    fprintf(fid,'  Samples per stage: %d\n',samplesPerStage);
    fprintf(fid,'  Relative pursuer gain error: %.16e\n',pursuerGainError);
    fprintf(fid,'  Relative evader gain error: %.16e\n',evaderGainError);
    fprintf(fid,'  Maximum Bellman residual: %.16e\n',max(bellmanResidual));
    fprintf(fid,'  Minimum feature rank: %d of %d\n',min(featureRank),nFeatures);

    fprintf(fid,'\nOutput files:\n');
    fprintf(fid,'  01_capture_margin_across_evader_modes.png\n');
    fprintf(fid,'  02_capture_time_across_evader_modes.png\n');
    fprintf(fid,'  03_pursuer_saturation_percentage.png\n');
    fprintf(fid,'  04_trajectory_straight_evader.png\n');
    fprintf(fid,'  05_separation_straight_evader.png\n');
    fprintf(fid,'  06_trajectory_limited_saddle_evader.png\n');
    fprintf(fid,'  07_separation_limited_saddle_evader.png\n');
    fprintf(fid,'  08_trajectory_saddle_evader.png\n');
    fprintf(fid,'  09_separation_saddle_evader.png\n');
    fprintf(fid,'  10_fitted_q_bellman_residual.png\n');
    fprintf(fid,'  11_fitted_q_gain_recovery_error.png\n');
    fprintf(fid,'  nominal_evader_mode_summary_wide.csv\n');
    fprintf(fid,'  nominal_evader_mode_summary_compact.csv\n');
    fprintf(fid,'  nominal_evader_mode_summary_long.csv\n');
    fprintf(fid,'  nominal_evader_mode_summary_vertical_table.tex\n');

    fclose(fid);

end

fprintf('\n================================================\n');
fprintf('Simulation and figure export complete.\n');
fprintf('PNG, CSV, LaTeX table, and manifest files saved in:\n%s\n',outputFolder);
fprintf('================================================\n');
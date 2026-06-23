%% ================================================================
% PEV Fitted-Q Recovery with Comments.m
%
% Chapter 5: Stage-Wise Fitted Quadratic Q-function Recovery
%
% This script reconstructs the finite-horizon Riccati saddle policy from
% sampled nominal transitions and exports only Chapter 5 diagnostics:
%   Fig. 5.1: Bellman residuals
%   Fig. 5.2: stage-wise gain-recovery error
%   Fig. 5.3: full-saddle separation histories
%
% Chapter 4 figures are intentionally NOT generated here.
%% ================================================================

clear;
clc;
close all;

% The fixed seed makes the sampled-transition experiment repeatable and keeps the Chapter 5 recovery diagnostics reproducible.
rng(11,'twister');

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
% Chapter 5 deliberately retains the same two-vehicle state used in Chapters 3 and 4; the comparison changes the route to the policy, not the task being solved.
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

% These matrices reproduce the Chapter 4 nominal linear game so that the Riccati solution remains the analytical reference for fitted-Q recovery.
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

% The Chapter 4 saddle construction is repeated here only to generate the exact reference gains; it is not treated as an additional competing controller.
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

% The learning experiment keeps the nominal finite-horizon game fixed and asks whether sampled transitions can reconstruct its stage-dependent saddle policy.
%% ================================================================
% 8. FITTED QUADRATIC Q-LEARNING RECOVERY
%% ================================================================

% A symmetric quadratic form in the 12 augmented state-action variables has 78 unique coefficients, which defines the regression problem at every stage.
nFeatures = nz*(nz+1)/2;

% The sample batch is intentionally larger than the quadratic basis so the code can check that the fitted feature matrix has full rank.
samplesPerStage = max(850,7*nFeatures);
lambda = 1e-10;

Plearn  = zeros(nx,nx,N+1);
KpLearn = zeros(nup,nx,N);
KeLearn = zeros(nue,nx,N);

bellmanResidual = zeros(N,1);
featureRank = zeros(N,1);
featureCondition = zeros(N,1);

% As in the Riccati recursion, fitted-Q works backward from the terminal objective because each stage depends on the value remaining after the successor state.
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

        % The nominal linear model generates transition samples, while the regression below uses the sampled tuples to recover the quadratic action-value representation.
        xnext = A*x + Bp*up + Be*ue;

        target = ...
            x'*Q*x + ...
            up'*Rp*up - ...
            ue'*Re*ue + ...
            xnext'*Pnext*xnext;

        % The feature vector contains every unique quadratic term needed to represent the joint state--pursuer-action--evader-action Q-function.
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

    % Rank and conditioning are recorded as implementation diagnostics: full rank is required for identifiability, while conditioning indicates numerical sensitivity of the least-squares solve.
    featureRank(k) = rank(Phi);
    featureCondition(k) = cond(Phi'*Phi + lambda*eye(nFeatures));

    if featureRank(k) < nFeatures
        error('Feature matrix rank deficient at stage %d.',k);
    end

    theta = (Phi'*Phi + lambda*eye(nFeatures))\(Phi'*Y);

    % Rebuild the symmetric quadratic matrix so its state and joint-action blocks can be interpreted exactly as in the Chapter 5 derivation.
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

    % These curvature checks confirm that the learned quadratic form still defines a local zero-sum saddle before gains are extracted.
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

    % This in-sample residual is a recovery diagnostic for the matched nominal experiment; it is not presented as a held-out generalization result.
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

% Both gain sequences are now evaluated on the same constrained nonlinear rollout. This is a nominal policy-recovery check, not a separate model-based-versus-RL competition.
%% ================================================================
% 9. FULL-SADDLE CLOSED-LOOP RECOVERY CHECK
%% ================================================================

controllerNames = ["Riccati saddle policy","Fitted-Q recovered policy"];
numControllers = 2;

captured = false(numControllers,1);
captureStep = N*ones(numControllers,1);
captureTime = N*Ts*ones(numControllers,1);
minDistance = nan(numControllers,1);
Xcheck = zeros(nx,N+1,numControllers);

bike = @(z,u,L) [ ...
    z(4)*cos(z(3)); ...
    z(4)*sin(z(3)); ...
    z(4)/L*tan(u(2)); ...
    u(1)];

for controller = 1:numControllers

    X = zeros(nx,N+1);
    X(:,1) = x0;

    for k = 1:N
        currentState = X(:,k);
        errorState = currentState - xbar;

        % The only switch in this loop is the gain sequence. Initial state, constraints, nonlinear propagation, and capture condition remain identical.
        if controller == 1
            KpCurrent = KpExact(:,:,k);
            KeCurrent = KeExact(:,:,k);
        else
            KpCurrent = KpLearn(:,:,k);
            KeCurrent = KeLearn(:,:,k);
        end

        upRaw = -KpCurrent*errorState;
        ueRaw = -KeCurrent*errorState;

        up = [min(max(upRaw(1),-amaxP),amaxP); ...
              min(max(upRaw(2),-deltaMaxP),deltaMaxP)];
        ue = [min(max(ueRaw(1),-amaxE),amaxE); ...
              min(max(ueRaw(2),-deltaMaxE),deltaMaxE)];

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

        X(:,k+1) = [zPnext;zEnext];

        separation = norm(zPnext(1:2)-zEnext(1:2));
        if separation <= captureRadius
            captured(controller) = true;
            captureStep(controller) = k;
            captureTime(controller) = k*Ts;
            break;
        end
    end

    used = 1:(captureStep(controller)+1);
    d = sqrt((X(1,used)-X(5,used)).^2 + (X(2,used)-X(6,used)).^2);
    minDistance(controller) = min(d);
    Xcheck(:,:,controller) = X;
end

% Stage-wise errors show whether recovery changes with time-to-go, while the aggregate errors below summarize agreement over the whole horizon.
gainErrorP = zeros(N,1);
gainErrorE = zeros(N,1);
for k = 1:N
    gainErrorP(k) = norm(KpLearn(:,:,k)-KpExact(:,:,k),'fro') / max(norm(KpExact(:,:,k),'fro'),eps);
    gainErrorE(k) = norm(KeLearn(:,:,k)-KeExact(:,:,k),'fro') / max(norm(KeExact(:,:,k),'fro'),eps);
end

fprintf('================================================\n');
fprintf('CHAPTER 5: FITTED-Q RECOVERY RESULTS\n');
fprintf('================================================\n');
fprintf('Relative pursuer gain error = %.6e\n',pursuerGainError);
fprintf('Relative evader gain error  = %.6e\n',evaderGainError);
fprintf('Maximum Bellman residual    = %.6e\n',max(bellmanResidual));
fprintf('Full-saddle Riccati capture = %.2f s, d_min = %.4f m\n',captureTime(1),minDistance(1));
fprintf('Full-saddle fitted-Q capture= %.2f s, d_min = %.4f m\n',captureTime(2),minDistance(2));

% These diagnostics provide the evidence used in Chapter 5: regression consistency, gain recovery, and common-rollout agreement with the Chapter 4 reference.
%% ================================================================
% 10. CHAPTER 5 FIGURE EXPORT
%% ================================================================

scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    scriptDir = pwd;
else
    scriptDir = fileparts(scriptPath);
end
outputFolder = fullfile(scriptDir,'PEV_Ch5_FittedQ_Recovery_Output_v11');
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

figSize = [100 100 1050 700];
exportResolution = 300;

% Figure 5.1: Bellman residuals.
fig = figure('Color','w','InvertHardcopy','off','Position',figSize);
hResidual = semilogy(1:N,bellmanResidual,'LineWidth',1.8);
xlabel('Finite-horizon stage','Color','k');
ylabel('Relative Bellman residual','Color','k');
title('Bellman Residuals for Backward Fitted-Q Recovery','Color','k');
grid on;
box on;
set(gca,'Color','w','XColor','k','YColor','k', ...
    'FontName','Times New Roman','FontSize',12,'LineWidth',1.0);
lgd = legend(hResidual,'Relative Bellman residual','Location','northeast');
set(lgd,'Color','w','TextColor','k','EdgeColor','k', ...
    'FontName','Times New Roman','FontSize',10.5);
set(fig,'PaperPositionMode','auto');
print(fig,fullfile(outputFolder,'fig5_1_bellman_residuals.png'), ...
    '-dpng',sprintf('-r%d',exportResolution));
close(fig);

% Figure 5.2: gain recovery errors.
fig = figure('Color','w','InvertHardcopy','off','Position',figSize);
hP = semilogy(1:N,gainErrorP,'-','LineWidth',1.8);
hold on;
hE = semilogy(1:N,gainErrorE,'--','LineWidth',1.8);
xlabel('Finite-horizon stage','Color','k');
ylabel('Relative gain error','Color','k');
title('Relative Gain Error Between Fitted-Q and Riccati Saddle Gains','Color','k');
grid on;
box on;
set(gca,'Color','w','XColor','k','YColor','k', ...
    'FontName','Times New Roman','FontSize',12,'LineWidth',1.0);
lgd = legend([hP hE],{'Pursuer gain error','Evader gain error'}, ...
    'Location','northeast');
set(lgd,'Color','w','TextColor','k','EdgeColor','k', ...
    'FontName','Times New Roman','FontSize',10.5);
set(fig,'PaperPositionMode','auto');
print(fig,fullfile(outputFolder,'fig5_2_gain_recovery_error.png'), ...
    '-dpng',sprintf('-r%d',exportResolution));
close(fig);

% Figure 5.3: source-of-record full-saddle separation history.
lastExact = captureStep(1);
lastLearn = captureStep(2);
idxExact = 1:(lastExact+1);
idxLearn = 1:(lastLearn+1);
tExact = (0:lastExact)*Ts;
tLearn = (0:lastLearn)*Ts;

dExact = sqrt((Xcheck(1,idxExact,1)-Xcheck(5,idxExact,1)).^2 + ...
              (Xcheck(2,idxExact,1)-Xcheck(6,idxExact,1)).^2);
dLearn = sqrt((Xcheck(1,idxLearn,2)-Xcheck(5,idxLearn,2)).^2 + ...
              (Xcheck(2,idxLearn,2)-Xcheck(6,idxLearn,2)).^2);

fig = figure('Color','w','InvertHardcopy','off','Position',figSize);
hRiccati = plot(tExact,dExact,'-','LineWidth',1.9);
hold on;
markerSpacing = max(1,floor(numel(tLearn)/18));
markerIdx = unique([1:markerSpacing:numel(tLearn), numel(tLearn)]);
hFittedQ = plot(tLearn(markerIdx),dLearn(markerIdx),'o', ...
    'LineStyle','none','MarkerSize',5,'LineWidth',1.1);
hCap = yline(captureRadius,':','LineWidth',1.5);
hCapture = plot(tExact(end),dExact(end),'x','MarkerSize',10,'LineWidth',1.8);
xlabel('Time [s]','Color','k');
ylabel('Separation distance [m]','Color','k');
title('Separation Histories Under the Riccati and Fitted-Q Recovered Policies', ...
    'Color','k');
grid on;
box on;
set(gca,'Color','w','XColor','k','YColor','k', ...
    'FontName','Times New Roman','FontSize',12,'LineWidth',1.0);
lgd = legend([hRiccati hFittedQ hCap hCapture], ...
    {'Riccati saddle-policy separation', ...
     'Fitted-Q recovered-policy samples', ...
     'Capture radius', ...
     'First sampled capture'}, ...
    'Location','best');
set(lgd,'Color','w','TextColor','k','EdgeColor','k', ...
    'FontName','Times New Roman','FontSize',10.5);
set(fig,'PaperPositionMode','auto');
print(fig,fullfile(outputFolder,'fig5_3_full_saddle_separation.png'), ...
    '-dpng',sprintf('-r%d',exportResolution));
close(fig);

% The CSV is retained so the numerical recovery evidence behind the Chapter 5 figures and table can be checked directly.
diagnosticTable = table( ...
    nFeatures,samplesPerStage,min(featureRank),max(featureCondition), ...
    pursuerGainError,evaderGainError,max(bellmanResidual), ...
    captureTime(1),captureTime(2),minDistance(1),minDistance(2), ...
    'VariableNames',{'UniqueQuadraticCoefficients','SamplesPerStage','MinimumFeatureRank', ...
    'MaximumFeatureCondition','RelativePursuerGainError','RelativeEvaderGainError', ...
    'MaximumBellmanResidual','RiccatiCaptureTime_s','FittedQCaptureTime_s', ...
    'RiccatiMinimumSeparation_m','FittedQMinimumSeparation_m'});
writetable(diagnosticTable,fullfile(outputFolder,'chapter5_fittedq_diagnostics.csv'));

fid = fopen(fullfile(outputFolder,'run_manifest.txt'),'w');
if fid ~= -1
    fprintf(fid,'Chapter 5 fitted-Q recovery run manifest\n');
    fprintf(fid,'Seed: rng(11,''twister'')\n');
    fprintf(fid,'Fitted-Q function class: exact 12-by-12 symmetric quadratic form (78 coefficients).\n');
    fprintf(fid,'Training transitions are sampled from the nominal discrete linear game.\n');
    fprintf(fid,'Figure 5.3 is a recovery check, not an independent RL-versus-model-based benchmark.\n');
    fprintf(fid,'Riccati full-saddle first sampled capture: %.6f s.\n',captureTime(1));
    fprintf(fid,'Fitted-Q full-saddle first sampled capture: %.6f s.\n',captureTime(2));
    fclose(fid);
end

fprintf('\nChapter 5 v11 fitted-Q figures and CSV diagnostics saved to:\n%s\n',outputFolder);

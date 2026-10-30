function R = drone_swarm_rescue_v2(userOpts)
% DRONE_SWARM_RESCUE_V2  Search-and-rescue UAV swarm simulation (version 2).
%
% What is modelled (see the report for the theory behind each block):
%   * Terrain, buildings (cylinders), a ground control station (GCS) and charging pads
%   * A missing person whose location is uncertain: Bayesian prior around a
%     Last Known Position (LKP)
%   * Downward camera: footprint grows with altitude, detection probability falls
%   * Bayesian search: the belief map is updated after every "no detection"
%   * Four search strategies: 'bayes' | 'lawnmower' | 'spiral' | 'random'
%   * Radio links: log-distance path loss + shadowing + building blockage
%     -> SNR -> packet error rate (no fixed "comm radius")
%   * Multi-hop reporting: the alert travels one hop per time step; the finder acts
%     as a data mule if it is cut off; informed drones build a relay chain to the GCS
%   * Energy: speed-dependent power draw, energy-aware return-to-charge
%   * Random hardware failures (exponential lifetime), crash-site inspection and
%     hand-over of the failed drone's unfinished search path
%   * Graph metrics: GCS-reachability, algebraic connectivity (lambda_2)
%
% Usage
%   R = drone_swarm_rescue_v2();                                   % dialog + live dashboard
%   R = drone_swarm_rescue_v2(struct('searchPattern','lawnmower','seed',3));
%   R = drone_swarm_rescue_v2(struct('visualize',false));          % headless (Monte Carlo)
%
% Works in MATLAB R2018b+ (no toolboxes needed) and GNU Octave 7+.
% Authors: Niyathi Vedagiri (23BRS1132), Anna Sunny George (23BRS1321)

if nargin < 1, userOpts = struct(); end
P = defaultParams();
P.useDialog = (nargin == 0);
P = mergeStruct(P, userOpts);
if P.useDialog
    P = paramDialog(P);
    if isempty(P), R = []; disp('Simulation cancelled.'); return; end
end
P.numDrones = max(2, round(P.numDrones));
setSeed(P.seed);
M = modes();

%% ------------------------------------------------------------------ world
W   = buildWorld(P);
N   = P.numDrones;
dt  = P.dt;
nSteps = round(P.simTime / dt);
rfp = P.searchAlt * tand(P.camFOVdeg/2);          % footprint radius at search altitude
P.rangeDD  = linkRange(P, 0);                      % drone-drone nominal range (PER ~ 5 %)
P.rangeGCS = linkRange(P, P.gcsGain);              % drone-GCS nominal range
P.relaySpacing = 0.8 * P.rangeDD;

%% ------------------------------------------------------------------ swarm state
S.pos = zeros(N,3);
ang = linspace(0, pi/2, N)';  rad = 10 + 12*rand(N,1);
S.pos(:,1) = W.base(1) + rad.*cos(ang);
S.pos(:,2) = W.base(2) + rad.*sin(ang);
S.pos(:,3) = groundAt(W, S.pos(:,1), S.pos(:,2)) + 0.3;
S.vel      = zeros(N,3);
S.Emax     = P.batteryJ;
S.E        = P.batteryJ * (P.socInitRange(1) + diff(P.socInitRange)*rand(N,1));
S.mode     = M.SEARCH * ones(N,1);
S.alive    = true(N,1);
S.informed = false(N,1);
S.goal     = S.pos;
S.goalTimer= zeros(N,1);
S.charger  = ones(N,1);
S.inspectPos = nan(N,2);
S.slotPos  = nan(N,2);
S.orbitAng = 2*pi*rand(N,1);
S.wp       = planSearch(P, W, rfp);
S.wpIdx    = ones(N,1);
S.q        = W.prior;                 % un-normalised "not yet detected" mass
S.covered  = false(size(W.prior));
S.detected = false;  S.finder = 0;  S.gcsInformed = false;
S.events   = {};

T_detect = NaN; T_report = NaN; T_all = NaN; covAtDetect = NaN;
nextAssign = 0;

%% ------------------------------------------------------------------ history
H.t        = (1:nSteps)*dt;
H.pos      = zeros(N,3,nSteps);
H.mode     = zeros(N,nSteps);
H.soc      = zeros(N,nSteps);
H.inf      = false(N,nSteps);
H.cov      = zeros(1,nSteps);
H.POS      = zeros(1,nSteps);
H.baseConn = zeros(1,nSteps);
H.informed = zeros(1,nSteps);
H.lam2     = zeros(1,nSteps);
H.nLinks   = zeros(1,nSteps);
H.alive    = zeros(1,nSteps);

if P.verbose
    fprintf('\n=== UAV rescue swarm v2 | %d drones | strategy: %s | seed %d ===\n', N, upper(P.searchPattern), P.seed);
    fprintf('Footprint radius %.1f m | drone-drone range %.0f m | drone-GCS range %.0f m\n', rfp, P.rangeDD, P.rangeGCS);
end

V = [];
if P.visualize, V = initVis(P, W, S, rfp); end

%% ------------------------------------------------------------------ main loop
A = false(N+1); Am = A; Qlink = zeros(N+1); livePath = zeros(0,3);
shadow = P.shadowSigma*randn(N+1); shadow = triu(shadow,1); shadow = shadow + shadow';
lastStep = nSteps;
for k = 1:nSteps
    t = k*dt;

    % --- 1. random hardware failures (exponential lifetime, rate 1/MTBF)
    pFail = 1 - exp(-dt / P.mtbf);
    for i = find(S.alive & S.mode ~= M.CHARGE)'
        if rand < pFail
            [S, msg] = failDrone(S, P, M, i, 'hardware fault');
            S = logEv(S, t, msg, P);
        end
    end

    % --- 2. energy: drain, depletion, energy-aware return, charging
    for i = find(S.alive)'
        if S.mode(i) == M.CHARGE
            S.E(i) = min(S.Emax, S.E(i) + P.chargeW*dt);
            if S.E(i) >= 0.97*S.Emax
                if S.detected && S.informed(i)
                    S.mode(i) = M.OVERWATCH; nextAssign = 0;
                else
                    S.mode(i) = M.SEARCH; S.goalTimer(i) = 0;
                end
                S = logEv(S, t, sprintf('D%d recharged, back in the air', i), P);
            end
            continue
        end
        S.E(i) = S.E(i) - powerDraw(P, S.vel(i,:))*dt;
        if S.E(i) <= 0
            S.E(i) = 0;
            [S, msg] = failDrone(S, P, M, i, 'battery depleted');
            S = logEv(S, t, msg, P);
            continue
        end
        if S.mode(i) ~= M.RETURN
            [dc, ci] = nearestCharger(W, S.pos(i,1:2));
            Ereturn = powerDraw(P, [P.vmax 0 0]) * dc / P.vmax;
            if S.E(i) < P.returnMargin*Ereturn + P.reserveFrac*S.Emax
                S.mode(i) = M.RETURN; S.charger(i) = ci;
                S = logEv(S, t, sprintf('D%d low energy (%.0f%%) -> returning to pad %d', i, 100*S.E(i)/S.Emax, ci), P);
            end
        end
    end

    % --- 3. mission logic after a detection (relay chain / data mule / loiter)
    if S.detected && t >= nextAssign
        S = assignRoles(S, P, W, M);
        nextAssign = t + P.assignEvery;
    end

    % --- 4. guidance: goal for every drone, then dynamics
    [S, ev] = updateGoals(S, P, W, M, rfp, dt);
    for e = 1:numel(ev), S = logEv(S, t, ev{e}, P); end
    S = stepDynamics(S, P, W, M, dt);

    % --- 5. sensing: footprint coverage, Bayesian update, detection
    for i = find(S.alive & S.mode ~= M.CHARGE)'
        agl = S.pos(i,3) - groundAt(W, S.pos(i,1), S.pos(i,2));
        if agl < P.minSenseAGL, continue; end
        r  = agl * tand(P.camFOVdeg/2);
        pdLook = P.pd0 * exp(-agl / P.pdScale);
        pStep  = 1 - (1 - pdLook)^(dt / P.lookInterval);
        inFp = (W.CX - S.pos(i,1)).^2 + (W.CY - S.pos(i,2)).^2 <= r^2;
        S.covered(inFp) = true;
        if ~S.detected
            S.q(inFp) = S.q(inFp) * (1 - pStep);
            if norm(W.target(1:2) - S.pos(i,1:2)) <= r && rand < pStep
                S.detected = true; S.finder = i; S.informed(i) = true;
                T_detect = t; covAtDetect = 100*mean(S.covered(:));
                nextAssign = 0;
                S = logEv(S, t, sprintf('*** TARGET DETECTED by D%d (alt %.0f m AGL) ***', i, agl), P);
            end
        end
    end

    % --- 6. radio graph, reachability, algebraic connectivity
    [A, Am, Qlink, shadow] = commGraph(P, W, S.pos, S.alive, shadow);
    on = [S.alive; true];
    [reach, parent] = bfsFrom(Am, N+1);        % topology = links with PER < 50 %
    lam2 = algebraicConnectivity(Am, on);

    % --- 7. multi-hop alert propagation (one hop per step)
    if S.detected
        known  = [S.informed; S.gcsInformed] & on;
        heard  = any(A(:, known), 2) & on;
        newInf = heard(1:N) & ~S.informed;
        for i = find(newInf)'
            S.informed(i) = true;
        end
        if heard(N+1) && ~S.gcsInformed
            S.gcsInformed = true; T_report = t; nextAssign = 0;
            S = logEv(S, t, sprintf('*** GCS RECEIVED ALERT (%.1f s after detection) ***', t - T_detect), P);
        end
        if isnan(T_all) && all(S.informed(S.alive))
            T_all = t;
            S = logEv(S, t, 'All surviving drones informed', P);
        end
    end

    % live link from the drone over the target to the GCS (for display)
    livePath = zeros(0,3);
    if S.gcsInformed
        src = find(S.alive & (S.mode == M.LOITER), 1);
        if ~isempty(src) && reach(src)
            nodes = [S.pos; gcsXYZ(P, W)];
            u = src; livePath = nodes(u,:);
            while u ~= N+1 && parent(u) > 0
                u = parent(u); livePath(end+1,:) = nodes(u,:); %#ok<AGROW>
            end
        end
    end

    % --- 8. metrics
    H.pos(:,:,k) = S.pos;   H.mode(:,k) = S.mode;  H.soc(:,k) = 100*S.E/S.Emax;
    H.inf(:,k)   = S.informed;
    H.cov(k)     = 100*mean(S.covered(:));
    H.POS(k)     = 100*(1 - sum(S.q(:)));
    nAlive       = sum(S.alive);
    H.alive(k)   = nAlive;
    H.baseConn(k)= 100*sum(reach(1:N) & S.alive)/max(nAlive,1);
    H.informed(k)= 100*sum(S.informed & S.alive)/max(nAlive,1);
    H.lam2(k)    = lam2;
    H.nLinks(k)  = nnz(triu(Am));

    % --- 9. visualisation
    if P.visualize && (mod(k, P.drawEvery) == 0 || k == nSteps)
        V = updateVis(V, P, W, S, H, k, Am, Qlink, livePath);
        for ts = P.saveSnapshotsAt(:)'
            if abs(t - ts) < dt*P.drawEvery/2 + 1e-9
                saveFig(V.fig, sprintf('%s_t%03d.png', P.snapshotPrefix, round(ts)));
            end
        end
    end

    if P.stopOnReport && S.gcsInformed
        lastStep = k; break
    end
end

%% ------------------------------------------------------------------ wrap-up
if lastStep < nSteps
    H = trimHistory(H, lastStep);
end
if P.visualize && isfield(V,'vw') && ~isempty(V.vw)
    try, close(V.vw); catch, end
end

R.params = P; R.world = W; R.hist = H; R.events = S.events;
R.T_detect = T_detect; R.T_report = T_report; R.T_allInformed = T_all;
R.found = S.detected; R.reported = S.gcsInformed; R.finder = S.finder;
R.coverageAtDetect = covAtDetect; R.finalCoverage = H.cov(end); R.finalPOS = H.POS(end);
R.numFailed = sum(~S.alive); R.footprintRadius = rfp;
R.belief = S.q / max(sum(S.q(:)), eps); R.covered = S.covered;
R.meanLam2 = mean(H.lam2); R.meanBaseConn = mean(H.baseConn);

if P.verbose
    fprintf('\n--- Summary ---\n');
    fprintf('Target detected : %s', yesno(R.found));
    if R.found, fprintf(' at t = %.1f s by D%d (coverage then %.1f %%)', T_detect, S.finder, covAtDetect); end
    fprintf('\nGCS alerted     : %s', yesno(R.reported));
    if R.reported, fprintf(' at t = %.1f s (report latency %.1f s)', T_report, T_report - T_detect); end
    fprintf('\nFinal coverage  : %.1f %% | cumulative POS %.1f %%\n', R.finalCoverage, R.finalPOS);
    fprintf('Failed drones   : %d | mean GCS-connected %.1f %% | mean lambda2 %.2f\n', R.numFailed, R.meanBaseConn, R.meanLam2);
end

if P.enableCSVLog, writeCSV(P, H); end
if P.visualize && P.showSummary
    fs = summaryFigure(R, P);
    if P.saveSummary, saveFig(fs, P.summaryFile); end
end
end

%% ========================================================================
%% Parameters
%% ========================================================================
function P = defaultParams()
P.numDrones     = 10;
P.simTime       = 420;      % s
P.dt            = 0.25;     % s
P.areaXY        = 400;      % m (square)
P.maxAlt        = 90;       % m
P.searchAlt     = 25;       % m above ground while searching
P.relayAlt      = 40;       % m above ground for relays / mule
P.vmax          = 6;        % m/s horizontal
P.vzMax         = 3;        % m/s vertical
P.amax          = 4;        % m/s^2
P.tauV          = 0.8;      % s velocity time constant
P.slowRadius    = 12;       % m, start slowing down near a goal
% sensing
P.camFOVdeg     = 60;       % camera field of view
P.pd0           = 0.95;     % detection prob. of one look at zero altitude
P.pdScale       = 120;      % m, altitude scale of detection decay
P.lookInterval  = 2;        % s between statistically independent looks
P.minSenseAGL   = 6;        % m, below this the camera is not used
% search prior
P.cellSize      = 10;       % m belief-grid resolution
P.lkpSigma      = 60;       % m spread of the LKP Gaussian
P.priorMix      = 0.7;      % weight of the LKP Gaussian vs uniform
P.searchPattern = 'bayes';  % bayes | lawnmower | spiral | random
P.sweepOverlap  = 0.9;      % lane spacing = overlap * footprint diameter
P.bayesLambda   = 80;       % m, distance discount in the greedy Bayes planner
P.replanTime    = 12;       % s
% radio (2.4 GHz mesh, e.g. XBee-class)
P.txPower       = 10;       % dBm
P.PL0           = 40;       % dB at 1 m
P.pathLossExp   = 2.7;
P.noiseFloor    = -95;      % dBm
P.snr50         = 10;       % dB where PER = 50 %
P.perSlope      = 0.8;      % 1/dB
P.shadowSigma   = 3;        % dB log-normal shadowing
P.shadowTau     = 4;        % s correlation time of shadowing
P.obstacleLoss  = 15;       % dB extra loss if a building blocks the line of sight
P.gcsGain       = 6;        % dB extra antenna gain at the GCS
P.gcsMastHeight = 6;        % m
% energy (scaled so that energy management is visible in a 7-minute mission)
P.batteryJ      = 60000;    % J
P.Phover        = 160;      % W
P.dragCoef      = 0.3;      % P = Phover*(1 + c*(v/vmax)^2) + Pcomm
P.Pcomm         = 5;        % W
P.chargeW       = 1200;     % W (fast charge / battery swap)
P.socInitRange  = [0.55 1.0];
P.reserveFrac   = 0.06;
P.returnMargin  = 1.3;
% failures and avoidance
P.mtbf          = 1500;     % s mean time between hardware failures (per drone)
P.numObstacles  = 6;
P.enableInspection = true;
P.sepDist       = 10;  P.kSep = 6;
P.obsInfluence  = 12;  P.kObs = 8;
% mission logic
P.assignEvery   = 2;        % s
% output
P.visualize     = true;
P.figVisible    = 'on';
P.drawEvery     = 2;
P.trailLen      = 80;
P.saveSnapshotsAt = [];
P.snapshotPrefix  = 'swarm';
P.showSummary   = true;
P.saveSummary   = true;
P.summaryFile   = 'swarm_v2_summary.png';
P.recordVideo   = false;
P.videoFile     = 'swarm_v2.mp4';
P.enableCSVLog  = false;
P.csvFile       = 'swarm_v2_log.csv';
P.stopOnReport  = false;
P.seed          = 7;
P.verbose       = true;
end

function P = paramDialog(P)
f = {'numDrones','simTime','searchPattern','areaXY','searchAlt','vmax','mtbf', ...
     'numObstacles','seed','recordVideo','enableCSVLog'};
prompts = {'Number of drones', 'Simulated time (s)', ...
    'Search strategy: bayes | lawnmower | spiral | random', 'Area side (m)', ...
    'Search altitude above ground (m)', 'Max speed (m/s)', ...
    'Mean time between failures per drone (s)', 'Number of buildings', ...
    'Random seed', 'Record MP4 video (0/1)', 'Write CSV log (0/1)'};
def = cell(size(f));
for k = 1:numel(f)
    v = P.(f{k});
    if ischar(v), def{k} = v; else, def{k} = num2str(v); end
end
a = inputdlg(prompts, 'UAV Rescue Swarm v2 - parameters', 1, def);
if isempty(a), P = []; return; end
for k = 1:numel(f)
    if ischar(P.(f{k})), P.(f{k}) = lower(strtrim(a{k}));
    else, P.(f{k}) = str2double(a{k}); end
end
P.recordVideo = logical(P.recordVideo); P.enableCSVLog = logical(P.enableCSVLog);
end

function M = modes()
M.SEARCH = 1; M.RETURN = 2; M.CHARGE = 3; M.INSPECT = 4; M.MULE = 5;
M.RELAY = 6; M.LOITER = 7; M.OVERWATCH = 8; M.FAILED = 9;
end

function C = modeColors()
C = [0.30 0.65 1.00;   % search    blue
     1.00 0.62 0.15;   % return    orange
     1.00 0.90 0.25;   % charge    yellow
     0.78 0.50 1.00;   % inspect   violet
     1.00 0.40 0.80;   % mule      pink
     0.15 1.00 0.85;   % relay     teal
     0.25 1.00 0.25;   % loiter    green
     0.70 1.00 0.45;   % overwatch lime
     1.00 0.22 0.22];  % failed    red
end

function n = modeNames()
n = {'Search','Return','Charging','Inspect','Data mule','Relay','Loiter (beacon)','Overwatch','Failed'};
end

%% ========================================================================
%% World
%% ========================================================================
function W = buildWorld(P)
L = P.areaXY;
W.L = L;
W.base = [20 20];
W.chargers = [20 20; L-25 L-25];
% terrain: sum of Gaussian hills  [cx cy amplitude sigma]
nH = 5;
W.hills = [L*rand(nH,2), 8 + 18*rand(nH,1), 45 + 60*rand(nH,1)];
g = linspace(0, L, 81);
[W.GX, W.GY] = meshgrid(g, g);
W.GZ = groundAt(W, W.GX, W.GY);
% belief grid
cx = (P.cellSize/2):P.cellSize:L;
W.cx = cx;
[W.CX, W.CY] = meshgrid(cx, cx);
% last known position and prior
W.lkp = L*(0.45 + 0.4*rand(1,2));
G = exp(-((W.CX - W.lkp(1)).^2 + (W.CY - W.lkp(2)).^2) / (2*P.lkpSigma^2));
G = G / sum(G(:));
prior = P.priorMix*G + (1 - P.priorMix)/numel(G);
% true target drawn from the prior (so the prior is "honest")
c = cumsum(prior(:)); kk = find(c >= rand*c(end), 1);
tx = min(max(W.CX(kk) + (rand-0.5)*P.cellSize, 1), L-1);
ty = min(max(W.CY(kk) + (rand-0.5)*P.cellSize, 1), L-1);
W.target = [tx, ty, groundAt(W, tx, ty)];
% buildings (cylinders), kept away from base, chargers and target
W.obs = struct('c', {}, 'r', {}, 'base', {}, 'top', {});
tries = 0;
while numel(W.obs) < P.numObstacles && tries < 500
    tries = tries + 1;
    c0 = 30 + (L-60)*rand(1,2); r0 = 8 + 10*rand;
    ok = norm(c0 - W.base) > 70 && norm(c0 - W.chargers(2,:)) > 40 && ...
         norm(c0 - W.target(1:2)) > r0 + 25;
    for j = 1:numel(W.obs)
        ok = ok && norm(c0 - W.obs(j).c) > r0 + W.obs(j).r + 20;
    end
    if ok
        b0 = groundAt(W, c0(1), c0(2));
        W.obs(end+1) = struct('c', c0, 'r', r0, 'base', b0, 'top', b0 + 18 + 35*rand); %#ok<AGROW>
    end
end
% the target cannot be inside a building
for j = 1:numel(W.obs)
    inside = (W.CX - W.obs(j).c(1)).^2 + (W.CY - W.obs(j).c(2)).^2 < W.obs(j).r^2;
    prior(inside) = 0;
end
W.prior = prior / sum(prior(:));
end

function z = groundAt(W, x, y)
z = zeros(size(x));
for h = 1:size(W.hills,1)
    z = z + W.hills(h,3) * exp(-((x - W.hills(h,1)).^2 + (y - W.hills(h,2)).^2) / (2*W.hills(h,4)^2));
end
end

function p = gcsXYZ(P, W)
p = [W.base, groundAt(W, W.base(1), W.base(2)) + P.gcsMastHeight];
end

function [d, idx] = nearestCharger(W, xy)
dd = sqrt(sum((W.chargers - xy).^2, 2));
[d, idx] = min(dd);
end

%% ========================================================================
%% Search planning
%% ========================================================================
function wp = planSearch(P, W, rfp)
N = P.numDrones; L = W.L; sw = 2*rfp*P.sweepOverlap; margin = 0.6*rfp;
wp = cell(N,1);
switch lower(P.searchPattern)
    case 'lawnmower'   % boustrophedon decomposition: one vertical strip per drone
        stripW = L/N;
        for i = 1:N
            x0 = (i-1)*stripW; nl = max(1, ceil(stripW/sw)); lane = stripW/nl;
            xs = x0 + lane/2 + (0:nl-1)*lane;
            w = zeros(0,2);
            for kk = 1:nl
                if mod(kk,2) == 1, ys = [margin, L-margin]; else, ys = [L-margin, margin]; end
                w = [w; xs(kk) ys(1); xs(kk) ys(2)]; %#ok<AGROW>
            end
            wp{i} = w;
        end
    case 'spiral'      % N interleaved Archimedean arms centred on the LKP
        b = sw*N/(2*pi);
        corners = [0 0; L 0; 0 L; L L];
        rmax = max(sqrt(sum((corners - W.lkp).^2, 2)));
        for i = 1:N
            phi = 2*pi*(i-1)/N; th = 0; w = zeros(0,2);
            while b*th < rmax
                r = b*th;
                pnt = W.lkp + r*[cos(th+phi), sin(th+phi)];
                if all(pnt >= margin) && all(pnt <= L-margin), w(end+1,:) = pnt; end %#ok<AGROW>
                th = th + min(0.5, 15/max(r, 5));
            end
            wp{i} = w;
        end
    case 'random'
        for i = 1:N, wp{i} = margin + (L-2*margin)*rand(30,2); end
    otherwise          % 'bayes' : goals are chosen online
        for i = 1:N, wp{i} = zeros(0,2); end
end
end

function g = bayesGoal(i, S, W, P, M, rfp)
% greedy information-seeking goal: argmax  b(c) * exp(-dist/lambda),
% skipping cells already claimed by other searching drones
b = S.q / max(sum(S.q(:)), eps);
d = sqrt((W.CX - S.pos(i,1)).^2 + (W.CY - S.pos(i,2)).^2);
score = b .* exp(-d / P.bayesLambda);
others = find(S.alive & S.mode == M.SEARCH);
for j = others'
    if j == i, continue; end
    dj = (W.CX - S.goal(j,1)).^2 + (W.CY - S.goal(j,2)).^2;
    score(dj < (2*rfp)^2) = 0;
end
[~, kk] = max(score(:));
g = [W.CX(kk), W.CY(kk)];
end

%% ========================================================================
%% Mission logic
%% ========================================================================
function [S, msg] = failDrone(S, P, M, i, why)
S.alive(i) = false; S.mode(i) = M.FAILED;
msg = sprintf('D%d FAILED (%s)', i, why);
% hand over unfinished waypoints and (optionally) inspect the crash site
cand = find(S.alive & S.mode == M.SEARCH);
if isempty(cand), return; end
d = sqrt(sum((S.pos(cand,1:2) - S.pos(i,1:2)).^2, 2));
[~, m] = min(d); j = cand(m);
rest = S.wp{i}(min(S.wpIdx(i), size(S.wp{i},1)+1):end, :);
S.wp{j} = [S.wp{j}(S.wpIdx(j):end, :); rest];
S.wpIdx(j) = 1;
if P.enableInspection
    S.mode(j) = M.INSPECT; S.inspectPos(j,:) = S.pos(i,1:2);
    msg = sprintf('%s -> D%d inspects crash site and takes over its area', msg, j);
elseif ~isempty(rest)
    msg = sprintf('%s -> D%d takes over its area', msg, j);
end
end

function slots = relaySlots(P, W)
% points on the straight line target -> GCS, spaced 0.8 * drone-drone range,
% until the last one is comfortably inside the GCS range
tgt = W.target(1:2); v = W.base - tgt; Ltot = norm(v); u = v / max(Ltot, eps);
slots = zeros(0,2); d = 0;
while Ltot - d > 0.8*P.rangeGCS
    d = d + P.relaySpacing;
    slots(end+1,:) = tgt + min(d, Ltot)*u; %#ok<AGROW>
end
end

function S = assignRoles(S, P, W, M)
avail = find(S.alive & S.informed & ~ismember(S.mode, [M.RETURN M.CHARGE M.FAILED M.INSPECT]));
if isempty(avail), return; end
if ~S.gcsInformed
    % finder becomes a data mule (store-carry-forward) towards the GCS
    if ismember(S.finder, avail)
        S.mode(S.finder) = M.MULE;
        avail(avail == S.finder) = [];
    end
else
    % the informed drone closest to the target stays as a beacon
    d = sqrt(sum((S.pos(avail,1:2) - W.target(1:2)).^2, 2));
    [~, m] = min(d);
    S.mode(avail(m)) = M.LOITER; avail(m) = [];
end
slots = relaySlots(P, W);
for s = 1:size(slots,1)
    if isempty(avail), break; end
    d = sqrt(sum((S.pos(avail,1:2) - slots(s,:)).^2, 2));
    [~, m] = min(d);
    S.mode(avail(m)) = M.RELAY; S.slotPos(avail(m),:) = slots(s,:);
    avail(m) = [];
end
S.mode(avail) = M.OVERWATCH;
end

function [S, ev] = updateGoals(S, P, W, M, rfp, dt)
ev = {};
N = size(S.pos,1); tgt = W.target(1:2);
for i = 1:N
    p = S.pos(i,:);
    switch S.mode(i)
        case M.SEARCH
            useBayes = strcmpi(P.searchPattern, 'bayes') || S.wpIdx(i) > size(S.wp{i},1);
            if useBayes
                S.goalTimer(i) = S.goalTimer(i) - dt;
                if S.goalTimer(i) <= 0 || norm(S.goal(i,1:2) - p(1:2)) < 0.5*rfp
                    S.goal(i,1:2) = bayesGoal(i, S, W, P, M, rfp);
                    S.goalTimer(i) = P.replanTime;
                end
                gxy = S.goal(i,1:2);
            else
                gxy = S.wp{i}(S.wpIdx(i),:);
                if norm(gxy - p(1:2)) < max(3, 0.25*rfp)
                    S.wpIdx(i) = S.wpIdx(i) + 1;
                    if S.wpIdx(i) <= size(S.wp{i},1), gxy = S.wp{i}(S.wpIdx(i),:); end
                end
            end
            S.goal(i,:) = [gxy, groundAt(W, p(1), p(2)) + P.searchAlt];
        case M.RETURN
            pad = W.chargers(S.charger(i),:) + 3*[cos(2*pi*i/N), sin(2*pi*i/N)];
            hd = norm(pad - p(1:2));
            zg = groundAt(W, pad(1), pad(2));
            if hd < 8, gz = zg + 0.3; else, gz = groundAt(W, p(1), p(2)) + P.searchAlt; end
            S.goal(i,:) = [pad, gz];
            if hd < 3 && p(3) - zg < 1.5
                S.mode(i) = M.CHARGE; S.vel(i,:) = 0; S.pos(i,:) = [pad, zg + 0.3];
                ev{end+1} = sprintf('D%d landed on pad %d, charging', i, S.charger(i)); %#ok<AGROW>
            end
        case M.INSPECT
            gxy = S.inspectPos(i,:);
            S.goal(i,:) = [gxy, groundAt(W, gxy(1), gxy(2)) + 12];
            if norm(gxy - p(1:2)) < 4
                S.mode(i) = M.SEARCH; S.goalTimer(i) = 0;
                ev{end+1} = sprintf('D%d inspected crash site, resumes search', i); %#ok<AGROW>
            end
        case M.MULE
            S.goal(i,:) = [W.base, groundAt(W, W.base(1), W.base(2)) + P.relayAlt];
        case M.RELAY
            gxy = S.slotPos(i,:);
            S.goal(i,:) = [gxy, groundAt(W, gxy(1), gxy(2)) + P.relayAlt];
        case M.LOITER
            S.goal(i,:) = [tgt, W.target(3) + 15];
        case M.OVERWATCH
            S.orbitAng(i) = S.orbitAng(i) + (0.6*P.vmax/30)*dt;
            gxy = tgt + 30*[cos(S.orbitAng(i)), sin(S.orbitAng(i))];
            S.goal(i,:) = [gxy, groundAt(W, gxy(1), gxy(2)) + 25];
    end
end
end

%% ========================================================================
%% Dynamics: velocity tracking + artificial potential fields
%% ========================================================================
function S = stepDynamics(S, P, W, M, dt)
N = size(S.pos,1);
nearPad = false(N,1);
for i = 1:N
    nearPad(i) = nearestCharger(W, S.pos(i,1:2)) < 15 && ismember(S.mode(i), [M.RETURN M.CHARGE]);
end
air = S.alive & ~nearPad;
for i = 1:N
    if S.mode(i) == M.FAILED                         % falls to the ground
        zg = groundAt(W, S.pos(i,1), S.pos(i,2));
        if S.pos(i,3) > zg + 0.3
            S.vel(i,:) = [0.9*S.vel(i,1:2), -8];
            S.pos(i,:) = S.pos(i,:) + S.vel(i,:)*dt;
            S.pos(i,3) = max(S.pos(i,3), zg + 0.3);
        else
            S.vel(i,:) = 0;
        end
        continue
    end
    if S.mode(i) == M.CHARGE, S.vel(i,:) = 0; continue; end
    p = S.pos(i,:); g = S.goal(i,:);
    e = g(1:2) - p(1:2); d = norm(e);
    vxy = P.vmax * min(1, d/P.slowRadius) * e / max(d, 1e-6);
    vz  = max(-P.vzMax, min(P.vzMax, 0.8*(g(3) - p(3))));
    vdes = [vxy, vz];
    a = (vdes - S.vel(i,:)) / P.tauV;
    if air(i)
        % separation (repulsive potential between drones)
        for j = find(air)'
            if j == i, continue; end
            r = p - S.pos(j,:); dd = norm(r);
            if dd < P.sepDist && dd > 1e-3
                a = a + P.kSep * (P.sepDist - dd)/P.sepDist * r/dd;
            end
        end
        % buildings: radial repulsion + tangential slide (avoids local minima)
        for kk = 1:numel(W.obs)
            o = W.obs(kk);
            if p(3) < o.top + 4
                r = p(1:2) - o.c; dd = norm(r); clr = dd - o.r;
                if clr < P.obsInfluence
                    n = r / max(dd, 1e-6);
                    push = P.kObs * (P.obsInfluence - clr) / P.obsInfluence;
                    tng = [-n(2), n(1)];
                    if dot(tng, vdes(1:2)) < 0, tng = -tng; end
                    a(1:2) = a(1:2) + push*n + 0.6*push*tng;
                end
            end
        end
    end
    an = norm(a); if an > P.amax, a = a/an*P.amax; end
    v = S.vel(i,:) + a*dt;
    sxy = norm(v(1:2)); if sxy > P.vmax, v(1:2) = v(1:2)/sxy*P.vmax; end
    v(3) = max(-P.vzMax, min(P.vzMax, v(3)));
    p = p + v*dt;
    p(1:2) = min(max(p(1:2), 0), W.L);
    zg = groundAt(W, p(1), p(2));
    if p(3) < zg + 0.3, p(3) = zg + 0.3; v(3) = max(v(3), 0); end
    p(3) = min(p(3), P.maxAlt);
    for kk = 1:numel(W.obs)                          % hard constraint
        o = W.obs(kk); r = p(1:2) - o.c; dd = norm(r);
        if dd < o.r + 0.5 && p(3) < o.top
            p(1:2) = o.c + (o.r + 0.5) * r / max(dd, 1e-6);
        end
    end
    S.pos(i,:) = p; S.vel(i,:) = v;
end
end

function Pw = powerDraw(P, v)
Pw = P.Phover * (1 + P.dragCoef * (norm(v(1:2))/P.vmax)^2) + P.Pcomm;
end

%% ========================================================================
%% Radio model and graph theory
%% ========================================================================
function R = linkRange(P, extraGain)
% distance at which mean SNR = snr50 + 3.7 dB  (PER ~ 5 %)
margin = P.snr50 + log(19)/P.perSlope;
R = 10^((P.txPower + extraGain - P.PL0 - P.noiseFloor - margin) / (10*P.pathLossExp));
end

function [A, Am, Q, sh] = commGraph(P, W, pos, alive, sh)
% A  : packets actually delivered this step (Bernoulli draw with 1-PER)
% Am : topology used for metrics (links whose current PER < 50 %)
% sh : log-normal shadowing, AR(1) correlated in time (Gudmundson model)
N = size(pos,1); n = N + 1;
X = [pos; gcsXYZ(P, W)];
on = double([alive; true]);
dx = X(:,1) - X(:,1)'; dy = X(:,2) - X(:,2)'; dz = X(:,3) - X(:,3)';
D  = max(sqrt(dx.^2 + dy.^2 + dz.^2), 1);
rho = exp(-P.dt / P.shadowTau);
w = P.shadowSigma * randn(n); w = triu(w,1); w = w + w';
sh = rho*sh + sqrt(1 - rho^2)*w;
PL = P.PL0 + 10*P.pathLossExp*log10(D) + sh + P.obstacleLoss*losBlocked(X, W);
gain = zeros(n); gain(n,:) = P.gcsGain; gain(:,n) = P.gcsGain; gain(n,n) = 0;
SNR = P.txPower + gain - PL - P.noiseFloor;
PER = 1 ./ (1 + exp(P.perSlope*(SNR - P.snr50)));
U = rand(n); U = triu(U,1); U = U + U';
A = (U > PER) & ((on*on') > 0);
A(logical(eye(n))) = false;
Am = (PER < 0.5) & ((on*on') > 0);
Am(logical(eye(n))) = false;
Q = (1 - PER) .* Am;
end

function B = losBlocked(X, W)
n = size(X,1); B = false(n);
x1 = X(:,1); y1 = X(:,2); z1 = X(:,3);
DX = x1' - x1; DY = y1' - y1; DZ = z1' - z1;
L2 = max(DX.^2 + DY.^2, 1e-9);
for k = 1:numel(W.obs)
    c = W.obs(k).c;
    T  = ((c(1) - x1).*DX + (c(2) - y1).*DY) ./ L2;
    T  = min(max(T, 0), 1);
    dist = sqrt((x1 + T.*DX - c(1)).^2 + (y1 + T.*DY - c(2)).^2);
    zAt  = z1 + T.*DZ;
    B = B | (dist < W.obs(k).r & zAt < W.obs(k).top);
end
end

function [reach, parent] = bfsFrom(A, src)
n = size(A,1); reach = false(n,1); parent = zeros(n,1);
reach(src) = true; queue = src;
while ~isempty(queue)
    u = queue(1); queue(1) = [];
    nb = find(A(u,:) & ~reach');
    reach(nb) = true; parent(nb) = u; queue = [queue, nb]; %#ok<AGROW>
end
end

function l2 = algebraicConnectivity(A, on)
idx = find(on); As = double(A(idx, idx));
if numel(idx) < 2, l2 = 0; return; end
Lp = diag(sum(As,2)) - As;
ev = sort(real(eig(Lp)));
l2 = max(ev(2), 0);
end

%% ========================================================================
%% Visualisation
%% ========================================================================
function V = initVis(P, W, S, rfp)
N = P.numDrones; C = modeColors(); names = modeNames();
bg = [0.055 0.063 0.085]; axbg = [0.085 0.095 0.125]; fg = [0.86 0.89 0.95];
V.bg = bg; V.axbg = axbg; V.fg = fg;
V.fig = figure('Name', 'UAV Rescue Swarm v2', 'NumberTitle', 'off', 'Color', bg, ...
    'Position', [20 30 1600 900], 'Visible', P.figVisible, 'InvertHardcopy', 'off');

% ---------- header
hx = axes('Parent', V.fig, 'Position', [0 0.935 1 0.065], 'Visible', 'off'); hold(hx, 'on');
xlim(hx, [0 1]); ylim(hx, [0 1]);
text(hx, 0.015, 0.55, 'UAV Search-and-Rescue Swarm', 'Color', fg, 'FontSize', 17, 'FontWeight', 'bold');
text(hx, 0.015, 0.08, sprintf('%d drones  |  strategy: %s  |  footprint r = %.0f m  |  radio range %.0f m (drone) / %.0f m (GCS)', ...
    N, upper(P.searchPattern), rfp, P.rangeDD, P.rangeGCS), 'Color', [0.6 0.66 0.78], 'FontSize', 10);
V.hHud = text(hx, 0.985, 0.4, '', 'Color', fg, 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');

% ---------- 3-D scene
ax = axes('Parent', V.fig, 'Position', [0.035 0.06 0.50 0.85], 'Color', axbg, ...
    'XColor', fg, 'YColor', fg, 'ZColor', fg, 'FontSize', 9);
hold(ax, 'on'); grid(ax, 'on');
surf(ax, W.GX, W.GY, W.GZ, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax, terrainMap());
caxis(ax, [min(W.GZ(:)) - 2, max(W.GZ(:)) + 8]);
th = linspace(0, 2*pi, 29);
for kk = 1:numel(W.obs)
    o = W.obs(kk);
    [cx, cy, cz] = cylinder(o.r, 28);
    surf(ax, cx + o.c(1), cy + o.c(2), o.base - 2 + cz*(o.top - o.base + 2), ...
        'FaceColor', [0.52 0.55 0.62], 'EdgeColor', 'none');
    patch(ax, o.c(1) + o.r*cos(th), o.c(2) + o.r*sin(th), o.top*ones(size(th)), [0.72 0.75 0.82], 'EdgeColor', 'none');
end
g = gcsXYZ(P, W);
plot3(ax, [g(1) g(1)], [g(2) g(2)], [g(3) - P.gcsMastHeight, g(3)], '-', 'Color', [0.9 0.9 0.9], 'LineWidth', 2);
plot3(ax, g(1), g(2), g(3), '^', 'MarkerSize', 12, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0.2 0.2 0.2]);
text(ax, g(1) + 6, g(2) - 6, g(3) + 10, 'GCS', 'Color', 'w', 'FontWeight', 'bold');
for c = 1:size(W.chargers,1)
    zc = groundAt(W, W.chargers(c,1), W.chargers(c,2));
    plot3(ax, W.chargers(c,1), W.chargers(c,2), zc + 0.2, 's', 'MarkerSize', 13, ...
        'MarkerFaceColor', [1 0.85 0.1], 'MarkerEdgeColor', [0.4 0.3 0]);
    if c > 1, text(ax, W.chargers(c,1) - 30, W.chargers(c,2), zc + 8, sprintf('Pad %d', c), 'Color', [1 0.9 0.4]); end
end
zl = groundAt(W, W.lkp(1), W.lkp(2));
plot3(ax, W.lkp(1), W.lkp(2), zl + 0.5, 'o', 'MarkerSize', 9, 'LineWidth', 2, 'Color', [1 0.55 0.1]);
text(ax, W.lkp(1) + 6, W.lkp(2), zl + 6, 'LKP', 'Color', [1 0.6 0.2], 'FontWeight', 'bold');
V.hTgtTrue = plot3(ax, W.target(1), W.target(2), W.target(3) + 0.5, 'x', 'MarkerSize', 9, 'LineWidth', 1.5, 'Color', [0.95 0.95 0.95]);
V.hTgt = plot3(ax, W.target(1), W.target(2), W.target(3) + 1, 'p', 'MarkerSize', 20, ...
    'MarkerFaceColor', [1 0.15 0.15], 'MarkerEdgeColor', 'w', 'Visible', 'off');
V.hTgtRing = plot3(ax, W.target(1) + 18*cos(th), W.target(2) + 18*sin(th), W.target(3) + 1 + 0*th, '-', ...
    'Color', [1 0.3 0.3], 'LineWidth', 2, 'Visible', 'off');
V.lnkCol = [0.25 0.95 0.45; 0.95 0.85 0.25; 0.95 0.35 0.30];
for q = 1:3
    V.hLink(q) = plot3(ax, NaN, NaN, NaN, '-', 'Color', V.lnkCol(q,:), 'LineWidth', 1.6 - 0.4*q);
end
V.hPath = plot3(ax, NaN, NaN, NaN, '-', 'Color', [0.2 1 1], 'LineWidth', 3.5);
for i = 1:N
    V.hTrail(i) = plot3(ax, NaN, NaN, NaN, '-', 'Color', [0.55 0.75 1], 'LineWidth', 0.7);
    V.hFp(i)    = plot3(ax, NaN, NaN, NaN, '-', 'Color', [0.45 0.95 1], 'LineWidth', 1.1);
end
V.hDr = scatter3(ax, S.pos(:,1), S.pos(:,2), S.pos(:,3), 70*ones(N,1), C(S.mode,:), 'filled', 'MarkerEdgeColor', [0.05 0.05 0.05]);
for i = 1:N
    V.hLbl(i) = text(ax, S.pos(i,1), S.pos(i,2), S.pos(i,3) + 5, sprintf('%d', i), 'Color', 'w', 'FontSize', 8, 'FontWeight', 'bold');
end
hl = plot3(ax, NaN, NaN, NaN, 'o', 'MarkerFaceColor', C(1,:), 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
for m = 2:9
    hl(m) = plot3(ax, NaN, NaN, NaN, 'o', 'MarkerFaceColor', C(m,:), 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
end
hl(10) = plot3(ax, NaN, NaN, NaN, '-', 'Color', [0.2 1 1], 'LineWidth', 3);
hl(11) = plot3(ax, NaN, NaN, NaN, '-', 'Color', [0.45 0.95 1]);
try
    lg = legend(ax, hl, [names, {'Live link to GCS', 'Camera footprint'}], 'Location', 'northwest');
    set(lg, 'TextColor', fg, 'Color', axbg, 'EdgeColor', [0.3 0.33 0.4], 'FontSize', 8);
catch
end
xlim(ax, [0 W.L]); ylim(ax, [0 W.L]); zlim(ax, [0 P.maxAlt]);
daspect(ax, [1 1 0.45]); view(ax, -38, 34);
xlabel(ax, 'x (m)'); ylabel(ax, 'y (m)'); zlabel(ax, 'z (m)');
try, light(ax, 'Position', [-300 -500 600], 'Style', 'infinite'); lighting(ax, 'gouraud'); catch, end
V.ax3 = ax;

% ---------- belief map
axb = axes('Parent', V.fig, 'Position', [0.575 0.545 0.19 0.355]);
V.hBel = imagesc(axb, W.cx, W.cx, log10(W.prior + 1e-12));
set(axb, 'YDir', 'normal'); hold(axb, 'on'); axis(axb, 'equal'); axis(axb, 'tight');
colormap(axb, beliefMap()); caxis(axb, [-5.5 -1.5]);
for kk = 1:numel(W.obs)
    o = W.obs(kk); plot(axb, o.c(1) + o.r*cos(th), o.c(2) + o.r*sin(th), '-', 'Color', [0.75 0.78 0.85]);
end
plot(axb, W.base(1), W.base(2), '^', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
plot(axb, W.lkp(1), W.lkp(2), 'o', 'Color', [1 0.6 0.2], 'LineWidth', 1.5);
V.hBelFp = plot(axb, NaN, NaN, '-', 'Color', [0.55 1 1], 'LineWidth', 0.8);
V.hBelDr = scatter(axb, S.pos(:,1), S.pos(:,2), 22*ones(N,1), C(S.mode,:), 'filled', 'MarkerEdgeColor', 'k');
V.hBelTgt = plot(axb, W.target(1), W.target(2), 'p', 'MarkerSize', 13, 'MarkerFaceColor', [1 0.2 0.2], 'MarkerEdgeColor', 'w', 'Visible', 'off');
styleAx(axb, V, 'Belief  log_{10} P(target in cell)');
try
    cb = colorbar(axb); set(cb, 'Color', fg, 'FontSize', 8);
catch
end

% ---------- battery panel
axe = axes('Parent', V.fig, 'Position', [0.815 0.545 0.165 0.355]); hold(axe, 'on');
xlim(axe, [0 100]); ylim(axe, [0.4 N + 0.6]); set(axe, 'YDir', 'reverse');
for i = 1:N
    patch(axe, [0 100 100 0], [i-0.33 i-0.33 i+0.33 i+0.33], [0.16 0.18 0.23], 'EdgeColor', 'none');
    V.hBat(i) = patch(axe, [0 1 1 0], [i-0.33 i-0.33 i+0.33 i+0.33], [0.3 0.9 0.4], 'EdgeColor', 'none');
    V.hBatTxt(i) = text(axe, 98, i, '', 'Color', 'w', 'FontSize', 8, 'HorizontalAlignment', 'right', 'FontWeight', 'bold');
end
set(axe, 'YTick', 1:N, 'YTickLabel', arrayfun(@(i) sprintf('D%d', i), 1:N, 'UniformOutput', false));
styleAx(axe, V, 'Battery (%) and mode');

% ---------- time series
V.axN = axes('Parent', V.fig, 'Position', [0.575 0.30 0.19 0.17]); hold(V.axN, 'on');
V.hConn = plot(V.axN, NaN, NaN, '-', 'Color', [0.3 0.8 1], 'LineWidth', 1.8);
V.hInf  = plot(V.axN, NaN, NaN, '-', 'Color', [0.3 1 0.4], 'LineWidth', 1.8);
xlim(V.axN, [0 P.simTime]); ylim(V.axN, [0 105]);
styleAx(V.axN, V, 'Network (%)');
try, lg = legend(V.axN, {'GCS-connected (2 s avg)', 'Informed'}, 'Location', 'southeast'); set(lg, 'TextColor', fg, 'Color', axbg, 'FontSize', 7); catch, end

V.axS = axes('Parent', V.fig, 'Position', [0.815 0.30 0.165 0.17]); hold(V.axS, 'on');
V.hCov = plot(V.axS, NaN, NaN, '-', 'Color', [0.45 0.65 1], 'LineWidth', 1.8);
V.hPOS = plot(V.axS, NaN, NaN, '-', 'Color', [1 0.6 0.2], 'LineWidth', 1.8);
xlim(V.axS, [0 P.simTime]); ylim(V.axS, [0 105]);
styleAx(V.axS, V, 'Search progress (%)');
try, lg = legend(V.axS, {'Area covered', 'Prob. of success'}, 'Location', 'southeast'); set(lg, 'TextColor', fg, 'Color', axbg, 'FontSize', 7); catch, end

V.axL = axes('Parent', V.fig, 'Position', [0.575 0.06 0.19 0.16]); hold(V.axL, 'on');
V.hLam = plot(V.axL, NaN, NaN, '-', 'Color', [0.8 0.55 1], 'LineWidth', 1.6);
xlim(V.axL, [0 P.simTime]); ylim(V.axL, [-0.3 4]);
styleAx(V.axL, V, 'Algebraic connectivity \lambda_2');
xlabel(V.axL, 'time (s)');

axg = axes('Parent', V.fig, 'Position', [0.8 0.035 0.19 0.2], 'Visible', 'off'); hold(axg, 'on');
xlim(axg, [0 1]); ylim(axg, [0 1]);
text(axg, 0, 1, 'Event log', 'Color', fg, 'FontWeight', 'bold', 'FontSize', 10, 'VerticalAlignment', 'top');
V.hLog = text(axg, 0, 0.86, '', 'Color', [0.75 0.8 0.9], 'FontSize', 7.5, 'VerticalAlignment', 'top', 'Interpreter', 'none');

V.vw = [];
if P.recordVideo
    try
        V.vw = VideoWriter(P.videoFile, 'MPEG-4'); V.vw.FrameRate = 24; open(V.vw);
    catch
        warning('Video recording not available on this platform.'); V.vw = [];
    end
end
drawnow;
end

function styleAx(ax, V, ttl)
set(ax, 'Color', V.axbg, 'XColor', V.fg, 'YColor', V.fg, 'FontSize', 8, 'Box', 'on');
grid(ax, 'on');
try, set(ax, 'GridColor', [0.5 0.55 0.65]); catch, end
title(ax, ttl, 'Color', V.fg, 'FontSize', 10, 'FontWeight', 'bold');
end

function V = updateVis(V, P, W, S, H, k, A, Q, livePath)
M = modes(); C = modeColors(); short = {'SRCH','RTN','CHG','INSP','MULE','RLY','LOIT','OVW','FAIL'};
N = size(S.pos,1); t = H.t(k);
cols = C(S.mode,:);
sz = 60*ones(N,1); sz(S.informed) = 110;
set(V.hDr, 'XData', S.pos(:,1), 'YData', S.pos(:,2), 'ZData', S.pos(:,3), 'CData', cols, 'SizeData', sz);
X = [S.pos; gcsXYZ(P, W)];
[I, J] = find(triu(A));
qv = Q(sub2ind(size(Q), I, J));
bins = {qv >= 0.9, qv >= 0.6 & qv < 0.9, qv < 0.6};
for q = 1:3
    sel = bins{q};
    if any(sel)
        ii = I(sel); jj = J(sel); m = numel(ii);
        xs = [X(ii,1) X(jj,1) nan(m,1)]'; ys = [X(ii,2) X(jj,2) nan(m,1)]'; zs = [X(ii,3) X(jj,3) nan(m,1)]';
        set(V.hLink(q), 'XData', xs(:), 'YData', ys(:), 'ZData', zs(:));
    else
        set(V.hLink(q), 'XData', NaN, 'YData', NaN, 'ZData', NaN);
    end
end
if size(livePath,1) >= 2
    set(V.hPath, 'XData', livePath(:,1), 'YData', livePath(:,2), 'ZData', livePath(:,3));
else
    set(V.hPath, 'XData', NaN, 'YData', NaN, 'ZData', NaN);
end
th = linspace(0, 2*pi, 25); fpx = []; fpy = [];
k0 = max(1, k - P.trailLen);
for i = 1:N
    tr = squeeze(H.pos(i,:,k0:k));
    if size(tr,1) == 3 && size(tr,2) > 1
        set(V.hTrail(i), 'XData', tr(1,:), 'YData', tr(2,:), 'ZData', tr(3,:));
    end
    zg = groundAt(W, S.pos(i,1), S.pos(i,2)); agl = S.pos(i,3) - zg;
    if S.alive(i) && agl >= P.minSenseAGL && S.mode(i) ~= M.CHARGE
        r = agl*tand(P.camFOVdeg/2);
        cxr = S.pos(i,1) + r*cos(th); cyr = S.pos(i,2) + r*sin(th);
        set(V.hFp(i), 'XData', cxr, 'YData', cyr, 'ZData', groundAt(W, cxr, cyr) + 0.6);
        fpx = [fpx, cxr, NaN]; fpy = [fpy, cyr, NaN]; %#ok<AGROW>
    else
        set(V.hFp(i), 'XData', NaN, 'YData', NaN, 'ZData', NaN);
    end
    set(V.hLbl(i), 'Position', S.pos(i,:) + [0 0 5]);
    soc = max(0, 100*S.E(i)/S.Emax);
    set(V.hBat(i), 'Vertices', [0 i-0.33; soc i-0.33; soc i+0.33; 0 i+0.33], 'Faces', 1:4, 'FaceColor', socColor(soc/100));
    set(V.hBatTxt(i), 'String', sprintf('%s  %3.0f%%', short{S.mode(i)}, soc), 'Color', C(S.mode(i),:));
end
if S.detected
    set([V.hTgt V.hTgtRing V.hBelTgt], 'Visible', 'on'); set(V.hTgtTrue, 'Visible', 'off');
end
b = S.q / max(sum(S.q(:)), eps);
set(V.hBel, 'CData', log10(b + 1e-12));
set(V.hBelDr, 'XData', S.pos(:,1), 'YData', S.pos(:,2), 'CData', cols);
if isempty(fpx), fpx = NaN; fpy = NaN; end
set(V.hBelFp, 'XData', fpx, 'YData', fpy);
tt = H.t(1:k);
set(V.hConn, 'XData', tt, 'YData', smooth2s(H.baseConn(1:k), P.dt));
set(V.hInf,  'XData', tt, 'YData', H.informed(1:k));
set(V.hCov,  'XData', tt, 'YData', H.cov(1:k));
set(V.hPOS,  'XData', tt, 'YData', H.POS(1:k));
set(V.hLam,  'XData', tt, 'YData', H.lam2(1:k));
ev = S.events; ev = ev(max(1, numel(ev) - 8):end);
for e = 1:numel(ev), if numel(ev{e}) > 58, ev{e} = [ev{e}(1:56) '..']; end, end
set(V.hLog, 'String', ev);
hud = sprintf('t = %5.1f s   |   alive %d/%d   |   coverage %4.1f%%', t, sum(S.alive), N, H.cov(k));
if S.detected, hud = [hud, sprintf('   |   FOUND @ %.0f s', firstTrue(H.inf, S.finder, H.t))]; end
if S.gcsInformed, hud = [hud, '   |   GCS ALERTED']; end
set(V.hHud, 'String', hud);
drawnow;
if isfield(V, 'vw') && ~isempty(V.vw)
    writeVideo(V.vw, getframe(V.fig));
end
end

function y = smooth2s(x, dt)
w = max(1, round(2/dt)); y = filter(ones(1,w)/w, 1, x);
y(1:min(w,end)) = cumsum(x(1:min(w,end))) ./ (1:min(w,numel(x)));
end

function tf = firstTrue(INF, i, t)
k = find(INF(i,:), 1); if isempty(k), tf = NaN; else, tf = t(k); end
end

function c = socColor(s)
s = min(max(s, 0), 1);
if s > 0.5, c = [2*(1 - s), 0.85, 0.35]; else, c = [0.95, 1.7*s, 0.25]; end
c = min(max(c, 0), 1);
end

function m = terrainMap()
k = [0 0.35 0.65 1];
c = [0.10 0.26 0.16; 0.20 0.42 0.22; 0.50 0.50 0.30; 0.72 0.66 0.52];
x = linspace(0, 1, 128)';
m = [interp1(k, c(:,1), x), interp1(k, c(:,2), x), interp1(k, c(:,3), x)];
end

function m = beliefMap()
k = [0 0.3 0.6 0.85 1];
c = [0.03 0.03 0.10; 0.25 0.05 0.45; 0.80 0.20 0.35; 1.00 0.65 0.15; 1.00 1.00 0.75];
x = linspace(0, 1, 128)';
m = [interp1(k, c(:,1), x), interp1(k, c(:,2), x), interp1(k, c(:,3), x)];
end

function fs = summaryFigure(R, P)
H = R.hist; W = R.world; N = P.numDrones; C = modeColors(); names = modeNames();
fg = [0.86 0.89 0.95]; V.axbg = [0.085 0.095 0.125]; V.fg = fg;
fs = figure('Name', 'Swarm v2 - mission summary', 'NumberTitle', 'off', 'Color', [0.055 0.063 0.085], ...
    'Position', [60 60 1500 820], 'Visible', P.figVisible, 'InvertHardcopy', 'off');
th = linspace(0, 2*pi, 40);
a1 = subplot(2,3,1); hold(a1, 'on');
imagesc(a1, W.cx, W.cx, log10(R.belief + 1e-12)); set(a1, 'YDir', 'normal');
colormap(a1, beliefMap()); caxis(a1, [-5.5 -1.5]); axis(a1, 'equal'); axis(a1, 'tight');
lc = lines(N);
for i = 1:N, plot(a1, squeeze(H.pos(i,1,:)), squeeze(H.pos(i,2,:)), '-', 'Color', [lc(i,:) ], 'LineWidth', 0.8); end
for kk = 1:numel(W.obs), o = W.obs(kk); plot(a1, o.c(1) + o.r*cos(th), o.c(2) + o.r*sin(th), 'w-'); end
plot(a1, W.target(1), W.target(2), 'p', 'MarkerSize', 14, 'MarkerFaceColor', [1 0.2 0.2], 'MarkerEdgeColor', 'w');
plot(a1, W.base(1), W.base(2), '^', 'MarkerFaceColor', 'w', 'MarkerSize', 9);
styleAx(a1, V, 'Flight paths over final belief');

a2 = subplot(2,3,2); hold(a2, 'on');
plot(a2, H.t, H.cov, '-', 'Color', [0.45 0.65 1], 'LineWidth', 2);
plot(a2, H.t, H.POS, '-', 'Color', [1 0.6 0.2], 'LineWidth', 2);
markEvents(a2, R); ylim(a2, [0 105]); xlim(a2, [0 H.t(end)]);
styleAx(a2, V, 'Coverage and probability of success (%)');
try, lg = legend(a2, {'Coverage', 'POS'}, 'Location', 'southeast'); set(lg, 'TextColor', fg, 'Color', V.axbg); catch, end

a3 = subplot(2,3,3); hold(a3, 'on');
plot(a3, H.t, smooth2s(H.baseConn, P.dt), '-', 'Color', [0.3 0.8 1], 'LineWidth', 1.6);
plot(a3, H.t, H.informed, '-', 'Color', [0.3 1 0.4], 'LineWidth', 2);
markEvents(a3, R); ylim(a3, [0 105]); xlim(a3, [0 H.t(end)]);
styleAx(a3, V, 'GCS-connected and informed drones (%)');
try, lg = legend(a3, {'GCS-connected (2 s avg)', 'Informed'}, 'Location', 'southeast'); set(lg, 'TextColor', fg, 'Color', V.axbg); catch, end

a4 = subplot(2,3,4); hold(a4, 'on');
for i = 1:N, plot(a4, H.t, H.soc(i,:), '-', 'Color', lc(i,:), 'LineWidth', 1.1); end
ylim(a4, [0 105]); styleAx(a4, V, 'Battery state of charge per drone (%)'); xlim(a4, [0 H.t(end)]); xlabel(a4, 'time (s)');

a5 = subplot(2,3,5);
image(a5, H.t, 1:N, ind2rgb(H.mode, C)); set(a5, 'YDir', 'normal');
styleAx(a5, V, 'Mode timeline'); xlabel(a5, 'time (s)'); ylabel(a5, 'drone');
set(a5, 'YTick', 1:N);

a6 = subplot(2,3,6); axis(a6, 'off'); xlim(a6, [0 1]); ylim(a6, [0 1]);
s = {sprintf('Strategy: %s   (seed %d)', upper(P.searchPattern), P.seed), ...
     sprintf('Target detected: %s', tern(R.found, sprintf('t = %.1f s by D%d', R.T_detect, R.finder), 'no')), ...
     sprintf('Coverage at detection: %.1f %%', R.coverageAtDetect), ...
     sprintf('GCS alerted: %s', tern(R.reported, sprintf('t = %.1f s (latency %.1f s)', R.T_report, R.T_report - R.T_detect), 'no')), ...
     sprintf('All drones informed: %s', tern(~isnan(R.T_allInformed), sprintf('t = %.1f s', R.T_allInformed), 'no')), ...
     sprintf('Final coverage %.1f %%,  POS %.1f %%', R.finalCoverage, R.finalPOS), ...
     sprintf('Failed drones: %d', R.numFailed), ...
     sprintf('Mean GCS-connected %.1f %%,  mean lambda_2 = %.2f', R.meanBaseConn, R.meanLam2)};
text(a6, 0, 1.0, s, 'Color', fg, 'FontSize', 10.5, 'VerticalAlignment', 'top', 'Interpreter', 'none');
for m = 1:9
    patch(a6, [0 0.04 0.04 0] + 0.02 + 0.33*mod(m-1,3), [0 0 0.05 0.05] + 0.13 - 0.08*floor((m-1)/3), C(m,:), 'EdgeColor', 'none');
    text(a6, 0.07 + 0.33*mod(m-1,3), 0.155 - 0.08*floor((m-1)/3), names{m}, 'Color', fg, 'FontSize', 8);
end
end

function markEvents(ax, R)
yl = [0 105];
if ~isnan(R.T_detect), plot(ax, [R.T_detect R.T_detect], yl, '--', 'Color', [1 0.3 0.3]); end
if ~isnan(R.T_report), plot(ax, [R.T_report R.T_report], yl, '--', 'Color', [0.2 1 1]); end
end

function saveFig(f, fname)
try
    if exist('exportgraphics', 'file') == 2
        exportgraphics(f, fname, 'Resolution', 130, 'BackgroundColor', 'current');
    else
        print(f, fname, '-dpng', '-r110');
    end
catch err
    warning('Could not save %s: %s', fname, err.message);
end
end

%% ========================================================================
%% Utilities
%% ========================================================================
function S = logEv(S, t, msg, P)
line = sprintf('[%6.1f s] %s', t, msg);
S.events{end+1} = line;
if P.verbose, fprintf('%s\n', line); end
end

function H = trimHistory(H, k)
f = fieldnames(H);
for i = 1:numel(f)
    v = H.(f{i});
    if ndims(v) == 3, H.(f{i}) = v(:,:,1:k); else, H.(f{i}) = v(:,1:k); end
end
end

function writeCSV(P, H)
fid = fopen(P.csvFile, 'w');
fprintf(fid, 'time,drone,x,y,z,soc,mode,informed\n');
N = size(H.pos,1);
for k = 1:numel(H.t)
    for i = 1:N
        fprintf(fid, '%.2f,%d,%.2f,%.2f,%.2f,%.1f,%d,%d\n', H.t(k), i, H.pos(i,1,k), H.pos(i,2,k), H.pos(i,3,k), H.soc(i,k), H.mode(i,k), H.inf(i,k));
    end
end
fclose(fid);
fprintf('CSV log written to %s\n', P.csvFile);
end

function P = mergeStruct(P, U)
f = fieldnames(U);
for i = 1:numel(f), P.(f{i}) = U.(f{i}); end
end

function setSeed(s)
try, rng(s); catch, rand('state', s); randn('state', s); end %#ok<RAND>
end

function s = yesno(b)
if b, s = 'yes'; else, s = 'no'; end
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end

function T = run_monte_carlo(nRuns, strategies, outPrefix, seedBase)
% RUN_MONTE_CARLO  Compare search strategies over many random scenarios.
%   T = run_monte_carlo(30)      % 30 scenarios x 4 strategies (headless)
%
% The same seed gives the same terrain, buildings, LKP and casualty for every
% strategy (common random numbers), so the comparison is paired and fair.
% Results are written to <outPrefix>.csv and plotted.

if nargin < 1, nRuns = 30; end
if nargin < 2, strategies = {'bayes', 'lawnmower', 'spiral', 'random'}; end
if nargin < 3, outPrefix = 'mc_results'; end
if nargin < 4, seedBase = 1000; end

nS = numel(strategies);
Tdet = nan(nRuns, nS); Trep = nan(nRuns, nS); Cov = nan(nRuns, nS); Fail = nan(nRuns, nS);
for s = 1:nS
    for r = 1:nRuns
        R = drone_swarm_rescue_v2(struct('visualize', false, 'verbose', false, ...
            'seed', seedBase + r, 'searchPattern', strategies{s}, 'stopOnReport', true));
        Tdet(r,s) = R.T_detect; Trep(r,s) = R.T_report;
        Cov(r,s) = R.coverageAtDetect; Fail(r,s) = R.numFailed;
        fprintf('%-10s run %2d: detect %6.1f s  report %6.1f s\n', strategies{s}, r, R.T_detect, R.T_report);
    end
end

% ---- CSV
fid = fopen([outPrefix '.csv'], 'w');
fprintf(fid, 'strategy,seed,T_detect,T_report,coverage_at_detect,failed\n');
for s = 1:nS
    for r = 1:nRuns
        fprintf(fid, '%s,%d,%.2f,%.2f,%.2f,%d\n', strategies{s}, seedBase + r, Tdet(r,s), Trep(r,s), Cov(r,s), Fail(r,s));
    end
end
fclose(fid);

% ---- summary table
fprintf('\n%-10s %8s %14s %14s %10s\n', 'strategy', 'found %', 'T_detect (s)', 'latency (s)', 'cov@det %');
for s = 1:nS
    ok = ~isnan(Tdet(:,s));
    lat = Trep(:,s) - Tdet(:,s);
    fprintf('%-10s %8.0f %8.1f+-%4.1f %8.1f+-%4.1f %10.1f\n', strategies{s}, 100*mean(ok), ...
        mean(Tdet(ok,s)), ci95(Tdet(ok,s)), mean(lat(~isnan(lat))), ci95(lat(~isnan(lat))), mean(Cov(ok,s)));
end
T.strategies = strategies; T.Tdet = Tdet; T.Trep = Trep; T.Cov = Cov; T.Fail = Fail;

% ---- plots: mean +- 95 % CI and empirical CDF of detection time
f = figure('Name', 'Monte Carlo comparison', 'Color', 'w', 'Position', [80 80 1100 420]);
subplot(1,2,1); hold on; grid on;
m = zeros(1,nS); c = zeros(1,nS);
for s = 1:nS
    x = Tdet(~isnan(Tdet(:,s)), s); m(s) = mean(x); c(s) = ci95(x);
end
bar(1:nS, m, 0.6, 'FaceColor', [0.3 0.55 0.9]);
he = errorbar(1:nS, m, c, 'k.'); set(he, 'LineWidth', 1.5);
set(gca, 'XTick', 1:nS, 'XTickLabel', strategies);
ylabel('time to detect (s)'); title(sprintf('Mean detection time (\\pm 95%% CI, %d scenarios)', nRuns));
subplot(1,2,2); hold on; grid on;
for s = 1:nS
    x = sort(Tdet(:,s)); x(isnan(x)) = Inf;
    stairs([0; x(isfinite(x))], [0; (1:sum(isfinite(x)))'/nRuns], 'LineWidth', 2);
end
xlabel('time (s)'); ylabel('P(found by t)'); title('Empirical CDF of detection time');
legend(strategies, 'Location', 'southeast');
try, print(f, [outPrefix '.png'], '-dpng', '-r120'); catch, end
end

function c = ci95(x)
x = x(~isnan(x));
if numel(x) < 2, c = 0; else, c = 1.96*std(x)/sqrt(numel(x)); end
end

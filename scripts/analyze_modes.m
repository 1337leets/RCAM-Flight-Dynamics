function analyze_modes(Alat, Alon)
%ANALYZE_MODES  Identify the flight modes from the eigenvalues.
%
%   analyze_modes(Alat, Alon)
%
%   Prints the lateral and longitudinal modes with their natural
%   frequency, damping ratio, time constant and period, and draws a pole
%   map.

fprintf('\n==================== LATERAL MODES ====================\n');
lam_lat = eig(Alat);
print_modes(lam_lat, "lat");

fprintf('\n=================== LONGITUDINAL MODES =================\n');
lam_lon = eig(Alon);
print_modes(lam_lon, "lon");

% Pole map
figure;
tiledlayout(1,2);

nexttile;
plot(real(lam_lat), imag(lam_lat), 'x', 'LineWidth', 1.5, 'MarkerSize', 10);
grid on; axis equal;
xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
title('Lateral poles');

nexttile;
plot(real(lam_lon), imag(lam_lon), 'x', 'LineWidth', 1.5, 'MarkerSize', 10);
grid on; axis equal;
xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
title('Longitudinal poles');

end

function print_modes(lam, kind)
tol = 1e-6;

n = numel(lam);
wn   = nan(n,1);
zeta = nan(n,1);
tau  = nan(n,1);
T    = nan(n,1);
mode = strings(n,1);

isComplex = abs(imag(lam)) > tol;
realIdx   = find(~isComplex);
cplxIdx   = find(isComplex);

% Basic metrics
for i = 1:n
    sigma = real(lam(i));
    omega = abs(imag(lam(i)));
    wn(i) = abs(lam(i));

    if omega > tol
        % Complex pole
        zeta(i) = -sigma / wn(i);
        T(i) = 2*pi / omega;
    else
        % Real pole
        zeta(i) = NaN;
        if abs(sigma) > tol
            tau(i) = -1 / sigma;   % positive time constant when sigma < 0
        end
    end
end

% Mode labels
if kind == "lat"
    % The complex pair is the Dutch roll
    for i = cplxIdx.'
        mode(i) = "Dutch roll";
    end

    % Of the remaining real poles, the largest |lambda| is roll subsidence
    % and the smallest is the spiral mode
    if ~isempty(realIdx)
        [~, ord] = sort(abs(real(lam(realIdx))), 'descend');
        realSorted = realIdx(ord);

        if numel(realSorted) >= 1
            mode(realSorted(1)) = "Roll subsidence";
        end
        if numel(realSorted) >= 2
            mode(realSorted(end)) = "Spiral";
        end

        % Label anything in between
        for k = 2:max(1, numel(realSorted)-1)
            if mode(realSorted(k)) == ""
                mode(realSorted(k)) = "Real mode";
            end
        end
    end

else % kind == "lon"
    % Complex pairs: larger wn is the short period, smaller wn the phugoid
    if ~isempty(cplxIdx)
        % group the complex poles into conjugate pairs
        used = false(n,1);
        pairs = {};

        for i = cplxIdx.'
            if used(i), continue; end
            li = lam(i);
            j = find(~used & abs(lam - conj(li)) < 1e-6, 1);
            if isempty(j)
                pairs{end+1} = i; %#ok<AGROW>
                used(i) = true;
            else
                pairs{end+1} = [i j]; %#ok<AGROW>
                used([i j]) = true;
            end
        end

        pairWn = zeros(numel(pairs),1);
        for p = 1:numel(pairs)
            idx = pairs{p};
            pairWn(p) = abs(lam(idx(1)));
        end

        % Larger wn is the short period, smaller wn the phugoid
        [~, ord] = sort(pairWn, 'descend');
        if ~isempty(ord)
            idx = pairs{ord(1)};
            mode(idx) = "Short period";
        end
        if numel(ord) >= 2
            idx = pairs{ord(end)};
            mode(idx) = "Phugoid";
        end

        % Generic label for any remaining complex pair
        for p = 1:numel(pairs)
            idx = pairs{p};
            if all(mode(idx) == "")
                mode(idx) = "Complex mode";
            end
        end
    end

    % Real poles are only reported
    for i = realIdx.'
        mode(i) = "Real mode";
    end
end

% Fill in anything left unlabeled
for i = 1:n
    if mode(i) == ""
        if isComplex(i)
            mode(i) = "Complex mode";
        else
            mode(i) = "Real mode";
        end
    end
end

% Table output
fprintf('%-18s %-16s %10s %10s %12s %12s\n', ...
    'Mode', 'lambda', 'wn', 'zeta', 'tau', 'T');
fprintf('%s\n', repmat('-',1,80));

for i = 1:n
    if isnan(zeta(i))
        zetaStr = '   -';
    else
        zetaStr = sprintf('%10.4f', zeta(i));
    end

    if isnan(tau(i))
        tauStr = '      -';
    else
        tauStr = sprintf('%12.4f', tau(i));
    end

    if isnan(T(i))
        TStr = '      -';
    else
        TStr = sprintf('%12.4f', T(i));
    end

    fprintf('%-18s %s %10.4f %s %s %s\n', ...
        mode(i), lambda_str(lam(i)), wn(i), zetaStr, tauStr, TStr);
end

% Hint for the yaw damper design
if kind == "lat"
    dr_idx = find(mode == "Dutch roll", 1);
    if ~isempty(dr_idx) && ~isnan(zeta(dr_idx))
        fprintf('\nDutch roll damping ratio (zeta) ~ %.4f\n', zeta(dr_idx));
        fprintf('Compare with the CS-25 / MIL-F-8785C handling qualities limits.\n');
    end
end

end

function s = lambda_str(lam)
if abs(imag(lam)) < 1e-6
s = sprintf('%+10.4f', real(lam));
else
s = sprintf('%+7.4f%+7.4fj', real(lam), imag(lam));
end
end

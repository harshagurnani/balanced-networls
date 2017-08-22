clear();
figure();

% Params
nNeurons = 50;
dt = 20e-3;
tTotal = 1000;
lambda =1;
sigma = 0;
mu = 0.1;
signalFrequency = 0.001;

% Create readout weights of uniform length
x  = linspace(0, 2*pi, nNeurons+1);
x(end) = [];
D(1, :) = sin(x) * 10/nNeurons;
D(2, :) = cos(x) * 10/nNeurons;

% Create arrays
t = (dt : dt : tTotal)';
nT = numel(t);
V = zeros(nNeurons, nT);
V(:,1)=unifrnd(0,1,nNeurons,1);
R = zeros(nNeurons, nT);
xH = zeros(nT, 1);

% Create 2D signal
omega = t * signalFrequency * 2*pi;
X = [sin(omega+pi) cos(omega)]*0.5;
X(1:200, :) = 0;

isSpike = false(nNeurons, 1);
isSpikeAll = false(nNeurons, nT);

omega = D'*D + lambda*mu*eye(nNeurons);
vThresh = diag(omega)/2;
xdot = (X( 2:end,:)-X(1:end-1,:)) / dt;
for s = 2:nT
    dx = xdot(s-1,:);% (X(s, :)-X(s-1, :)) / dt;
    c = dx + X(s);
    
    noise = sigma*randn(nNeurons, 1);
    dVdT = -V(:, s-1) + (c*D)' - (omega)*isSpikeAll(:,s-1)/dt ;%+ noise;
    dRdT = -R(:, s-1) + lambda*isSpikeAll(:,s-1)/dt;
    V(:, s) = V(:, s-1) + dt*dVdT;
    R(:, s) = R(:, s-1) + dt*dRdT;
    
    isSpike = V(:, s) >= vThresh;
    indxspk = find( V(:, s) >= vThresh);
    if ~isempty(indxspk)
    isSpikeAll( indxspk(unidrnd(length(indxspk))), s) = true;%isSpike;
    end
end

XHat = D*R;

%% Plots

clf();

% Signal
ax(1) = subplot(3, 2, 1);
plot(t, X);
xlabel('time / s');
ylabel('x');
title('Signal');

% Readout weights
subplot(3, 2, 2);
x = D(1, :);
y = D(2, :);
line(x, y, 'color', 'k', 'lineStyle', 'none', 'marker', 'o');
axis equal
xlabel('W(x_1)');
ylabel('W(x_2)');
title('Weights');

% Voltage
ax(end+1) = subplot(3, 2, 3);
cla();
y = V';
plot(t, y)
xlabel('Time / s');
ylabel('V');

% Spikes
ax(end+1) = subplot(3, 2, 4);
for n = 1:nNeurons
    tSpike = find(isSpikeAll(n, :)) * dt;
    line(tSpike, n*ones(size(tSpike)), 'color', 'k', 'lineStyle', 'none', 'marker', '.', 'markerSize', 10);
end
title('Spikes');

% Readout
ax(end+1) = subplot(3, 2, 5);
plot(t, XHat)
title('Readouts');

% Readout error
clear h
subplot(3, 2, 6);
cla();
h(1) = line(X(:, 1), X(:, 2), 'color', 'k');
hold on
col = rg.color.dataToRgb(1:nT, [], 'jet');
h(2) = scatter(XHat(1, :), XHat(2, :), 5, col);
axis equal
title('Readout error');
legend(h, {'X', 'Xhat'});
xlabel('x_1');
ylabel('x_2');

linkaxes(ax, 'x');

%%

figure();
plot(t, R', 'color', 'k');
title('Rates');
xlabel('Time');
ylabel('Rate');
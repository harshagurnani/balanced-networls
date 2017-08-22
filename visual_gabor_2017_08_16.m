clear();
figure();
close all

% Params
nRfPositions = [5 15];
nRfAngles = 8;
gaborWidths = [10 4];
nNeurons = sum(nRfPositions.^2) * nRfAngles;
plotRfs = false;
dt = 1e-5;
tTotal = 1;
dtSignal = 1/30;
signalNSmooth = 2;
imageWidth = 31;

deadNeurons = [];
% signalFrequency = 1;

% % Lambda 1
lambda = 50;
sigma = .0004;
mu = 0.0001;

% Create gabor readout weights
if plotRfs, figure('colormap', rg.color.maps.seismic); end
% sply = floor(sqrt(nNeurons));
% splx = ceil(nNeurons / sply);
sply = 5;
splx = 5;

D = zeros(imageWidth^2, nNeurons);
clear muRF
n = 0;
imageHalfWidth = floor(imageWidth/2);

for w = 1:numel(gaborWidths)
    
    x = linspace(1, imageWidth, nRfPositions(w)+1);
    x = x(1:end-1)+(x(2)-x(1))/2;
    [xx,yy] = meshgrid(x);
    
    for p = 1:nRfPositions(w)^2
        for o = 1:nRfAngles
            n = n+1;
            h = rg.images.gabor( ...
                'wavelength', gaborWidths(w), ...
                'orientation', o*2*pi/nRfAngles, ...
                'center', [xx(p), yy(p)]-imageHalfWidth, ...
                'size', imageWidth*[1 1], ...
                'aspectRatio', 1, ...
                'frequencyBandwidth', 2.25, ...
                'phi', pi/2);
            h = real(h);
            h(abs(h)<0.1) = 0;
            %     [muRF(1), muRF(2)] = ind2sub(imageWidth*[1 1], rem(n-1, imageWidth^2)+1);
            %     muRF = ceil(rand(1, 2)*imageWidth);
            %     h = rg.helpers.normpdf2(xx, yy, muRF, 1);
            %
            
            if plotRfs
                ax = subplot(sply, splx, n);
                rg.plot.resizeobj(ax, 1.4);
                imagesc(real(h));
                caxis([-1 1]);
                axis off equal
            end
            
            D(:, n) = h(:) / 400;
        end
    end
end

drawnow();

%%

clear XFrames

neuronIsDead = false(nNeurons, 1);
neuronIsDead(deadNeurons) = true;

% Create arrays
t = (dt : dt : tTotal)';
tSignal = (dtSignal : dtSignal : tTotal)';
nT = numel(t);
nTSignal = numel(tSignal);
V = zeros(nNeurons, nT);
R = zeros(nNeurons, nT);
isSpikeAll = false(nNeurons, nT);

img = imread('~/Downloads/20170817_154831-1.jpg');
img = img(1:800, 300:1100, :);
img = imresize(img, 0.1);
img = squeeze(img(:, :, 1));
imgMean = mean(img(:));
X = zeros(nTSignal, imageWidth^2);
for s = 1:nTSignal
    tmp = imrotate(img, s * 360/nTSignal, 'crop');
    tmp = imresize(tmp, imageWidth*[1 1]);
    X(s, :) = tmp(:);
    XFrames(s, 1) = struct('cdata', uint8(tmp), 'colormap', gray(255));
end
X = (X-0.5*imgMean)/255;

% % Create 2D signal
% X = randn(nTSignal, imageWidth^2);
% for p1 = 1:imageWidth
%     for p2 = 1:imageWidth
%         idx = (p1-1)*imageWidth + p2;
%         X(:, idx) = rg.helpers.gsmooth(X(:, idx), [], signalNSmooth);
%     end
% end
% 
% h = fspecial('gaussian', [10 10], 0.7);
% for s = 1:nTSignal
%    tmp = X(s, :);
%    tmp = reshape(tmp, [imageWidth, imageWidth]);
%    tmp = imfilter(tmp, h);
%    tmp = tmp-min(tmp(:));
%    tmp = tmp/max(tmp(:));
%    X(s, :) = tmp(:);
%    XFrames(s, 1) = struct('cdata', uint8(tmp*255), 'colormap', gray(255));
% end
XInds = rg.helpers.binarySearch(t, tSignal);

% Compute recurrent weight matrix and thresholds
RW = [D'*D + mu*eye(nNeurons)];
vThresh = diag(RW)/2;
% DT = D';

tic();

for s = 2:10000
    
    xDot = (X(XInds(s), :) - X(XInds(s-1), :)) / dt;
    newSignal = ~all(xDot==0);
    c = xDot + X(XInds(s), :)*lambda;
    
    % Detect spikes
    isSpike = [V(:, s-1) >= vThresh];    
    isSpikeAll(:, s-1) = isSpike;
    noise = sigma*randn(nNeurons, 1)/sqrt(dt);
    
    % New method: summed columns
    term1 = -lambda*V(:, s-1);
%     term2 = DT*c';

    % Only recompute c*D if the signal has changed
    if newSignal || s==2
        term2 = (c*D)'; %bottleneck
    end
    
%     term3 = sum((RW(:, isSpike)), 2)/dt; %bottleneck
%     spikeInds = find(isSpike);
%     term3 = rg.helpers.matrixSum(RW, spikeInds, 2) / dt;
%     term3 = rg.helpers.matrixSumLog(RW, isSpike, 2) / dt;
    term3 = RW*isSpike / dt;
%     term3 = RW*isSpike / dt;
    dVdT = term1 + term2 + term3 + noise;
%     dVdT = -lambda*V(:, s-1) + D'*c' - sum((RW(:, isSpike))/dt, 2) + noise;

    % Old method: M-V multiplication for RW*isSpike
%     dVdT = -lambda*V(:, s-1) + D'*c' - (RW)*isSpike/dt + noise;
    dRdT = -lambda*R(:, s-1) + isSpike/dt;
    V(:, s) = V(:, s-1) + dt*dVdT;
    R(:, s) = R(:, s-1) + dt*dRdT;
    
end

toc

% Compute mean rate
meanRateCount = numel(find(isSpikeAll(:))) / nNeurons / tTotal;
meanRateCalc = mean(R(:));
fprintf('Mean spike rate was %.2f, calculated rate was %.2f\n', meanRateCount, meanRateCalc);

XHat = D*R;
clear XHatFrames
inds = rg.helpers.binarySearch(tSignal, t);
for s = 1:nTSignal
    idx = find(XInds==s, 1, 'last');
    tmp = reshape(XHat(:, idx), imageWidth*[1 1]);
    XHatFrames(s, 1) = struct('cdata', uint8(tmp*255 + 0.5*imgMean), 'colormap', gray(255));
end

fprintf('That took %.2f seconds\n', toc());

%% Plot an example neuron

figure();

neuronInds = [1001 1005];
colors = {'k', 'r'};

subplot(4, 4, 1);
img = reshape(D(:, neuronInds(1)), imageWidth*[1 1]);
imagesc(img);

subplot(4, 4, 2);
img = reshape(D(:, neuronInds(2)), imageWidth*[1 1]);
imagesc(img);

subplot(4, 1, 3);
y = V(neuronInds, :)';
line(t, y(:, 1), 'color', colors{1});
line(t, y(:, 2), 'color', colors{2});
title('Voltage');
hCursor(1) = line([0 0], get(gca, 'ylim'), 'color', 'b');
for n = 1:2
    y = vThresh(neuronInds(n));
    line(get(gca, 'xlim'), y*[1 1], 'color', colors{n}, 'lineStyle', '--');
end

subplot(4, 1, 4);
y = R(neuronInds, :)';
line(t, y(:, 1), 'color', colors{1});
line(t, y(:, 2), 'color', colors{2});
title('Rate');
for n = 1:2
    tSpike = find(isSpikeAll(neuronInds(n), :))*dt;
    x = [tSpike; tSpike; nan(size(tSpike))];
    y = [zeros(size(tSpike)); ones(size(tSpike)); nan(size(tSpike))];
    line(x(:), y(:), 'color', colors{n});
end
hCursor(2) = line([0 0], get(gca, 'ylim'), 'color', 'b');

clear ax
ax(1) = subplot(4, 4, 3);
ax(2) = subplot(4, 4, 4);
img1 = imagesc(ax(1), XFrames(1).cdata);
img2 = imagesc(ax(2), XHatFrames(1).cdata);
axis(ax(1), 'square');
axis(ax(2), 'square');
title(ax(1), 'X');
title(ax(2), 'XHat');
colormap(ax(1), gray());

colormap(ax(2), gray());
clims{1} = [80 255];
clims{2} = [120 230];

set(ax(1), 'clim', clims{1});
set(ax(2), 'clim', clims{2});
for f = 1:nTSignal
    img1.CData = XFrames(f).cdata();
    img2.CData = XHatFrames(f).cdata();
    hCursor(1).XData = tSignal(f)*[1 1];
    hCursor(2).XData = tSignal(f)*[1 1];
    drawnow();
    pause(5 * 1/30);
end

%%

figure();

subplot(1, 3, 1);
img = imresize(V, [1000, 1000]);
imagesc(t, 1:nNeurons, img);
axis xy
title('Voltages');
xlabel('Time / s');
ylabel('Neuron');


subplot(1, 3, 2);
img = imresize(R, [1000, 1000]);
imagesc(t, 1:nNeurons, img);
axis xy
title('Rates');
xlabel('Time / s');
ylabel('Neuron');

subplot(1, 3, 3);
y = mean(R, 2);
barh(y)
ylim([0 nNeurons]);
title('Rates');
xlabel('Mean rate');
ylabel('Neuron');

%% Spikes

figure();
img = double(isSpikeAll);
img = imresize(nNeurons, 2000);
imagesc(img);
colormap gray
set(gca, 'clim', [0 1]);


%% Plots

figure();

subplot(4, 2, 1);
imagesc(D);
xlabel('Pixel');
ylabel('Neuron');
axis equal
colormap(gca, rg.color.maps.seismic());
set(gca, 'clim', [-1 1]*.003);
title('Weights');

%%

% Signal
ax(1) = subplot(4, 2, 3);
ax(2) = subplot(4, 2, 4);

img1 = imagesc(ax(1), XFrames(1).cdata);
img2 = imagesc(ax(2), XHatFrames(1).cdata);
axis(ax(1), 'square');
axis(ax(2), 'square');
title(ax(1), 'X');
title(ax(2), 'XHat');
colormap(ax(1), gray());

colormap(ax(2), gray());
clim = [80 255];
set(ax, 'clim', clim);
for f = 1:nTSignal
    img1.CData = XFrames(f).cdata();
    img2.CData = XHatFrames(f).cdata();
    drawnow();
    pause(1/30);
end

%% Tuning curves

figure('colormap', hot());

suptitle('Tuning curves');

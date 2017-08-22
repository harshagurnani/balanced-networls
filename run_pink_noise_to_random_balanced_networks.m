function [] = run_pink_noise_to_random_balanced_networks( varargin )

params = parse_arguments( varargin );
% sIMULATION PARAMETERS

% Parameters to test over
meanN = 0;%[0 :0.5:2];
nParams = length(meanN);        %No. of parametric values

N = 60;                 % No. of neurons
% NX = 2;                 % Input dimensions
% mu = 0.02;              % Sparseness/Reset/regularity
lambda_d = 1;           % Timescale of firing rate decay
T = 100;                % ms, Simulation time
dt = 0.025;             %ms (by default)
NRepeats = 100;

% Sinusoidal input
freq = rand(NX,1);
phase = rand(NX,1)*2*pi;
input_x = 2*sin(freq*(1:T) + phase);

% input_x = randn(NX,T+100);
% for jj=1:NX
%     input_x(jj,:) = smooth()
% % plot(input_x');
% input_x = [sin(0.2*(1:T));sin(0.3*(1:T)+pi/5); sin(0.05*(1:T))] ;

% Optimise decoder weights
W = rand(NX, N/2)-0.5;       % Random decoding weights 
while rank(W)<NX
    W = rand(NX, N/2)-0.5;
end
W = [W -W; meanN*zeros(1,N)];

for jj=1:N
    W(:,jj) = W(:,jj)/norm(W(:,jj));
end
W = W/N*sqrt(NX);

% Calculate statistics
FR = nan(nParams, N);               %Firing rates
FanoF = FR;                         %Fano factors
COV = FR;
% SC = ones( nParams, N, N);           %Spike Cross-correlation
c=0;


for ampNoise = meanN
    c=c+1;
    [ Allspikes, Allrates, FullInput, read_out] = run_nD_random_network( N, W, mu, lambda_d,T, input_x, ampNoise, NRepeats );
    FR(c,:) = (mean(mean(Allspikes,3),2))' ;
    FanoF(c, :) = ( nanstd( (squeeze(sum(Allspikes,2)))' ).^2 ./  (nanmean( (squeeze(sum(Allspikes,2)))' ) +eps) )';
    for nrn=1:N
        tmp=nan(NRepeats,1);
        for trial =1:NRepeats
            tmp(trial) = find_COV(Allspikes(nrn,:,trial), dt);
        end
        COV(c, nrn) = nanmean(tmp);
        
    end
%     for nrn1 = 1:N-1
%         for nrn2=nrn1+1:N
%             tmp=0;
%             for trial = 1:NRepeats
%                 tmp = xcorr( Allrates(nrn1, :, trial), Allrates(nrn2,:,trial));
%                 SC(c, nrn1, nrn2) = mean( cov() ) corr
%             end
%         end
%     end
end
%% Plotting
% mf = mean(FanoF,2);
% figure;
% for jj=1:N
%     scatter(meanN, FanoF(:,jj));hold on
% end
% 
% scatter(meanN, mf', 'k', 'filled' )
% ylim([0 1.2])
end

function PStruct = parse_arguments()

end

function [cov] = find_COV( spikes, dt)
        cov = NaN;
        spiketimes = find(spikes==1)*dt;
        if ~isempty(spiketimes)
        isi = (spiketimes(2:end) - spiketimes(1:end-1));
        ub = max(prctile(isi,90),5);
        isi = isi(isi<=ub);
        cov = var(isi)/mean(isi);
        end
end

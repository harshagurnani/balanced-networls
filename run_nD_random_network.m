function [ AllSpikes, AllRates, X_in, EstimateX] = run_nD_random_network( N, W, mu, lambda_d,T, input_x, ampNoise, nRepeats, varargin )
%% Balanced network representing a NX input space. Driven by Input-dependent Excitation and Recurrent inhibition.
   % N - Number of neurons
   % W - Readout weight
   % mu - Sparsifying parameter
   % T - total sim time in simulation (ms)
   % input_x - input
   
   
   Tx = size(input_x,2);
   NX = size(input_x,1);
   
   
   tstep_x = T/Tx;   % Time step of x
   dt = 0.010;       % ms
   tbins  = floor(T/dt);
   
   
   %% Calculate inputs
   x=nan(NX+1, tbins);
   cn = dsp.ColoredNoise('SamplesPerFrame',Tx);
   
   for jj=1:NX
      x(jj,:) = interp1(0:tstep_x:T-tstep_x,  input_x(jj,:),   0:dt:T-dt,  'spline');    %Interpolated to integration time intervals
   end
   
%    freq = exp([log(0.01):log(1.2):log(25)]);     %kHz
% %    freq = 0.01:0.05:25;
% 
%    amps = 1./freq * max(freq);
%    x(NX+1,:) = zeros(1,tbins);
%    
   
   
   x_dot = (x(:,2:end)-x(:,1:end-1))/dt;                                   %First derivative of x
   
%    W = [W; 0.2*ones(1,N)];      % Inject equal shared noise in all neurons
   
   %% Parameters
   %W = Readout Weight matrix.  Size =  Nx x N
   Thresh = nan(N,1);   % Threshold for firing
   for nrn = 1:N
       Thresh(nrn) = mu/lambda_d^2/2 + W(:,nrn)'*W(:,nrn)/ 2;
   end
   rec_W = W'*W  +  mu*eye(N); % Recurrent weights
   rec_W = [rec_W zeros(N,1)];
%    if ~isempty(varargin)
%        killNrn = varargin{1};
%        killNrn(killNrn<1 | killNrn > N) = N+1;
%    else
%        killNrn = N+1;
%    end
   %% Initialise
%#####    Vm = nan(N, tbins);          % Membrane potential
%#####    Vm = nan(N,1);
   AllSpikes = nan(N, tbins, nRepeats);
   AllRates = nan(N, tbins, nRepeats);
   EstimateX = nan(NX, tbins, nRepeats);
   
   for trial = 1:nRepeats
       spikes = false(N, tbins);    % is spike?
       rates = zeros(N, tbins);     % Firing rates
%#####        Vm(:,1) = unifrnd(0, 1, [N,1]);
       Vm = unifrnd(0, 1, [N,1]);

       % Filter x to make it pink noise
    %    for ff = 1:size(freq,2)
    %    x(NX+1,:) = x(NX+1,:)+ amps(ff)*sin(2*pi*freq(ff)*(0:dt:T-dt)+rand*2*pi);
    %    end
    %    x(NX+1,:) = x(NX+1,:)/prctile(x(NX+1,:),95)*ampNoise;
    
       pinknoise = step(cn);
       x(NX+1,:) = ampNoise* interp1(0:tstep_x:T-tstep_x,  pinknoise,   0:dt:T-dt,  'spline');
       x_dot(NX+1,:) = (x(NX+1,2:end) - x(NX+1,1:end-1))/dt;
       
       %% Integrate
       for tt = 2:tbins
           input = W'*( x(:,tt) + x_dot(:,tt-1) );
           rec_syn = sum(rec_W(:,[spikes(:,tt-1); true]'),2);
    %        rec_syn = rec_W * spikes(:,tt-1);
%######            Vm(:,tt) = Vm(:,tt-1) + dt*(-lambda_d*Vm(:,tt-1)  + input) - rec_syn  ;
           Vm = Vm + dt*(-Vm + input) - rec_syn /lambda_d ;

           % Only 1 neuron spikes
%####            spks = find(Vm(:,tt) > Thresh(:));
           spks = find(Vm > Thresh(:));
           
           if ~isempty(spks)
             doSpike = spks(unidrnd(length(spks)));
    %          if ~any(doSpike == killNrn)
             spikes( doSpike, tt) = true; 
    %          end
           end
    %      spikes(Vm(:,tt) > Thresh(:),tt)=1;

           rates(:,tt) = rates(:,tt-1) + spikes(:,tt-1) - dt*lambda_d*rates(:,tt-1);
       end

       read_out = W*rates;

       AllSpikes(:,:, trial) = spikes;
       AllRates(:,:,trial) = rates;
       EstimateX(:,:, trial) = read_out(1:NX,:);
   
   end
   X_in = x(1:NX,:);
%    read_out_vec = reshape(read_out, [NX,NY,tbins]);
  
   %% Plotting
%    figure()
%    
%    subplot(2,2,1)
%    cmap=colormap(jet);
%    xc = floor(64/(NX+1));
%    xc=16;
%    for jj=1:min(3,NX)
%        
% %        plot(dt:dt:T-dt,  x_dot(jj,:),'color',cmap(xc*jj,:),'linestyle','--'); hold on
%        plot(0:dt:T-dt,  mean(EstimateX(jj,:,:),3)/prctile(mean(EstimateX(jj,:,:),3),95),'Linestyle','--','color',cmap(xc*jj,:)); hold on
%        plot(0:dt:T-dt,  x(jj,:)/prctile(x(jj,:),95),'color',cmap(xc*jj,:),'linewidth',2); hold on
%    end
%    xlim([0,T])
   
% %    plot_x = reshape(x, [NX,NY,tbins]);
% %    cin = [prctile(plot_x(:),5),prctile(plot_x(:),95)];
% %    cout = [min(prctile(read_out_vec(:,300:end),5)),max(prctile(read_out_vec(:,300:end),95))];
% % %    toplot = 
% %    
% %    for tt=1:Tx-1
% % %        tt
% %    subplot(3,3,1)
% %    colormap('jet')
% %    imagesc( plot_x(:,:,tt/dt),cin); 
% %    pause(0.0005)
% %    subplot(3,3,2)
% %    colormap('jet')
% %    imagesc(read_out_vec(:,:,tt/dt),cout);
% %    pause(0.0005)
% %    end
% %    title('Input x and Estimate')
% % %    legend(h, '$\hat{x}$', '$x$')
% %    set(legend,'Interpreter','latex')
% 
% 
%    
% %    subplot(3,3,2)
% %    for nrn=1:N
% %       scatter(W(1,nrn),W(2,nrn),'b');hold on 
% %    end
% %    for jj=1:length(killNrn)
% %       deadnrn = killNrn(jj);
% %       if deadnrn>0 && deadnrn<N+1
% %           scatter(W(1,deadnrn),W(2,deadnrn),'k','filled');hold on
% %       end
% %    end
% %    xlim([min(W(:))-0.2*abs(min(W(:))),max(W(:))*1.2])
% %    ylim([min(W(:))-0.2*abs(min(W(:))),max(W(:))*1.2])
% %    title('Read out weights')
%    
%    subplot(2,3,3)
%    imagesc(rec_W)
%    title('Recurrent weights')
%    
%    rate_m = max(rates(:));
%    subplot(2,2,3)
%    for nrn=1:N
%        plot(0:dt:T-dt,  rate_m*(nrn-1)+rates(nrn,:)); hold on
%    end
%    title('Firing rates')
%    xlim([0,T])
%    ylim([0,rate_m*nrn])
%    
%    subplot(2,2,4)
%    for nrn=1:N
%        scatter(dt*find(spikes(nrn,:)==1), nrn*ones(1,sum(spikes(nrn,:)==1)),4,'filled');hold on
%    end
%    xlim([0,T])
%    title('Spike raster')
%    
% %    subplot(3,1,3)
% %    cmap=colormap(jet);
% %    xc = floor(64/(NX+1));
% %    for jj=1:NX
% %        
% % %        plot(dt:dt:T-dt,  x_dot(jj,:),'color',cmap(xc*jj,:),'linestyle','--'); hold on
% %        plot(0:dt:T-dt,  read_out(jj,:),'Linestyle','--','color',cmap(xc*jj,:)); hold on
% %        plot(0:dt:T-dt,  x(jj,:),'color',cmap(xc*jj,:),'linewidth',2); hold on
% %    end
% %    xlim([0,T])
%    
% %    figure();
% %    nr = floor(sqrt(N));
% %    nc = nr+1;
% %    c=[min(TC(:)), max(TC(:))];
% %    for nrn=1:N
% %       subplot(nr,nc,nrn) 
% %       imagesc(TC(:,:,nrn));
% %    end
%    
% end
% 

function [ Vm,spikes] = run_nD_hetero_network( N, W, mu, T, input_x, varargin )
%% Balanced network representing a Nx-Variables. Driven by Input-dependent Excitation and Recurrent inhibition.
   % N - Number of neurons
   % W - Readout weight
   % mu - Sparsifying parameter
   % T - total sim time in simulation (ms)
   % input_x - input
   
   
   Tx = size(input_x,2);
   Nx = size(input_x,1);
   
   tstep_x = T/Tx;   % Time step of x
   dt = 0.025;       % ms
   tbins  = floor(T/dt);
   
   
   %% Calculate inputs
   x=nan(Nx, tbins);
   for jj=1:Nx
      x(jj,:) = interp1(0:tstep_x:T-tstep_x,  input_x(jj,:),   0:dt:T-dt,  'spline');    %Interpolated to integration time intervals
   end
   
   x_dot = (x(:,2:end)-x(:,1:end-1))/dt;                                   %First derivative of x
   
   
   %% Parameters
   %W = Readout Weight matrix.  Size =  Nx x N
   Thresh = nan(N,1);   % Threshold for firing
   lambda_d = 0.5;        % Exp kernel for rates
   for nrn = 1:N
       Thresh(nrn) = mu*lambda_d^2/2 + W(:,nrn)'*W(:,nrn)/ 2;
   end
   rec_W = W'*W  +  mu*eye(N)*lambda_d; % Recurrent weights
   
   if ~isempty(varargin)
       killNrn = varargin{1};
       killNrn(killNrn<1 | killNrn > N) = N+1;
   else
       killNrn = N+1;
   end
   %% Initialise
   Vm = nan(N, tbins);          % Membrane potential
   spikes = zeros(N, tbins);    % is spike?
   rates = zeros(N, tbins);     % Firing rates
   Vm(:,1) = unifrnd(0, 1, [N,1]);
   
   %% Integrate
   for tt = 2:tbins
       input = W'*( x(:,tt) + x_dot(:,tt-1) );
       rec_syn = rec_W*spikes(:,tt-1);
       
       Vm(:,tt) = Vm(:,tt-1) + dt*(-Vm(:,tt-1)  + input) - rec_syn  ;
       
       % Only 1 neuron spikes
       spks = find(Vm(:,tt) > Thresh(:));
       if ~isempty(spks)
         doSpike = spks(unidrnd(length(spks)));
         if ~any(doSpike == killNrn)
         spikes( doSpike, tt) = 1; 
         end
       end
%      spikes(Vm(:,tt) > Thresh(:),tt)=1;

       rates(:,tt) = rates(:,tt-1) + spikes(:,tt-1) - dt*lambda_d*rates(:,tt-1);
   end
   
   read_out = W*rates;
   
   %% Tuning curves
   xm = min(x')-eps; xM = max(x');
   nbins = 20; 
   dim1Bins = linspace(xm(1), xM(1), nbins+1);
   dim2Bins = linspace(xm(2), xM(2), nbins+1);
   TC = nan(nbins, nbins, N);
   
   for ii = 1:nbins
       for jj=1:nbins
           t_ind = find( x(1,:) > dim1Bins(ii) & x(1,:) <= dim1Bins(ii+1) & x(2,:) > dim2Bins(jj) & x(2,:) <= dim2Bins(jj+1) );
           for nrn = 1:N
               TC(ii,jj, nrn) = sum(spikes(nrn,t_ind))./(dt*length(t_ind)); 
           end
       end
   end
   
   %% Plotting
   figure()

   subplot(3,3,1)
   h(1)=plot(read_out(1,:),  read_out(2,:),'r--'); hold on
   h(2)=plot(x(1,:),  x(2,:),'k'); hold on
   title('Input x and Estimate')
   legend(h, '$\hat{x}$', '$x$')
   set(legend,'Interpreter','latex')


   
   subplot(3,3,2)
   for nrn=1:N
      scatter(W(1,nrn),W(2,nrn),'b');hold on 
   end
   for jj=1:length(killNrn)
      deadnrn = killNrn(jj);
      if deadnrn>0 && deadnrn<N+1
          scatter(W(1,deadnrn),W(2,deadnrn),'k','filled');hold on
      end
   end
   xlim([min(W(:))-0.2*abs(min(W(:))),max(W(:))*1.2])
   ylim([min(W(:))-0.2*abs(min(W(:))),max(W(:))*1.2])
   title('Read out weights')
   
   subplot(3,3,3)
   imagesc(rec_W)
   title('Recurrent weights')
   
   rate_m = max(rates(:));
   subplot(3,2,3)
   for nrn=1:N
       plot(0:dt:T-dt,  rate_m*(nrn-1)+rates(nrn,:)); hold on
   end
   title('Firing rates')
   xlim([0,T])
   ylim([0,rate_m*nrn])
   
   subplot(3,2,4)
   for nrn=1:N
       scatter(dt*find(spikes(nrn,:)==1), nrn*ones(1,sum(spikes(nrn,:)==1)),4,'filled');hold on
   end
   xlim([0,T])
   title('Spike raster')
   
   subplot(3,1,3)
   cmap=colormap(jet);
   xc = floor(64/(Nx+1));
   for jj=1:Nx
       
%        plot(dt:dt:T-dt,  x_dot(jj,:),'color',cmap(xc*jj,:),'linestyle','--'); hold on
       plot(0:dt:T-dt,  read_out(jj,:),'Linestyle','--','color',cmap(xc*jj,:)); hold on
       plot(0:dt:T-dt,  x(jj,:),'color',cmap(xc*jj,:),'linewidth',2); hold on
   end
   xlim([0,T])
   
   figure();
   nr = floor(sqrt(N));
   nc = nr+1;
   c=[min(TC(:)), max(TC(:))];
   for nrn=1:N
      subplot(nr,nc,nrn) 
      imagesc(TC(:,:,nrn));
   end
   
end


function [ Vm,spikes] = run_1D_hetero_network( N, lambda, mu, T, input_x, varargin )

   % N - Number of neurons
   % lambda - recurrent weight
   % t - total sim time in simulation (ms)
   % x - input
   
   tstep_x = T/length(input_x);   %Time step of x
   dt = 0.025;      %ms
   tbins  = floor(T/dt);
   
   
   x = interp1(0:tstep_x:T-tstep_x,  input_x,   0:dt:T-dt,  'spline');    %Interpolated to integration time intervals
   if isempty(varargin)
        x_dot = (x(2:end)-x(1:end-1))/dt;                                   %First derivative of x
   else
        x_dot = interp1(0:tstep_x:T-eps,  varargin{1},   0:dt:T-eps,  'spline'); 

   end
   
   if length(lambda)==1
       W = unifrnd(-lambda,lambda,[1,N]);  %Readout = Weight matrix  
   else
       W=lambda;
       lambda = mean(abs(lambda(:)));
   end
   Thresh = W*W' / 2 + mu/2;
   
   
   Vm = nan(N, tbins);
   spikes = zeros(N, tbins);
   
   % Initialise
   Vm(:,1) = unifrnd(0, 1, [N,1])*lambda^2;
   
   for tt = 2:tbins
       input = W'*( x(tt) + x_dot(tt-1) );
       rec_syn = (W'*W + mu*eye(N))*spikes(:,tt-1);
       
       Vm(:,tt) = Vm(:,tt-1) + dt*(-Vm(:,tt-1)  + input) - rec_syn  ;
       
       for nrn=1:N
           if Vm(nrn,tt) > Thresh
              spikes(nrn,tt) = 1;
           end
       end
   end
   
end


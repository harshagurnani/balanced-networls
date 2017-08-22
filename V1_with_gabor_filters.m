N2 =1200;        %Number of neurons
N = 2*N2;
jump = sqrt(N2/16)+1;
% nx=size(xxx,2)-1;
% ny=size(xxx,1)-1;
X = -1.1:0.1:1.1;   %Location x
Y = -1.1:0.1:1.1;   %Location y

nx = length(X);
ny = length(Y);

%% Gabor filters
[centres_x, centres_y] = meshgrid([-1.1:2.2/jump:1.1]', -1.1:2.2/jump:1.1);%[X(unidrnd( nx, 1,N)); Y(unidrnd(ny, 1,N))];
orientns = [ 0 : 2*pi/8 : 2*pi-eps];        % unifrnd(0,1,1,N)*2*pi;
W = nan(nx*ny,N);
phsshift = [0 1];
gabor_centres = nan(N,2);

nrn=0;cx = size(centres_x,1)^2;
for onoff=phsshift
for ori = orientns
for centre_idx = 1:cx
    nrn = nrn+1;
    if nrn<=N2
    % Module 1
    W(:,nrn) = reshape((1-2*onoff)*imag(make_gabor(X,Y, centres_x(centre_idx), centres_y(centre_idx), 0.3, ori, 0.7)), [nx*ny,1]);
    gabor_centres(nrn,:) = [    centres_x(centre_idx)    centres_y(centre_idx) ];
    
%     W(abs(W(:,nrn))<0.8*max(abs(W(:,nrn))), nrn)=0.05*sign(W(abs(W(:,nrn))<0.8*max(abs(W(:,nrn))), nrn));

    W(abs(W(:,nrn))<0.5*max(abs(W(:,nrn))), nrn)=0;
    W(:,nrn) = W(:,nrn)/max(abs(W(:,nrn)));
    
    % Module 2
    nrn2 = N2+nrn;
    W(:,nrn2) = reshape((1-2*onoff)*imag(make_gabor(X,Y, centres_x(centre_idx), centres_y(centre_idx), 0.15, ori, 0.3)), [nx*ny,1]);
    gabor_centres(nrn2,:) = [    centres_x(centre_idx)    centres_y(centre_idx) ];
    
%     W(abs(W(:,nrn2))<0.8*max(abs(W(:,nrn2))),nrn2)=0.05*sign(W(abs(W(:,nrn2))<0.8*max(abs(W(:,nrn2))), nrn2));
    W(abs(W(:,nrn2))<0.5*max(abs(W(:,nrn2))),nrn2)=0;
    W(:,nrn2) = W(:,nrn2)/max(abs(W(:,nrn2)));
    else
        break;
    end
end
end
end
% figure;scatter(gabor_centres(:,1), gabor_centres(:,2))
    % Restrict receptive fileds to maximally responsive regions
% 
% nr=10;nc=10;figure;
% for nrn=1:100
% subplot(nr,nc,nrn)
% imagesc(reshape(W(:,nrn),nx,ny),[-0.5 0.5])
% end
%% Make input in pixel space
% Tx = 1000; ts =40;
% input = randn(nx,ny, Tx);
% input = reshape(input,[nx*ny,Tx]);
% for pxl=1:nx*ny
%     input(pxl,:) = (smooth(input(pxl,:)',ts));
% end
% input = reshape(input(:,ts/2+1:Tx-ts/2), [nx, ny, Tx-ts]);
% [xx,  yy]= meshgrid(X,Y);
% gaus = exp(-((xx.^2)+(yy.^2))/0.03);
% for tt=1:Tx-ts
%     input(:,:,tt) = conv2(gaus,input(:,:,tt), 'same');
% end
% input = input/max(max(max(input)));

% rectangle
% input = zeros(nx,ny,1);
% input(14:17,4:7)=0.5;
% input = repmat( input ,[1, 1, Tx]);

% face
input = repmat(2*(pface-mean(pface(:)))/max(pface(:)), [1,1,Tx]);

% input = xxx';
% input = input-mean(input(:));
% input = repmat( input ,[1, 1, Tx]);

% figure;
% colormap('jet')
% for t=1:Tx
%     imagesc(input(:,:,t),[-0.5 0.5])
%     pause(0.01)
% end
%%
mu = 0.2; TTotal = 1000;lambda=1;
[ Vm,spikes, rates, read_out] = run_V1_network( N, W/N, mu,lambda, TTotal, 5*input);

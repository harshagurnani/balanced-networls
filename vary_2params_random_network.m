muR = [0:0.02:0.25];
inputN = [1 2 5 10 20 30];
nMU = length(muR);
nNx = length(inputN);

lgnds = cell(nNx,1);

% 
% figure;hold on;
% cmap = colormap(jet);
% jump = floor(64/(length(muR)+1));

b2 = 1;
mmf =nan(nMU,nNx);
mmfstd = mmf;

mcov = nan(nMU, nNx);
mcovstd = mcov;
for NX = inputN
    b=1;
for mu = muR
   
    PinkNoise_to_Random_Balanced_network
    nNoiseA = length(meanN);
    mf = nanmean(FanoF');
    mfstd = nanstd(FanoF');
    
    mc = nanmean(COV');
    mcstd = nanstd(COV');
    
    mmf(b,b2) = mean(mf);
    mmfstd(b,b2) = mean(mfstd);
    mcov(b,b2) = mean(mc);
    mcovstd(b,b2) = mean(mcstd);
    
%     h(b)=fill([1:nNoiseA, flip(1:nNoiseA)],[mf-mfstd, flip(mf+mfstd)], cmap(b*jump,: )); hold on
%     h(b).FaceAlpha = 0.2;
%     h2(b)=plot(1:nNoiseA, mf, 'color',cmap(b*jump,:), 'linewidth',2);
    
 b = b+1;
end


lgnds{b2} = ['NX = ' num2str(NX)];
b2=b2+1;
end

%% plotting
figure;
for b2=1:nNx
   h2(b2)=errorbar(muR',mmf(:,b2), mmfstd(:,b2), 'o-') ;hold on;
end
xlim([-0.02 0.22])
legend(lgnds)
xlabel('mu')
ylabel('Fano factor')
title('Reliability across trials')

figure;
for b2=1:nNx
h3(b2)=errorbar(muR',mcov(:,b2), mcovstd(:,b2), 'o-') ;hold on; 
end
xlim([-0.02 0.22])
legend(lgnds)
xlabel('mu')
ylabel('CoV')
title('Regularity of spike trains')

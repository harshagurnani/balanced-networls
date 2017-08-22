muR = [0:0.02:0.2];
nMU = length(muR);
lgnds = cell(nMU,1);


figure;hold on;
cmap = colormap(jet);
jump = floor(64/(length(muR)+1));
b=1;

mmf =nan(nMU,1);
mmfstd = mmf;
for mu = muR
   
    PinkNoise_to_Random_Balanced_network
    nNoiseA = length(meanN);
    mf = nanmean(FanoF');
    mfstd = nanstd(FanoF');
    
    mmf(b) = mean(mf);
    mmfstd(b) = mean(mfstd);
    
%     h(b)=fill([1:nNoiseA, flip(1:nNoiseA)],[mf-mfstd, flip(mf+mfstd)], cmap(b*jump,: )); hold on
%     h(b).FaceAlpha = 0.2;
%     h2(b)=plot(1:nNoiseA, mf, 'color',cmap(b*jump,:), 'linewidth',2);
    
     lgnds{b} = num2str(muR(b));
 b = b+1;
end
errorbar(muR',mmf, mmfstd, 'o') ;
xlim([-0.02 0.22])
% legend(h2, lgnds)
%% Genzel et al. (2019) rat+mouse data

alphas = 0:0.1:1;
betas = -1:0.1:1;

tabDIstaT = zeros(length(alphas),length(betas));
tabDIstaL = zeros(length(alphas),length(betas));
tabDIoveT = zeros(length(alphas),length(betas));
tabDIoveL = zeros(length(alphas),length(betas));
tabDIranT = zeros(length(alphas),length(betas));
tabDIranL = zeros(length(alphas),length(betas));
tabDItest = zeros(length(alphas),length(betas));
tabDIlast = zeros(length(alphas),length(betas));

for aaa=1:length(alphas)
    for bbb=1:length(betas)
        [alphas(aaa) betas(bbb)]
        results = launchSimILonObjectSpaceTask('mou', 1, 1000, alphas(aaa), betas(bbb), 0, 0);
        tabDIstaT(aaa,bbb) = results(1);
        tabDItest(aaa,bbb) = results(2);
        tabDIlast(aaa,bbb) = results(3);
        tabDIstaL(aaa,bbb) = results(4);
        tabDIoveT(aaa,bbb) = results(5);
        tabDIoveL(aaa,bbb) = results(6);
        tabDIranT(aaa,bbb) = results(7);
        tabDIranL(aaa,bbb) = results(8);
    end
end
theMax = max(max(max(tabDIoveL)),max(max(tabDIoveT)));
theMin = min(min(min(tabDIoveL)),min(min(tabDIoveT)));

syms alpha beta

figure
subplot(2,3,1)
imagesc(tabDIranL')
axis square
set(gca,'YDir','normal') % ,'XTick',[]
xticks([1 6 11])
xticklabels([0 0.5 1])
xlabel(texlabel(alpha))
yticks([1 11 21])
yticklabels([-1 0 1])
ylabel(texlabel(beta))
c = colorbar;
%c.Label.String = 'D.I.';
caxis([theMin theMax])
title('D.I. random')
ylabel('last spl')

subplot(2,3,4)
imagesc(tabDIranT')
axis square
set(gca,'YDir','normal') % ,'XTick',[]
xticks([1 6 11])
xticklabels([0 0.5 1])
xlabel(texlabel(alpha))
yticks([1 11 21])
yticklabels([-1 0 1])
ylabel(texlabel(beta))
c = colorbar;
%c.Label.String = 'D.I.';
caxis([theMin theMax])
ylabel('test')

subplot(2,3,2)
imagesc(tabDIstaL')
axis square
set(gca,'YDir','normal') % ,'XTick',[]
xticks([1 6 11])
xticklabels([0 0.5 1])
xlabel(texlabel(alpha))
yticks([1 11 21])
yticklabels([-1 0 1])
ylabel(texlabel(beta))
c = colorbar;
%c.Label.String = 'D.I.';
caxis([theMin theMax])
title('D.I. stable')

subplot(2,3,5)
imagesc(tabDIstaT')
axis square
set(gca,'YDir','normal') % ,'XTick',[]
xticks([1 6 11])
xticklabels([0 0.5 1])
xlabel(texlabel(alpha))
yticks([1 11 21])
yticklabels([-1 0 1])
ylabel(texlabel(beta))
c = colorbar;
%c.Label.String = 'D.I.';
caxis([theMin theMax])

subplot(2,3,3)
imagesc(tabDIoveL')
axis square
set(gca,'YDir','normal') % ,'XTick',[]
xticks([1 6 11])
xticklabels([0 0.5 1])
xlabel(texlabel(alpha))
yticks([1 11 21])
yticklabels([-1 0 1])
ylabel(texlabel(beta))
c = colorbar;
%c.Label.String = 'D.I.';
caxis([theMin theMax])
title('D.I. overlap')

subplot(2,3,6)
imagesc(tabDIoveT')
axis square
set(gca,'YDir','normal') % ,'XTick',[]
xticks([1 6 11])
xticklabels([0 0.5 1])
xlabel(texlabel(alpha))
yticks([1 11 21])
yticklabels([-1 0 1])
ylabel(texlabel(beta))
c = colorbar;
%c.Label.String = 'D.I.';
caxis([theMin theMax])

% colormap gray


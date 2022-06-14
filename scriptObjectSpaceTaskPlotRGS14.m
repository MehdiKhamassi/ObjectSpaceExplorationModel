%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2022 : Lobato et al., in prep. (RGS14/CONTROL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% scripts for comparing CON and RGS
tabParamScoreCON = tabParamScore;
tabParamScoreRGS = tabParamScore;

% comparing alphas (CON/RGS)
num2str([mean(tabParamScoreCON(:,1)) mean(tabParamScoreRGS(:,1))])
[p,tbl,stats] = kruskalwallis([tabParamScoreCON(:,1);tabParamScoreRGS(:,1)],[ones(size(tabParamScoreCON(:,1)));ones(size(tabParamScoreRGS(:,1)))*2],'off')

% figure alpha (all 3 conditions together)
figure
plot(1-0.1+rand(size(tabParamScoreCON(:,1)))/5,tabParamScoreCON(:,1),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
hold on
plot(2-0.1+rand(size(tabParamScoreRGS(:,1)))/5,tabParamScoreRGS(:,1),'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
alpha 0.5
boxplot([tabParamScoreCON(:,1);tabParamScoreRGS(:,1)],[ones(size(tabParamScoreCON(:,1)));ones(size(tabParamScoreRGS(:,1)))*2],'Colors','k','Symbol','')
plot([1 2],[1.1 1.1],'k')
text(1.4,1.15,'**')
%axis([0.5 2.5 -0.05 1.2])
syms alpha
ylabel(texlabel(alpha))
xlabel('all 3 conditions (except RANDOM) grouped together')
title('grey circle: Control; black triangle: RGS14')
clear alpha
set(gca, 'YScale', 'log')

outlierLimit = 10000;

% comparing betas (CON/RGS)
num2str([mean(tabParamScoreCON(abs(tabParamScoreCON(:,2))<outlierLimit,2)) mean(tabParamScoreRGS(abs(tabParamScoreRGS(:,2))<outlierLimit,2))])
[size(tabParamScoreCON(:,2)) size(tabParamScoreCON(abs(tabParamScoreCON(:,2))<outlierLimit,2)) size(tabParamScoreRGS(:,2)) size(tabParamScoreRGS(abs(tabParamScoreRGS(:,2))<outlierLimit,2))]
[p,tbl,stats] = kruskalwallis([tabParamScoreCON(abs(tabParamScoreCON(:,2))<outlierLimit,2);tabParamScoreRGS(abs(tabParamScoreRGS(:,2))<outlierLimit,2)],[ones(size(tabParamScoreCON(abs(tabParamScoreCON(:,2))<outlierLimit,1)));ones(size(tabParamScoreRGS(abs(tabParamScoreRGS(:,2))<outlierLimit,1)))*2],'off')
[p,tbl,stats] = kruskalwallis([tabParamScoreCON(:,2);tabParamScoreRGS(:,2)],[ones(size(tabParamScoreCON(:,1)));ones(size(tabParamScoreRGS(:,1)))*2],'off')
% comparing abs(betas) (CON/RGS)
[p,tbl,stats] = kruskalwallis(abs([tabParamScoreCON(:,2);tabParamScoreRGS(:,2)]),[ones(size(tabParamScoreCON(:,1)));ones(size(tabParamScoreRGS(:,1)))*2],'off')

% figure beta (all 3 conditions together)
figure
plot(1-0.1+rand(size(abs(tabParamScoreCON(:,2))))/5,abs(tabParamScoreCON(:,2)),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
hold on
plot(2-0.1+rand(size(abs(tabParamScoreRGS(:,2))))/5,abs(tabParamScoreRGS(:,2)),'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
alpha 0.5
boxplot(abs([tabParamScoreCON(:,2);tabParamScoreRGS(:,2)]),[ones(size(tabParamScoreCON(:,1)));ones(size(tabParamScoreRGS(:,1)))*2],'Colors','k','Symbol','')
axis([0.5 2.5 -5 30])
plot([1 2],[20.1 20.1],'k')
text(1.4,25,'n.s.')
syms beta
ylabel(texlabel(abs(beta)))
xlabel('all 3 conditions (except RANDOM) grouped together')
title('grey circle: Control; black triangle: RGS14')
clear beta

% comparing alphas in overlapping condition (3)
num2str([mean(tabParamScoreCON(tabParamScoreCON(:,end)==3,1)) mean(tabParamScoreRGS(tabParamScoreRGS(:,end)==3,1))])
[p,tbl,stats] = kruskalwallis([tabParamScoreCON(tabParamScoreCON(:,end)==3,1);tabParamScoreRGS(tabParamScoreRGS(:,end)==3,1)],[ones(size(tabParamScoreCON(tabParamScoreCON(:,end)==3,1)));ones(size(tabParamScoreRGS(tabParamScoreRGS(:,end)==3,1)))*2],'off')
% comparing betas in overlapping condition (3)
num2str([mean(tabParamScoreCON(tabParamScoreCON(:,end)==3&abs(tabParamScoreCON(:,2))<outlierLimit,2)) mean(tabParamScoreRGS(tabParamScoreRGS(:,end)==3&abs(tabParamScoreRGS(:,2))<outlierLimit,2))])
[size(tabParamScoreCON(tabParamScoreCON(:,end)==3,2)) size(tabParamScoreCON(tabParamScoreCON(:,end)==3&abs(tabParamScoreCON(:,2))<outlierLimit,2)) size(tabParamScoreRGS(tabParamScoreRGS(:,end)==3,2)) size(tabParamScoreRGS(tabParamScoreRGS(:,end)==3&abs(tabParamScoreRGS(:,2))<outlierLimit,2))]
[p,tbl,stats] = kruskalwallis([tabParamScoreCON(tabParamScoreCON(:,end)==3&abs(tabParamScoreCON(:,2))<outlierLimit,2);tabParamScoreRGS(tabParamScoreRGS(:,end)==3&abs(tabParamScoreRGS(:,2))<outlierLimit,2)],[ones(size(tabParamScoreCON(tabParamScoreCON(:,end)==3&abs(tabParamScoreCON(:,2))<outlierLimit,1)));ones(size(tabParamScoreRGS(tabParamScoreRGS(:,end)==3&abs(tabParamScoreRGS(:,2))<outlierLimit,1)))*2],'off')
[p,tbl,stats] = kruskalwallis([tabParamScoreCON(tabParamScoreCON(:,end)==3,2);tabParamScoreRGS(tabParamScoreRGS(:,end)==3,2)],[ones(size(tabParamScoreCON(tabParamScoreCON(:,end)==3,1)));ones(size(tabParamScoreRGS(tabParamScoreRGS(:,end)==3,1)))*2],'off')
% comparing abs(betas) in overlapping condition (3)
[p,tbl,stats] = kruskalwallis(abs([tabParamScoreCON(tabParamScoreCON(:,end)==3,2);tabParamScoreRGS(tabParamScoreRGS(:,end)==3,2)]),[ones(size(tabParamScoreCON(tabParamScoreCON(:,end)==3,1)));ones(size(tabParamScoreRGS(tabParamScoreRGS(:,end)==3,1)))*2],'off')

% comparing alphas in stable condition (2)
num2str([mean(tabParamScoreCON(tabParamScoreCON(:,end)==2,1)) mean(tabParamScoreRGS(tabParamScoreRGS(:,end)==2,1))])
[p,tbl,stats] = kruskalwallis([tabParamScoreCON(tabParamScoreCON(:,end)==2,1);tabParamScoreRGS(tabParamScoreRGS(:,end)==2,1)],[ones(size(tabParamScoreCON(tabParamScoreCON(:,end)==2,1)));ones(size(tabParamScoreRGS(tabParamScoreRGS(:,end)==2,1)))*2],'off')
% comparing betas in stable condition (2)
num2str([mean(tabParamScoreCON(tabParamScoreCON(:,end)==2&abs(tabParamScoreCON(:,2))<outlierLimit,2)) mean(tabParamScoreRGS(tabParamScoreRGS(:,end)==2&abs(tabParamScoreRGS(:,2))<outlierLimit,2))])
[size(tabParamScoreCON(tabParamScoreCON(:,end)==2,2)) size(tabParamScoreCON(tabParamScoreCON(:,end)==2&abs(tabParamScoreCON(:,2))<outlierLimit,2)) size(tabParamScoreRGS(tabParamScoreRGS(:,end)==2,2)) size(tabParamScoreRGS(tabParamScoreRGS(:,end)==2&abs(tabParamScoreRGS(:,2))<outlierLimit,2))]
[p,tbl,stats] = kruskalwallis([tabParamScoreCON(tabParamScoreCON(:,end)==2&abs(tabParamScoreCON(:,2))<outlierLimit,2);tabParamScoreRGS(tabParamScoreRGS(:,end)==2&abs(tabParamScoreRGS(:,2))<outlierLimit,2)],[ones(size(tabParamScoreCON(tabParamScoreCON(:,end)==2&abs(tabParamScoreCON(:,2))<outlierLimit,1)));ones(size(tabParamScoreRGS(tabParamScoreRGS(:,end)==2&abs(tabParamScoreRGS(:,2))<outlierLimit,1)))*2],'off')
[p,tbl,stats] = kruskalwallis([tabParamScoreCON(tabParamScoreCON(:,end)==2,2);tabParamScoreRGS(tabParamScoreRGS(:,end)==2,2)],[ones(size(tabParamScoreCON(tabParamScoreCON(:,end)==2,1)));ones(size(tabParamScoreRGS(tabParamScoreRGS(:,end)==2,1)))*2],'off')
% comparing abs(betas) in stable condition (2)
[p,tbl,stats] = kruskalwallis(abs([tabParamScoreCON(tabParamScoreCON(:,end)==2,2);tabParamScoreRGS(tabParamScoreRGS(:,end)==2,2)]),[ones(size(tabParamScoreCON(tabParamScoreCON(:,end)==2,1)));ones(size(tabParamScoreRGS(tabParamScoreRGS(:,end)==2,1)))*2],'off')

% comparing alphas in random condition (1)
num2str([mean(tabParamScoreCON(tabParamScoreCON(:,end)==1,1)) mean(tabParamScoreRGS(tabParamScoreRGS(:,end)==1,1))])
[p,tbl,stats] = kruskalwallis([tabParamScoreCON(tabParamScoreCON(:,end)==1,1);tabParamScoreRGS(tabParamScoreRGS(:,end)==1,1)],[ones(size(tabParamScoreCON(tabParamScoreCON(:,end)==1,1)));ones(size(tabParamScoreRGS(tabParamScoreRGS(:,end)==1,1)))*2],'off')
% comparing betas in random condition (1)
num2str([mean(tabParamScoreCON(tabParamScoreCON(:,end)==1&abs(tabParamScoreCON(:,2))<outlierLimit,2)) mean(tabParamScoreRGS(tabParamScoreRGS(:,end)==1&abs(tabParamScoreRGS(:,2))<outlierLimit,2))])
[size(tabParamScoreCON(tabParamScoreCON(:,end)==1,2)) size(tabParamScoreCON(tabParamScoreCON(:,end)==1&abs(tabParamScoreCON(:,2))<outlierLimit,2)) size(tabParamScoreRGS(tabParamScoreRGS(:,end)==1,2)) size(tabParamScoreRGS(tabParamScoreRGS(:,end)==1&abs(tabParamScoreRGS(:,2))<outlierLimit,2))]
[p,tbl,stats] = kruskalwallis([tabParamScoreCON(tabParamScoreCON(:,end)==1,2);tabParamScoreRGS(tabParamScoreRGS(:,end)==1,2)],[ones(size(tabParamScoreCON(tabParamScoreCON(:,end)==1,1)));ones(size(tabParamScoreRGS(tabParamScoreRGS(:,end)==1,1)))*2],'off')
[p,tbl,stats] = kruskalwallis([tabParamScoreCON(tabParamScoreCON(:,end)==1&abs(tabParamScoreCON(:,2))<outlierLimit,2);tabParamScoreRGS(tabParamScoreRGS(:,end)==1&abs(tabParamScoreRGS(:,2))<outlierLimit,2)],[ones(size(tabParamScoreCON(tabParamScoreCON(:,end)==1&abs(tabParamScoreCON(:,2))<outlierLimit,1)));ones(size(tabParamScoreRGS(tabParamScoreRGS(:,end)==1&abs(tabParamScoreRGS(:,2))<outlierLimit,1)))*2],'off')
% comparing abs(betas) in random condition (1)
[p,tbl,stats] = kruskalwallis(abs([tabParamScoreCON(tabParamScoreCON(:,end)==1,2);tabParamScoreRGS(tabParamScoreRGS(:,end)==1,2)]),[ones(size(tabParamScoreCON(tabParamScoreCON(:,end)==1,1)));ones(size(tabParamScoreRGS(tabParamScoreRGS(:,end)==1,1)))*2],'off')

% figure alpha overlapping vs. stable vs. random
figure
plot(1-0.1+rand(size(tabParamScoreCON(tabParamScoreCON(:,end)==3,1)))/5,tabParamScoreCON(tabParamScoreCON(:,end)==3,1),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
hold on
plot(2-0.1+rand(size(tabParamScoreRGS(tabParamScoreRGS(:,end)==3,1)))/5,tabParamScoreRGS(tabParamScoreRGS(:,end)==3,1),'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
plot(3-0.1+rand(size(tabParamScoreCON(tabParamScoreCON(:,end)==2,1)))/5,tabParamScoreCON(tabParamScoreCON(:,end)==2,1),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
plot(4-0.1+rand(size(tabParamScoreRGS(tabParamScoreRGS(:,end)==2,1)))/5,tabParamScoreRGS(tabParamScoreRGS(:,end)==2,1),'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
plot(5-0.1+rand(size(tabParamScoreCON(tabParamScoreCON(:,end)==1,1)))/5,tabParamScoreCON(tabParamScoreCON(:,end)==1,1),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
plot(6-0.1+rand(size(tabParamScoreRGS(tabParamScoreRGS(:,end)==1,1)))/5,tabParamScoreRGS(tabParamScoreRGS(:,end)==1,1),'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
alpha 0.5
boxplot([tabParamScoreCON(tabParamScoreCON(:,end)==3,1);tabParamScoreRGS(tabParamScoreRGS(:,end)==3,1);tabParamScoreCON(tabParamScoreCON(:,end)==2,1);tabParamScoreRGS(tabParamScoreRGS(:,end)==2,1);tabParamScoreCON(tabParamScoreCON(:,end)==1,1);tabParamScoreRGS(tabParamScoreRGS(:,end)==1,1)],[ones(size(tabParamScoreCON(tabParamScoreCON(:,end)==3,1)));ones(size(tabParamScoreRGS(tabParamScoreRGS(:,end)==3,1)))*2;ones(size(tabParamScoreCON(tabParamScoreCON(:,end)==2,1)))*3;ones(size(tabParamScoreRGS(tabParamScoreRGS(:,end)==2,1)))*4;ones(size(tabParamScoreCON(tabParamScoreCON(:,end)==1,1)))*5;ones(size(tabParamScoreRGS(tabParamScoreRGS(:,end)==1,1)))*6],'Colors','k','Symbol','')
% no outliers
plot([1 2],[1.1 1.1],'k')
plot([3 4],[1.1 1.1],'k')
plot([5 6],[1.1 1.1],'k')
text(1.4,1.2,'n.s.') % text(1.4,1.2,'*')
text(3.4,1.2,'n.s.')
text(5.4,1.2,'n.s.')
axis([0.5 6.5 -0.05 1.3])
syms alpha
ylabel(texlabel(alpha))
xlabel('overlapping                         stable                         random')
title('grey circle: Control; black triangle: RGS14')
clear alpha

outlierLimit = 20; % 2000 

% figure abs(beta) overlapping vs. stable vs. random
figure
plot(1-0.1+rand(size(abs(tabParamScoreCON(tabParamScoreCON(:,end)==3,2))))/5,abs(tabParamScoreCON(tabParamScoreCON(:,end)==3,2)),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
hold on
plot(2-0.1+rand(size(abs(tabParamScoreRGS(tabParamScoreRGS(:,end)==3,2))))/5,abs(tabParamScoreRGS(tabParamScoreRGS(:,end)==3,2)),'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
plot(3-0.1+rand(size(abs(tabParamScoreCON(tabParamScoreCON(:,end)==2,2))))/5,abs(tabParamScoreCON(tabParamScoreCON(:,end)==2,2)),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
plot(4-0.1+rand(size(abs(tabParamScoreRGS(tabParamScoreRGS(:,end)==2,2))))/5,abs(tabParamScoreRGS(tabParamScoreRGS(:,end)==2,2)),'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
plot(5-0.1+rand(size(abs(tabParamScoreCON(tabParamScoreCON(:,end)==1,2))))/5,abs(tabParamScoreCON(tabParamScoreCON(:,end)==1,2)),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
plot(6-0.1+rand(size(abs(tabParamScoreRGS(tabParamScoreRGS(:,end)==1,2))))/5,abs(tabParamScoreRGS(tabParamScoreRGS(:,end)==1,2)),'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
alpha 0.5
boxplot(abs([tabParamScoreCON(tabParamScoreCON(:,end)==3&abs(tabParamScoreCON(:,2))<outlierLimit,2);tabParamScoreRGS(tabParamScoreRGS(:,end)==3&abs(tabParamScoreRGS(:,2))<outlierLimit,2);tabParamScoreCON(tabParamScoreCON(:,end)==2&abs(tabParamScoreCON(:,2))<outlierLimit,2);tabParamScoreRGS(tabParamScoreRGS(:,end)==2&abs(tabParamScoreRGS(:,2))<outlierLimit,2);tabParamScoreCON(tabParamScoreCON(:,end)==1&abs(tabParamScoreCON(:,2))<outlierLimit,2);tabParamScoreRGS(tabParamScoreRGS(:,end)==1&abs(tabParamScoreRGS(:,2))<outlierLimit,2)]),[ones(size(tabParamScoreCON(tabParamScoreCON(:,end)==3&abs(tabParamScoreCON(:,2))<outlierLimit,1)));ones(size(tabParamScoreRGS(tabParamScoreRGS(:,end)==3&abs(tabParamScoreRGS(:,2))<outlierLimit,1)))*2;ones(size(tabParamScoreCON(tabParamScoreCON(:,end)==2&abs(tabParamScoreCON(:,2))<outlierLimit,1)))*3;ones(size(tabParamScoreRGS(tabParamScoreRGS(:,end)==2&abs(tabParamScoreRGS(:,2))<outlierLimit,1)))*4;ones(size(tabParamScoreCON(tabParamScoreCON(:,end)==1&abs(tabParamScoreCON(:,2))<outlierLimit,1)))*5;ones(size(tabParamScoreRGS(tabParamScoreRGS(:,end)==1&abs(tabParamScoreRGS(:,2))<outlierLimit,1)))*6],'Colors','k','Symbol','')
%% if outlierLimit == 2000 M1
% plot([0 7],[1400 1400],'k')
% plot([0 7],[1390 1390],'k')
% % outliers : 0 5 0 0 5 2
% plot([3 3 3 3],[1420 1440 1460 1480],'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
% plot([4 4 4 4],[1420 1440 1460 1480],'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
% plot([5 5],[1420 1440],'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
% plot([6 6 6],[1420 1440 1460],'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
% plot([1 2],[1500 1500],'k')
% plot([3 4],[1500 1500],'k')
% plot([5 6],[1500 1500],'k')
% text(1.4,1525,'n.s.')
% text(3.4,1525,'n.s.')
% text(5.4,1525,'n.s.')
% axis([0.5 6.5 -5 1550])
%% if outlierLimit == 20 M1
% plot([0 7],[16.25 16.25],'k')
% plot([0 7],[16.75 16.75],'k')
% % outliers : 0 5 0 0 5 2
% plot([1]-0.1+rand(1,1)/5,[18]-0.1+rand(1,1),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
% plot([2 2]-0.1+rand(1,2)/5,[18 18.1]-0.1+rand(1,2),'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
% plot([3 3 3 3 3 3 3 3 3 3 3]-0.1+rand(1,11)/5,[17.6 17.7 17.8 17.9 18 18.1 18.2 18.3 18.4 18.5 18.6]-0.1+rand(1,11),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
% plot([4 4 4 4 4 4 4 4 4 4 4]-0.1+rand(1,11)/5,[17.6 17.7 17.8 17.9 18 18.1 18.2 18.3 18.4 18.5 18.6]-0.1+rand(1,11),'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
% plot([5 5 5]-0.1+rand(1,3)/5,[18 18.1 18.2]-0.1+rand(1,3),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
% plot([6 6]-0.1+rand(1,2)/5,[18 18.1]-0.1+rand(1,2),'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
% plot([1 2],[20 20],'k')
% plot([3 4],[20 20],'k')
% plot([5 6],[20 20],'k')
% text(1.4,21,'n.s.')
% text(3.4,21,'n.s.')
% text(5.4,21,'n.s.')
% axis([0.5 6.5 -1 23])
%% if outlierLimit == 20 M4 = 3 3 4 3 1 3
plot([0 7],[20.5 20.5],'k')
plot([0 7],[21 21],'k')
% outliers : 0 5 0 0 5 2
plot([1 1 1]-0.1+rand(1,3)/5,[21.5 21.6 21.7]-0.1+rand(1,3),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
plot([2 2 2]-0.1+rand(1,3)/5,[21.5 21.6 21.7]-0.1+rand(1,3),'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
plot([3 3 3 3]-0.1+rand(1,4)/5,[21.5 21.6 21.7 21.8]-0.1+rand(1,4),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
plot([4 4 4]-0.1+rand(1,3)/5,[21.5 21.6 21.7]-0.1+rand(1,3),'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
plot([5]-0.1+rand(1,1)/5,[21.5]-0.1+rand(1,1),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
plot([6 6 6]-0.1+rand(1,3)/5,[21.5 21.6 21.7]-0.1+rand(1,3),'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
plot([1 2],[23 23],'k')
plot([3 4],[23 23],'k')
plot([5 6],[23 23],'k')
text(1.4,23.5,'n.s.')
text(3.4,23.5,'n.s.')
text(5.4,23.5,'n.s.')
axis([0.5 6.5 -1 24])
syms beta
ylabel(texlabel(abs(beta)))
xlabel('overlapping                         stable                         random')
title('grey circle: Control; black triangle: RGS14')

%% figure 2D alpha x beta, C24 x R24, clusters
figure
theCond = 1; % 1 random ; 2 stable ; 3 overlapping
outlierLimit = 12;
theMax = outlierLimit; % 25 for rats 1.5 for mice
offset = theMax/3; space = theMax/40;
plot([-0.05 1.05],[0 0],'--','Color',[0.5 0.5 0.5])
hold on
plot([0.5 0.5],[-theMax theMax],'--','Color',[0.5 0.5 0.5])
plot(tabParamScoreRGS(tabParamScoreRGS(:,end)==theCond,1),tabParamScoreRGS(tabParamScoreRGS(:,end)==theCond,2),'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
plot(tabParamScoreCON(tabParamScoreCON(:,end)==theCond,1),tabParamScoreCON(tabParamScoreCON(:,end)==theCond,2),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
syms alpha beta
xlabel(texlabel(alpha))
ylabel(texlabel(beta))
switch (theCond)
    case 1
        title('grey circle: Control; black triangle: RGS14 (RANDOM)')
    case 2
        title('grey circle: Control; black triangle: RGS14 (STABLE)')
    case 3
        title('grey circle: Control; black triangle: RGS14 (OVERLAPPING)')
end
%title(['Fit of Model ' num2str(whichModel) ' on mouse data OVERLAPPING session'])
%axis([0 1 -0.8 1.7])
axis([-0.05 1.05 -theMax theMax])

% adding the outliers
plot([-0.05 1.05],[theMax-offset+space/2 theMax-offset+space/2],'-k')
plot([-0.05 1.05],[theMax-offset-space/2 theMax-offset-space/2],'-k')
plot([-0.05 1.05],[-theMax+offset+space/2 -theMax+offset+space/2],'-k')
plot([-0.05 1.05],[-theMax+offset-space/2 -theMax+offset-space/2],'-k')
outliersCON = tabParamScoreCON(tabParamScoreCON(:,end)==theCond&abs(tabParamScoreCON(:,2))>theMax-offset+space/2,1:2);
outliersRGS = tabParamScoreRGS(tabParamScoreRGS(:,end)==theCond&abs(tabParamScoreRGS(:,2))>theMax-offset+space/2,1:2);
outliersPosCON = outliersCON(outliersCON(:,2)>=0,:);
outliersPosRGS = outliersRGS(outliersRGS(:,2)>=0,:);
[boubouCON, indCON] = sort(outliersPosCON(:,2),1,'descend');
[boubouRGS, indRGS] = sort(outliersPosRGS(:,2),1,'descend');
outliersPosCON = outliersPosCON(indCON,:);
outliersPosRGS = outliersPosRGS(indRGS,:);
outliersNegCON = outliersCON(outliersCON(:,2)<0,:);
outliersNegRGS = outliersRGS(outliersRGS(:,2)<0,:);
[boubouCON, indCON] = sort(outliersNegCON(:,2),1,'descend');
[boubouRGS, indRGS] = sort(outliersNegRGS(:,2),1,'descend');
outliersNegCON = outliersNegCON(indCON,:);
outliersNegRGS = outliersNegRGS(indRGS,:);
nbOP = size(outliersPosCON,1);
for iii=1:nbOP
    if (outliersPosCON(iii,1)>=0.5)
        plot(outliersPosCON(iii,1),theMax-(offset-space/2)*iii/(nbOP+1),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
    else
        plot(outliersPosCON(iii,1),theMax-(offset-space/2)*iii/(nbOP+1),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
    end
end
nbOP = size(outliersPosRGS,1);
for iii=1:nbOP
    if (outliersPosRGS(iii,1)>=0.5)
        plot(outliersPosRGS(iii,1),theMax-(offset-space/2)*iii/(nbOP+1),'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
    else
        plot(outliersPosRGS(iii,1),theMax-(offset-space/2)*iii/(nbOP+1),'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
    end
end
nbON = size(outliersNegCON,1);
for iii=1:nbON
    if (outliersNegCON(iii,1)>=0.5)
        plot(outliersNegCON(iii,1),-theMax+(offset+space/2)*iii/(nbON+1),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
    else
        plot(outliersNegCON(iii,1),-theMax+(offset+space/2)*iii/(nbON+1),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
    end
end
nbON = size(outliersNegRGS,1);
for iii=1:nbON
    if (outliersNegRGS(iii,1)>=0.5)
        plot(outliersNegRGS(iii,1),-theMax+(offset+space/2)*iii/(nbON+1),'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
    else
        plot(outliersNegRGS(iii,1),-theMax+(offset+space/2)*iii/(nbON+1),'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
    end
end

%% kmeans_md K-means clustering
%[c s] = kmeans_md([tabParamScoreCON(tabParamScoreCON(:,end)==theCond,1:2);tabParamScoreRGS(tabParamScoreRGS(:,end)==theCond,1:2)],2)
% [c_a s_a] = kmeans([tabParamScoreCON(tabParamScoreCON(:,end)==theCond,1);tabParamScoreRGS(tabParamScoreRGS(:,end)==theCond,1)],2)
% [c_b s_b] = kmeans([tabParamScoreCON(tabParamScoreCON(:,end)==theCond,2);tabParamScoreRGS(tabParamScoreRGS(:,end)==theCond,2)],2)
[c_aCON s_aCON] = kmeans(tabParamScoreCON(tabParamScoreCON(:,end)==theCond,1)',2);
[c_bCON s_bCON] = kmeans(tabParamScoreCON(tabParamScoreCON(:,end)==theCond&abs(tabParamScoreCON(:,2))<=outlierLimit,2)',2);
[c_aRGS s_aRGS] = kmeans(tabParamScoreRGS(tabParamScoreRGS(:,end)==theCond,1)',2);
[c_bRGS s_bRGS] = kmeans(tabParamScoreRGS(tabParamScoreRGS(:,end)==theCond&abs(tabParamScoreRGS(:,2))<=outlierLimit,2)',2);
kmeans_alpha = [c_aCON ; c_aRGS]
kmeans_beta = num2str([c_bCON ; c_bRGS])

%% plot ellipses
eccentricity = 0.5;
numPoints = 300; % Less for a coarser ellipse, more for a finer resolution.
% Define parameters.
x1 = c_aCON(1)-0.2;
x2 = c_aCON(1)+0.2;
y1 = c_bCON(1);
y2 = c_bCON(1);
% Make equations:
a = (1/2) * sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2);
b = a * sqrt(1-eccentricity^2);
t = linspace(0, 2 * pi, numPoints); % Absolute angle parameter
X = a * cos(t);
Y = b * sin(t);
% Compute angles relative to (x1, y1).
angles = atan2(y2 - y1, x2 - x1);
x = (x1 + x2) / 2 + X * cos(angles) - Y * sin(angles);
y = (y1 + y2) / 2 + X * sin(angles) + Y * cos(angles);
% Plot the ellipse as a blue curve.
plot(x,y,'-', 'Color',[0.5 0.5 0.5], 'LineWidth', 1);	% Plot ellipse
% Define parameters.
x1 = c_aCON(2)-0.2;
x2 = c_aCON(2)+0.2;
y1 = c_bCON(2);
y2 = c_bCON(2);
% Make equations:
a = (1/2) * sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2);
b = a * sqrt(1-eccentricity^2);
t = linspace(0, 2 * pi, numPoints); % Absolute angle parameter
X = a * cos(t);
Y = b * sin(t);
% Compute angles relative to (x1, y1).
angles = atan2(y2 - y1, x2 - x1);
x = (x1 + x2) / 2 + X * cos(angles) - Y * sin(angles);
y = (y1 + y2) / 2 + X * sin(angles) + Y * cos(angles);
% Plot the ellipse as a blue curve.
plot(x,y,'-', 'Color',[0.5 0.5 0.5], 'LineWidth', 1);	% Plot ellipse
% Define parameters.
x1 = c_aRGS(1)-0.2;
x2 = c_aRGS(1)+0.2;
y1 = c_bRGS(1);
y2 = c_bRGS(1);
% Make equations:
a = (1/2) * sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2);
b = a * sqrt(1-eccentricity^2);
t = linspace(0, 2 * pi, numPoints); % Absolute angle parameter
X = a * cos(t);
Y = b * sin(t);
% Compute angles relative to (x1, y1).
angles = atan2(y2 - y1, x2 - x1);
x = (x1 + x2) / 2 + X * cos(angles) - Y * sin(angles);
y = (y1 + y2) / 2 + X * sin(angles) + Y * cos(angles);
% Plot the ellipse as a blue curve.
plot(x,y,'-k', 'LineWidth', 1);	% Plot ellipse
% Define parameters.
x1 = c_aRGS(2)-0.2;
x2 = c_aRGS(2)+0.2;
y1 = c_bRGS(2);
y2 = c_bRGS(2);
% Make equations:
a = (1/2) * sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2);
b = a * sqrt(1-eccentricity^2);
t = linspace(0, 2 * pi, numPoints); % Absolute angle parameter
X = a * cos(t);
Y = b * sin(t);
% Compute angles relative to (x1, y1).
angles = atan2(y2 - y1, x2 - x1);
x = (x1 + x2) / 2 + X * cos(angles) - Y * sin(angles);
y = (y1 + y2) / 2 + X * sin(angles) + Y * cos(angles);
% Plot the ellipse as a blue curve.
plot(x,y,'-k', 'LineWidth', 1);	% Plot ellipse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2019 : Schut et al 2020 (WT/KO)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% scripts for comparing WT and KO
tabParamScoreWT = tabParamScore;
tabParamScoreKO = tabParamScore;

% comparing alphas (WT/KO)
num2str([mean(tabParamScoreWT(:,1)) mean(tabParamScoreKO(:,1))])
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(:,1);tabParamScoreKO(:,1)],[ones(size(tabParamScoreWT(:,1)));ones(size(tabParamScoreKO(:,1)))*2],'off')

outlierLimit = 10000;

% comparing betas (WT/KO)
num2str([mean(tabParamScoreWT(abs(tabParamScoreWT(:,2))<outlierLimit,2)) mean(tabParamScoreKO(abs(tabParamScoreKO(:,2))<outlierLimit,2))])
[size(tabParamScoreWT(:,2)) size(tabParamScoreWT(abs(tabParamScoreWT(:,2))<outlierLimit,2)) size(tabParamScoreKO(:,2)) size(tabParamScoreKO(abs(tabParamScoreKO(:,2))<outlierLimit,2))]
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(abs(tabParamScoreWT(:,2))<outlierLimit,2);tabParamScoreKO(abs(tabParamScoreKO(:,2))<outlierLimit,2)],[ones(size(tabParamScoreWT(abs(tabParamScoreWT(:,2))<outlierLimit,1)));ones(size(tabParamScoreKO(abs(tabParamScoreKO(:,2))<outlierLimit,1)))*2],'off')
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(:,2);tabParamScoreKO(:,2)],[ones(size(tabParamScoreWT(:,1)));ones(size(tabParamScoreKO(:,1)))*2],'off')

% comparing alphas in overlapping condition (3)
num2str([mean(tabParamScoreWT(tabParamScoreWT(:,end)==3,1)) mean(tabParamScoreKO(tabParamScoreKO(:,end)==3,1))])
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(tabParamScoreWT(:,end)==3,1);tabParamScoreKO(tabParamScoreKO(:,end)==3,1)],[ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==3,1)));ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==3,1)))*2],'off')
% comparing betas in overlapping condition (3)
num2str([mean(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit,2)) mean(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit,2))])
[size(tabParamScoreWT(tabParamScoreWT(:,end)==3,2)) size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit,2)) size(tabParamScoreKO(tabParamScoreKO(:,end)==3,2)) size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit,2))]
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit,2);tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit,2)],[ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit,1)));ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit,1)))*2],'off')
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(tabParamScoreWT(:,end)==3,2);tabParamScoreKO(tabParamScoreKO(:,end)==3,2)],[ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==3,1)));ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==3,1)))*2],'off')
% for the figure of abs(beta) overlapping during bins of 5 trials
compteur = 1;
vectAlpha = []; % alpha
vectData = []; % beta
vectName = [];
% then after extracting data for each bin of 5 trials:
vectData = [vectData ; tabParamScoreWT(tabParamScoreWT(:,end)==3,2);tabParamScoreKO(tabParamScoreKO(:,end)==3,2)];
vectAlpha = [vectAlpha ; tabParamScoreWT(tabParamScoreWT(:,end)==3,1);tabParamScoreKO(tabParamScoreKO(:,end)==3,1)];
vectName = [vectName ; (ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==3,1)))*compteur) ; (ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==3,1)))*(compteur+1))];
compteur = compteur + 2;

% comparing alphas in stable condition (2)
num2str([mean(tabParamScoreWT(tabParamScoreWT(:,end)==2,1)) mean(tabParamScoreKO(tabParamScoreKO(:,end)==2,1))])
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(tabParamScoreWT(:,end)==2,1);tabParamScoreKO(tabParamScoreKO(:,end)==2,1)],[ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==2,1)));ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==2,1)))*2],'off')
% comparing betas in stable condition (2)
num2str([mean(tabParamScoreWT(tabParamScoreWT(:,end)==2&abs(tabParamScoreWT(:,2))<outlierLimit,2)) mean(tabParamScoreKO(tabParamScoreKO(:,end)==2&abs(tabParamScoreKO(:,2))<outlierLimit,2))])
[size(tabParamScoreWT(tabParamScoreWT(:,end)==2,2)) size(tabParamScoreWT(tabParamScoreWT(:,end)==2&abs(tabParamScoreWT(:,2))<outlierLimit,2)) size(tabParamScoreKO(tabParamScoreKO(:,end)==2,2)) size(tabParamScoreKO(tabParamScoreKO(:,end)==2&abs(tabParamScoreKO(:,2))<outlierLimit,2))]
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(tabParamScoreWT(:,end)==2&abs(tabParamScoreWT(:,2))<outlierLimit,2);tabParamScoreKO(tabParamScoreKO(:,end)==2&abs(tabParamScoreKO(:,2))<outlierLimit,2)],[ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==2&abs(tabParamScoreWT(:,2))<outlierLimit,1)));ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==2&abs(tabParamScoreKO(:,2))<outlierLimit,1)))*2],'off')
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(tabParamScoreWT(:,end)==2,2);tabParamScoreKO(tabParamScoreKO(:,end)==2,2)],[ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==2,1)));ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==2,1)))*2],'off')

% comparing alphas in random condition (1)
num2str([mean(tabParamScoreWT(tabParamScoreWT(:,end)==1,1)) mean(tabParamScoreKO(tabParamScoreKO(:,end)==1,1))])
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(tabParamScoreWT(:,end)==1,1);tabParamScoreKO(tabParamScoreKO(:,end)==1,1)],[ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==1,1)));ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==1,1)))*2],'off')
% comparing betas in random condition (1)
num2str([mean(tabParamScoreWT(tabParamScoreWT(:,end)==1&abs(tabParamScoreWT(:,2))<outlierLimit,2)) mean(tabParamScoreKO(tabParamScoreKO(:,end)==1&abs(tabParamScoreKO(:,2))<outlierLimit,2))])
[size(tabParamScoreWT(tabParamScoreWT(:,end)==1,2)) size(tabParamScoreWT(tabParamScoreWT(:,end)==1&abs(tabParamScoreWT(:,2))<outlierLimit,2)) size(tabParamScoreKO(tabParamScoreKO(:,end)==1,2)) size(tabParamScoreKO(tabParamScoreKO(:,end)==1&abs(tabParamScoreKO(:,2))<outlierLimit,2))]
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(tabParamScoreWT(:,end)==1,2);tabParamScoreKO(tabParamScoreKO(:,end)==1,2)],[ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==1,1)));ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==1,1)))*2],'off')
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(tabParamScoreWT(:,end)==1&abs(tabParamScoreWT(:,2))<outlierLimit,2);tabParamScoreKO(tabParamScoreKO(:,end)==1&abs(tabParamScoreKO(:,2))<outlierLimit,2)],[ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==1&abs(tabParamScoreWT(:,2))<outlierLimit,1)));ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==1&abs(tabParamScoreKO(:,2))<outlierLimit,1)))*2],'off')

% figure abs(beta) overlapping vs. stable vs. random
figure
plot(1-0.1+rand(size(abs(tabParamScoreWT(tabParamScoreWT(:,end)==3,2))))/5,abs(tabParamScoreWT(tabParamScoreWT(:,end)==3,2)),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
hold on
plot(2-0.1+rand(size(abs(tabParamScoreKO(tabParamScoreKO(:,end)==3,2))))/5,abs(tabParamScoreKO(tabParamScoreKO(:,end)==3,2)),'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
plot(3-0.1+rand(size(abs(tabParamScoreWT(tabParamScoreWT(:,end)==2,2))))/5,abs(tabParamScoreWT(tabParamScoreWT(:,end)==2,2)),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
plot(4-0.1+rand(size(abs(tabParamScoreKO(tabParamScoreKO(:,end)==2,2))))/5,abs(tabParamScoreKO(tabParamScoreKO(:,end)==2,2)),'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
plot(5-0.1+rand(size(abs(tabParamScoreWT(tabParamScoreWT(:,end)==1,2))))/5,abs(tabParamScoreWT(tabParamScoreWT(:,end)==1,2)),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
plot(6-0.1+rand(size(abs(tabParamScoreKO(tabParamScoreKO(:,end)==1,2))))/5,abs(tabParamScoreKO(tabParamScoreKO(:,end)==1,2)),'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
alpha 0.5
boxplot(abs([tabParamScoreWT(tabParamScoreWT(:,end)==3,2);tabParamScoreKO(tabParamScoreKO(:,end)==3,2);tabParamScoreWT(tabParamScoreWT(:,end)==2,2);tabParamScoreKO(tabParamScoreKO(:,end)==2,2);tabParamScoreWT(tabParamScoreWT(:,end)==1,2);tabParamScoreKO(tabParamScoreKO(:,end)==1,2)]),[ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==3,1)));ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==3,1)))*2;ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==2,1)))*3;ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==2,1)))*4;ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==1,1)))*5;ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==1,1)))*6],'Colors','k','Symbol','')
plot([0 7],[25.5 25.5],'k')
plot([0 7],[24.5 24.5],'k')
plot([2.5 2.5],[-5 40],'k')
plot([4.5 4.5],[-5 40],'k')
axis([0.5 6.5 -5 40])
syms beta
ylabel(texlabel(abs(beta)))
xlabel('overlapping                         stable                         random')
% outliers : 0 5 0 0 5 2
plot([2 2 2 2 2],[27 28 29 30 31],'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
plot([5 5 5 5 5],[27 28 29 30 31],'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
plot([6 6],[27 28],'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
plot([1 2],[35 35],'k')
plot([3 4],[35 35],'k')
plot([5 6],[35 35],'k')
text(1.4,37,'***')
text(3.4,37,'n.s.')
text(5.4,37,'n.s.')
title('grey circle: WT; black triangle: KO')

% figure alpha overlapping vs. stable vs. random
figure
plot(1-0.1+rand(size(tabParamScoreWT(tabParamScoreWT(:,end)==3,1)))/5,tabParamScoreWT(tabParamScoreWT(:,end)==3,1),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
hold on
plot(2-0.1+rand(size(tabParamScoreKO(tabParamScoreKO(:,end)==3,1)))/5,tabParamScoreKO(tabParamScoreKO(:,end)==3,1),'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
plot(3-0.1+rand(size(tabParamScoreWT(tabParamScoreWT(:,end)==2,1)))/5,tabParamScoreWT(tabParamScoreWT(:,end)==2,1),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
plot(4-0.1+rand(size(tabParamScoreKO(tabParamScoreKO(:,end)==2,1)))/5,tabParamScoreKO(tabParamScoreKO(:,end)==2,1),'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
plot(5-0.1+rand(size(tabParamScoreWT(tabParamScoreWT(:,end)==1,1)))/5,tabParamScoreWT(tabParamScoreWT(:,end)==1,1),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
plot(6-0.1+rand(size(tabParamScoreKO(tabParamScoreKO(:,end)==1,1)))/5,tabParamScoreKO(tabParamScoreKO(:,end)==1,1),'^k','MarkerFaceColor','k','MarkerEdgeColor','k')
alpha 0.5
boxplot([tabParamScoreWT(tabParamScoreWT(:,end)==3,1);tabParamScoreKO(tabParamScoreKO(:,end)==3,1);tabParamScoreWT(tabParamScoreWT(:,end)==2,1);tabParamScoreKO(tabParamScoreKO(:,end)==2,1);tabParamScoreWT(tabParamScoreWT(:,end)==1,1);tabParamScoreKO(tabParamScoreKO(:,end)==1,1)],[ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==3,1)));ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==3,1)))*2;ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==2,1)))*3;ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==2,1)))*4;ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==1,1)))*5;ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==1,1)))*6],'Colors','k','Symbol','')
plot([2.5 2.5],[-1 2],'k')
plot([4.5 4.5],[-1 2],'k')
axis([0.5 6.5 -0.05 1.3])
syms alpha
ylabel(texlabel(alpha))
xlabel('overlapping                         stable                         random')
% outliers : 0 5 0 0 5 2
plot([1 2],[1.1 1.1],'k')
plot([3 4],[1.1 1.1],'k')
plot([5 6],[1.1 1.1],'k')
text(1.4,1.2,'n.s.')
text(3.4,1.2,'n.s.')
text(5.4,1.2,'n.s.')
title('grey circle: WT; black triangle: KO')
clear alpha

% figure abs(beta) overlapping by bins of 5 trials
figure
data1 = vectData(1:sum(vectName==1),1);
name1 = vectName(1:sum(vectName==1),1);
plot(1-0.1+rand(sum(vectName==1),1)/5,abs(data1),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',2)
hold on
data2 = vectData(sum(vectName==1)+1:sum(vectName==1)+sum(vectName==2),1);
name2 = vectName(sum(vectName==1)+1:sum(vectName==1)+sum(vectName==2),1);
plot(2-0.1+rand(sum(vectName==2),1)/5,abs(data2),'^k','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',2)
data3 = vectData(sum(vectName==1)+sum(vectName==2)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3),1);
name3 = vectName(sum(vectName==1)+sum(vectName==2)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3),1);
plot(3-0.1+rand(sum(vectName==3),1)/5,abs(data3),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',2)
data4 = vectData(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4),1);
name4 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4),1);
plot(4-0.1+rand(sum(vectName==4),1)/5,abs(data4),'^k','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',2)
data5 = vectData(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5),1);
name5 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5),1);
plot(5-0.1+rand(sum(vectName==5),1)/5,abs(data5),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',2)
data6 = vectData(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6),1);
name6 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6),1);
plot(6-0.1+rand(sum(vectName==6),1)/5,abs(data6),'^k','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',2)
data7 = vectData(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7),1);
name7 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7),1);
plot(7-0.1+rand(sum(vectName==7),1)/5,abs(data7),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',2)
data8 = vectData(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7)+sum(vectName==8),1);
name8 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7)+sum(vectName==8),1);
plot(8-0.1+rand(sum(vectName==8),1)/5,abs(data8),'^k','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',2)
alpha 0.1
boxplot(abs(vectData),vectName,'Colors','k','Symbol','')
plot([2.5 2.5],[-20 900],'k')
plot([4.5 4.5],[-20 900],'k')
plot([6.5 6.5],[-20 900],'k')
axis([0.5 8.5 -20 900])
plot([0.5 8.5],[805 805],'k')
plot([0.5 8.5],[795 795],'k')
syms beta
ylabel(texlabel(abs(beta)))
xlabel('Bins of 5 trials')
% outliers : 3 3 1 2 1 3 0 1
plot([1 1 1],[830 850 870],'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',2)
plot([2 2 2],[830 850 870],'^k','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',2)
plot([3],[830],'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',2)
plot([4 4],[830 850],'^k','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',2)
plot([5],[830],'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',2)
plot([6 6 6],[830 850 870],'^k','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',2)
plot([8],[830],'^k','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',2)
plot([1 2],[520 520],'k')
plot([3 4],[520 520],'k')
plot([5 6],[520 520],'k')
plot([7 8],[520 520],'k')
text(1.4,550,'n.s.')
text(3.4,550,'n.s.')
text(5.4,550,'**')
text(7.4,550,'n.s.')
title('grey circle: WT; black triangle: KO')
[p,tbl,stats] = kruskalwallis(abs([data1;data2]),[name1;name2],'off')
[p,tbl,stats] = kruskalwallis(abs([data3;data4]),[name3;name4],'off')
[p,tbl,stats] = kruskalwallis(abs([data5;data6]),[name5;name6],'off')
[p,tbl,stats] = kruskalwallis(abs([data7;data8]),[name7;name8],'off')

% figure alpha overlapping by bins of 5 trials
figure
data1 = vectAlpha(1:sum(vectName==1),1);
name1 = vectName(1:sum(vectName==1),1);
plot(1-0.1+rand(sum(vectName==1),1)/5,abs(data1),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',2)
hold on
data2 = vectAlpha(sum(vectName==1)+1:sum(vectName==1)+sum(vectName==2),1);
name2 = vectName(sum(vectName==1)+1:sum(vectName==1)+sum(vectName==2),1);
plot(2-0.1+rand(sum(vectName==2),1)/5,abs(data2),'^k','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',2)
data3 = vectAlpha(sum(vectName==1)+sum(vectName==2)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3),1);
name3 = vectName(sum(vectName==1)+sum(vectName==2)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3),1);
plot(3-0.1+rand(sum(vectName==3),1)/5,abs(data3),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',2)
data4 = vectAlpha(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4),1);
name4 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4),1);
plot(4-0.1+rand(sum(vectName==4),1)/5,abs(data4),'^k','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',2)
data5 = vectAlpha(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5),1);
name5 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5),1);
plot(5-0.1+rand(sum(vectName==5),1)/5,abs(data5),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',2)
data6 = vectAlpha(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6),1);
name6 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6),1);
plot(6-0.1+rand(sum(vectName==6),1)/5,abs(data6),'^k','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',2)
data7 = vectAlpha(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7),1);
name7 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7),1);
plot(7-0.1+rand(sum(vectName==7),1)/5,abs(data7),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',2)
data8 = vectAlpha(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7)+sum(vectName==8),1);
name8 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7)+sum(vectName==8),1);
plot(8-0.1+rand(sum(vectName==8),1)/5,abs(data8),'^k','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',2)
alpha 0.1
boxplot(abs(vectAlpha),vectName,'Colors','k','Symbol','')
plot([2.5 2.5],[-20 900],'k')
plot([4.5 4.5],[-20 900],'k')
plot([6.5 6.5],[-20 900],'k')
axis([0.5 8.5 -0.05 1.3])
plot([0.5 8.5],[805 805],'k')
plot([0.5 8.5],[795 795],'k')
syms alpha
ylabel(texlabel(alpha))
xlabel('Bins of 5 trials')
% outliers : 3 3 1 2 1 3 0 1
plot([1 2],[1.1 1.1],'k')
plot([3 4],[1.1 1.1],'k')
plot([5 6],[1.1 1.1],'k')
plot([7 8],[1.1 1.1],'k')
text(1.4,1.2,'n.s.')
text(3.4,1.2,'n.s.')
text(5.4,1.2,'n.s.')
text(7.4,1.2,'n.s.')
title('grey circle: WT; black triangle: KO')
[p,tbl,stats] = kruskalwallis(abs([data1;data2]),[name1;name2],'off')
[p,tbl,stats] = kruskalwallis(abs([data3;data4]),[name3;name4],'off')
[p,tbl,stats] = kruskalwallis(abs([data5;data6]),[name5;name6],'off')
[p,tbl,stats] = kruskalwallis(abs([data7;data8]),[name7;name8],'off')
clear alpha

%% FIGURE BETA+ BETA- WITHOUT ELIMINATING OUTLIERS
data1 = vectData(1:sum(vectName==1),1);
name1 = vectName(1:sum(vectName==1),1);
data2 = vectData(sum(vectName==1)+1:sum(vectName==1)+sum(vectName==2),1);
name2 = vectName(sum(vectName==1)+1:sum(vectName==1)+sum(vectName==2),1);
data3 = vectData(sum(vectName==1)+sum(vectName==2)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3),1);
name3 = vectName(sum(vectName==1)+sum(vectName==2)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3),1);
data4 = vectData(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4),1);
name4 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4),1);
data5 = vectData(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5),1);
name5 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5),1);
data6 = vectData(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6),1);
name6 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6),1);
data7 = vectData(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7),1);
name7 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7),1);
data8 = vectData(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7)+sum(vectName==8),1);
name8 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7)+sum(vectName==8),1);
figure
plot([1 4],[0 0],'--')
hold on
errorfill(1:4, [mean(data1(data1>=0)) mean(data3(data3>=0)) mean(data5(data5>=0)) mean(data7(data7>=0))], [std(data1(data1>=0))/size(data1(data1>=0),1) std(data3(data3>=0))/size(data3(data3>=0),1) std(data5(data5>=0))/size(data5(data5>=0),1) std(data7(data7>=0))/size(data7(data7>=0),1)], '4')
errorfill(1:4, [mean(data1(data1<0)) mean(data3(data3<0)) mean(data3(data3<0)) mean(data3(data3<0))], [std(data1(data1<0))/size(data1(data1<0),1) std(data3(data3<0))/size(data3(data3<0),1) std(data5(data5<0))/size(data5(data5<0),1) std(data7(data7<0))/size(data7(data7<0),1)], '4')
errorfill(1:4, [mean(data2(data2>=0)) mean(data4(data4>=0)) mean(data6(data6>=0)) mean(data8(data8>=0))], [std(data2(data2>=0))/size(data2(data2>=0),1) std(data4(data4>=0))/size(data4(data4>=0),1) std(data6(data6>=0))/size(data6(data6>=0),1) std(data8(data8>=0))/size(data8(data8>=0),1)], 'k')
errorfill(1:4, [mean(data2(data2<0)) mean(data4(data4<0)) mean(data6(data6<0)) mean(data8(data8<0))], [std(data2(data2<0))/size(data2(data2<0),1) std(data4(data4<0))/size(data4(data4<0),1) std(data6(data6<0))/size(data6(data6<0),1) std(data8(data8<0))/size(data8(data8<0),1)], 'k')
syms beta
ylabel(texlabel(beta))
%yticks([-400 -300 -200 -100 0 100 200 300 400])
xticks(1:4)
xlabel('Bins of 5 trials')
title('grey circle: WT; black triangle: KO')

%% FIGURE WHEN ELIMINATING OUTLIERS abs(beta)>10000
% plots of beta separately for neophilic and neophobic
errorbar(4, mean(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2)), std(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2))/size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2),1), '+r','LineWidth',2)
hold on
errorbar(4, mean(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2)), std(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2))/size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2),1), '+r','LineWidth',2)
errorbar(4, mean(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2)), std(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2))/size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2),1), '+b','LineWidth',2)
errorbar(4, mean(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2)), std(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2))/size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2),1), '+b','LineWidth',2)
[mean(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2)) mean(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2)) mean(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2)) mean(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2))]
[std(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2))/size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2),1) std(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2))/size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2),1) std(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2))/size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2),1) std(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2))/size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2),1)]
[size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2),1) size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2),1) size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2),1) size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2),1)]
%errorbar(4,
%    3.6487   -2.3759    3.8058   -4.9023
%    0.3486    0.4087    0.5303    1.1640
%    18    14    15     7
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2);tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2)],[ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2)));ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2)))*2],'off')
%    {'Source'}    {'SS'        }    {'df'}    {'MS'      }    {'Chi-sq'  }    {'Prob>Chi-sq'}
%    {'Groups'}    {[   76.3889]}    {[ 1]}    {[ 76.3889]}    {[  0.8170]}    {[     0.3661]}
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2);tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2)],[ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2)));ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2)))*2],'off')
%    {'Source'}    {'SS'      }    {'df'}    {'MS'      }    {'Chi-sq'  }    {'Prob>Chi-sq'}
%    {'Groups'}    {[  9.0536]}    {[ 1]}    {[  9.0536]}    {[  0.2353]}    {[     0.6276]}
%errorbar(3, 
%    5.5603   -3.0936   84.8797 -234.5668
%    0.4240    0.5475   17.3750  130.9816
%    17    14    18     3
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2);tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2)],[ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2)));ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2)))*2],'off')
%    {'Source'}    {'SS'        }    {'df'}    {'MS'      }    {'Chi-sq'  }    {'Prob>Chi-sq'}
%    {'Groups'}    {[  132.2222]}    {[ 1]}    {[132.2222]}    {[  1.2593]}    {[     0.2618]}
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2);tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2)],[ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2)));ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2)))*2],'off')
%    {'Source'}    {'SS'      }    {'df'}    {'MS'      }    {'Chi-sq'  }    {'Prob>Chi-sq'}
%    {'Groups'}    {[ 14.5714]}    {[ 1]}    {[ 14.5714]}    {[  0.5714]}    {[     0.4497]}
%errorbar(2
%    23.6799  -65.6529  318.5156 -248.0037
%    4.8507   13.8769   73.1882   85.0626
%    20    11    18     5
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2);tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2)],[ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2)));ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2)))*2],'off')
%    {'Source'}    {'SS'        }    {'df'}    {'MS'      }    {'Chi-sq'  }    {'Prob>Chi-sq'}
%    {'Groups'}    {[  233.1722]}    {[ 1]}    {[233.1722]}    {[  1.8880]}    {[     0.1694]}
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2);tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2)],[ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2)));ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2)))*2],'off')
%    {'Source'}    {'SS'      }    {'df'}    {'MS'      }    {'Chi-sq'  }    {'Prob>Chi-sq'}
%    {'Groups'}    {[ 45.4545]}    {[ 1]}    {[ 45.4545]}    {[  2.0053]}    {[     0.1567]}
%errorbar(1,
%    51.4322 -210.8015  139.6058 -122.5171
%    11.9353   23.3381   22.3776   47.2119
%    15    14    15     6
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2);tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2)],[ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2)));ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2)))*2],'off')
%    {'Source'}    {'SS'        }    {'df'}    {'MS'      }    {'Chi-sq'  }    {'Prob>Chi-sq'}
%    {'Groups'}    {[   86.7000]}    {[ 1]}    {[ 86.7000]}    {[  1.1187]}    {[     0.2902]}
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2);tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2)],[ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2)));ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2)))*2],'off')
%    {'Source'}    {'SS'      }    {'df'}    {'MS'      }    {'Chi-sq'  }    {'Prob>Chi-sq'}
%    {'Groups'}    {[  0.9524]}    {[ 1]}    {[  0.9524]}    {[  0.0272]}    {[     0.8690]}
syms beta
ylabel(texlabel(beta))
%yticks([-400 -300 -200 -100 0 100 200 300 400])
xticks(1:4)
xlabel('Bins of 5 trials')
plot([1 4],[0 0],'--')
plot(1:4,[51.4322 23.6799 5.5603 3.6487], [0.5 0.5 0.5],'LineWidth',2)
plot(1:4,[-210.8015 -65.6529 -3.0936 -2.3759], [0.5 0.5 0.5],'LineWidth',2)
plot(1:4,[139.6058 318.5156 84.8797 3.8058], 'k','LineWidth',2)
plot(1:4,[-122.5171 -248.0037 -234.5668 -4.9023], 'k','LineWidth',2)
print('-bestfit','EHMT1_WT-KO_betaByBinsOf5trials','-dpdf')

%% FIGURE WHEN ELIMINATING OUTLIERS abs(beta)>10000
data1 = vectData(1:sum(vectName==1),1);
name1 = vectName(1:sum(vectName==1),1);
data2 = vectData(sum(vectName==1)+1:sum(vectName==1)+sum(vectName==2),1);
name2 = vectName(sum(vectName==1)+1:sum(vectName==1)+sum(vectName==2),1);
data3 = vectData(sum(vectName==1)+sum(vectName==2)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3),1);
name3 = vectName(sum(vectName==1)+sum(vectName==2)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3),1);
data4 = vectData(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4),1);
name4 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4),1);
data5 = vectData(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5),1);
name5 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5),1);
data6 = vectData(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6),1);
name6 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6),1);
data7 = vectData(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7),1);
name7 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7),1);
data8 = vectData(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7)+sum(vectName==8),1);
name8 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7)+sum(vectName==8),1);
figure
plot([1 4],[0 0],'--')
hold on
errorfill(1:4, [51.4322 23.6799 5.5603 3.6487], [11.9353 4.8507 0.4240 0.3486], '3')
errorfill(1:4, [-210.8015 -65.6529 -3.0936 -2.3759], [23.3381 13.8769 0.5475 0.4087], '3')
errorfill(1:4, [139.6058 318.5156 84.8797 3.8058], [22.3776 73.1882 17.3750 0.5303], 'k')
errorfill(1:4, [-122.5171 -248.0037 -234.5668 -4.9023], [47.2119 85.0626 130.9816 1.1640], 'k')
alpha 0.5
hold on
plot(1:4, [51.4322 23.6799 5.5603 3.6487], 'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
plot(1:4, [-210.8015 -65.6529 -3.0936 -2.3759], 'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
plot(1:4, [139.6058 318.5156 84.8797 3.8058], '^k','MarkerFaceColor','k','MarkerEdgeColor','k')
plot(1:4, [-122.5171 -248.0037 -234.5668 -4.9023], '^k','MarkerFaceColor','k','MarkerEdgeColor','k')
ylabel(texlabel(beta))
%yticks([-400 -300 -200 -100 0 100 200 300 400])
xticks(1:4)
xlabel('Bins of 5 trials')
title('grey circle: WT; black triangle: KO')

%% FIGURE WHEN ELIMINATING OUTLIERS abs(beta)>10000
data1 = vectData(1:sum(vectName==1),1);
name1 = vectName(1:sum(vectName==1),1);
data2 = vectData(sum(vectName==1)+1:sum(vectName==1)+sum(vectName==2),1);
name2 = vectName(sum(vectName==1)+1:sum(vectName==1)+sum(vectName==2),1);
data3 = vectData(sum(vectName==1)+sum(vectName==2)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3),1);
name3 = vectName(sum(vectName==1)+sum(vectName==2)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3),1);
data4 = vectData(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4),1);
name4 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4),1);
data5 = vectData(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5),1);
name5 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5),1);
data6 = vectData(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6),1);
name6 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6),1);
data7 = vectData(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7),1);
name7 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7),1);
data8 = vectData(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7)+sum(vectName==8),1);
name8 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7)+sum(vectName==8),1);
figure
plot([1 4],[0 0],'--')
hold on
errorfill(1:4, [mean(data1(data1>=0&abs(data1)<outlierLimit)) mean(data3(data3>=0&abs(data3)<outlierLimit)) mean(data5(data5>=0&abs(data5)<outlierLimit)) mean(data7(data7>=0&abs(data7)<outlierLimit))], [std(data1(data1>=0&abs(data1)<outlierLimit))/size(data1(data1>=0&abs(data1)<outlierLimit),1) std(data3(data3>=0&abs(data3)<outlierLimit))/size(data3(data3>=0&abs(data3)<outlierLimit),1) std(data5(data5>=0&abs(data5)<outlierLimit))/size(data5(data5>=0&abs(data5)<outlierLimit),1) std(data7(data7>=0&abs(data7)<outlierLimit))/size(data7(data7>=0&abs(data7)<outlierLimit),1)], '4')
errorfill(1:4, [mean(data1(data1<0&abs(data1)<outlierLimit)) mean(data3(data3<0&abs(data3)<outlierLimit)) mean(data5(data5<0&abs(data5)<outlierLimit)) mean(data7(data7<0&abs(data7)<outlierLimit))], [std(data1(data1<0&abs(data1)<outlierLimit))/size(data1(data1<0&abs(data1)<outlierLimit),1) std(data3(data3<0&abs(data3)<outlierLimit))/size(data3(data3<0&abs(data3)<outlierLimit),1) std(data5(data5<0&abs(data5)<outlierLimit))/size(data5(data5<0&abs(data5)<outlierLimit),1) std(data7(data7<0&abs(data7)<outlierLimit))/size(data7(data7<0&abs(data7)<outlierLimit),1)], '4')
errorfill(1:4, [mean(data2(data2>=0&abs(data2)<outlierLimit)) mean(data4(data4>=0&abs(data4)<outlierLimit)) mean(data6(data6>=0&abs(data6)<outlierLimit)) mean(data8(data8>=0&abs(data8)<outlierLimit))], [std(data2(data2>=0&abs(data2)<outlierLimit))/size(data2(data2>=0&abs(data2)<outlierLimit),1) std(data4(data4>=0&abs(data4)<outlierLimit))/size(data4(data4>=0&abs(data4)<outlierLimit),1) std(data6(data6>=0&abs(data6)<outlierLimit))/size(data6(data6>=0&abs(data6)<outlierLimit),1) std(data8(data8>=0&abs(data8)<outlierLimit))/size(data8(data8>=0&abs(data8)<outlierLimit),1)], 'k')
errorfill(1:4, [mean(data2(data2<0&abs(data2)<outlierLimit)) mean(data4(data4<0&abs(data4)<outlierLimit)) mean(data6(data6<0&abs(data6)<outlierLimit)) mean(data8(data8<0&abs(data8)<outlierLimit))], [std(data2(data2<0&abs(data2)<outlierLimit))/size(data2(data2<0&abs(data2)<outlierLimit),1) std(data4(data4<0&abs(data4)<outlierLimit))/size(data4(data4<0&abs(data4)<outlierLimit),1) std(data6(data6<0&abs(data6)<outlierLimit))/size(data6(data6<0&abs(data6)<outlierLimit),1) std(data8(data8<0&abs(data8)<outlierLimit))/size(data8(data8<0&abs(data8)<outlierLimit),1)], 'k')
alpha 0.5
hold on
plot(1:4, [mean(data1(data1>=0&abs(data1)<outlierLimit)) mean(data3(data3>=0&abs(data3)<outlierLimit)) mean(data5(data5>=0&abs(data5)<outlierLimit)) mean(data7(data7>=0&abs(data7)<outlierLimit))], 'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
plot(1:4, [mean(data1(data1<0&abs(data1)<outlierLimit)) mean(data3(data3<0&abs(data3)<outlierLimit)) mean(data5(data5<0&abs(data5)<outlierLimit)) mean(data7(data7<0&abs(data7)<outlierLimit))], 'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
plot(1:4, [mean(data2(data2>=0&abs(data2)<outlierLimit)) mean(data4(data4>=0&abs(data4)<outlierLimit)) mean(data6(data6>=0&abs(data6)<outlierLimit)) mean(data8(data8>=0&abs(data8)<outlierLimit))], '^k','MarkerFaceColor','k','MarkerEdgeColor','k')
plot(1:4, [mean(data2(data2<0&abs(data2)<outlierLimit)) mean(data4(data4<0&abs(data4)<outlierLimit)) mean(data6(data6<0&abs(data6)<outlierLimit)) mean(data8(data8<0&abs(data8)<outlierLimit))], '^k','MarkerFaceColor','k','MarkerEdgeColor','k')
syms beta
ylabel(texlabel(beta))
%yticks([-400 -300 -200 -100 0 100 200 300 400])
xticks(1:4)
xlabel('Bins of 5 trials')
title('grey circle: WT; black triangle: KO')

%% FIGURE ALPHA>0.5 ALPHA<0.5
data1 = vectAlpha(1:sum(vectName==1),1);
name1 = vectName(1:sum(vectName==1),1);
data2 = vectAlpha(sum(vectName==1)+1:sum(vectName==1)+sum(vectName==2),1);
name2 = vectName(sum(vectName==1)+1:sum(vectName==1)+sum(vectName==2),1);
data3 = vectAlpha(sum(vectName==1)+sum(vectName==2)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3),1);
name3 = vectName(sum(vectName==1)+sum(vectName==2)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3),1);
data4 = vectAlpha(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4),1);
name4 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4),1);
data5 = vectAlpha(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5),1);
name5 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5),1);
data6 = vectAlpha(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6),1);
name6 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6),1);
data7 = vectAlpha(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7),1);
name7 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7),1);
data8 = vectAlpha(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7)+sum(vectName==8),1);
name8 = vectName(sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7)+1:sum(vectName==1)+sum(vectName==2)+sum(vectName==3)+sum(vectName==4)+sum(vectName==5)+sum(vectName==6)+sum(vectName==7)+sum(vectName==8),1);
figure
plot([1 4],[0.5 0.5],'--')
hold on
errorfill(1:4, [mean(data1(data1>=0.5)) mean(data3(data3>=0.5)) mean(data5(data5>=0.5)) mean(data7(data7>=0.5))], [std(data1(data1>=0.5))/size(data1(data1>=0.5),1) std(data3(data3>=0.5))/size(data3(data3>=0.5),1) std(data5(data5>=0.5))/size(data5(data5>=0.5),1) std(data7(data7>=0.5))/size(data7(data7>=0.5),1)], '4')
errorfill(1:4, [mean(data1(data1<0.5)) mean(data3(data3<0.5)) mean(data3(data3<0.5)) mean(data3(data3<0.5))], [std(data1(data1<0.5))/size(data1(data1<0.5),1) std(data3(data3<0.5))/size(data3(data3<0.5),1) std(data5(data5<0.5))/size(data5(data5<0.5),1) std(data7(data7<0.5))/size(data7(data7<0.5),1)], '4')
errorfill(1:4, [mean(data2(data2>=0.5)) mean(data4(data4>=0.5)) mean(data6(data6>=0.5)) mean(data8(data8>=0.5))], [std(data2(data2>=0.5))/size(data2(data2>=0.5),1) std(data4(data4>=0.5))/size(data4(data4>=0.5),1) std(data6(data6>=0.5))/size(data6(data6>=0.5),1) std(data8(data8>=0.5))/size(data8(data8>=0.5),1)], 'k')
errorfill(1:4, [mean(data2(data2<0.5)) mean(data4(data4<0.5)) mean(data6(data6<0.5)) mean(data8(data8<0.5))], [std(data2(data2<0.5))/size(data2(data2<0.5),1) std(data4(data4<0.5))/size(data4(data4<0.5),1) std(data6(data6<0.5))/size(data6(data6<0.5),1) std(data8(data8<0.5))/size(data8(data8<0.5),1)], 'k')
alpha 0.5
hold on
plot(1:4, [mean(data1(data1>=0.5)) mean(data3(data3>=0.5)) mean(data5(data5>=0.5)) mean(data7(data7>=0.5))], 'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
plot(1:4, [mean(data1(data1<0.5)) mean(data3(data3<0.5)) mean(data3(data3<0.5)) mean(data3(data3<0.5))], 'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
plot(1:4, [mean(data2(data2>=0.5)) mean(data4(data4>=0.5)) mean(data6(data6>=0.5)) mean(data8(data8>=0.5))], '^k','MarkerFaceColor','k','MarkerEdgeColor','k')
plot(1:4, [mean(data2(data2<0.5)) mean(data4(data4<0.5)) mean(data6(data6<0.5)) mean(data8(data8<0.5))], '^k','MarkerFaceColor','k','MarkerEdgeColor','k')
syms alpha
ylabel(texlabel(alpha))
%yticks([-400 -300 -200 -100 0 100 200 300 400])
xticks(1:4)
xlabel('Bins of 5 trials')
title('grey circle: WT; black triangle: KO')
clear alpha

%% FIGURE WHEN ELIMINATING OUTLIERS abs(beta)>400
% plots of beta separately for neophilic and neophobic
errorbar(4, mean(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2)), std(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2))/size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2),1), '+r','LineWidth',2)
hold on
errorbar(4, mean(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2)), std(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2))/size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2),1), '+r','LineWidth',2)
errorbar(4, mean(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2)), std(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2))/size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2),1), '+b','LineWidth',2)
errorbar(4, mean(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2)), std(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2))/size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2),1), '+b','LineWidth',2)
[mean(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2)) mean(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2)) mean(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2)) mean(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2))]
[std(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2))/size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2),1) std(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2))/size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2),1) std(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2))/size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2),1) std(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2))/size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2),1)]
[size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2),1) size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2),1) size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2),1) size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2),1)]
%errorbar(4,
%    3.6487   -2.3759    3.8058   -4.9023
%    0.3486    0.4087    0.5303    1.1640
%    18    14    15     7
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2);tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2)],[ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2)));ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2)))*2],'off')
%    {'Source'}    {'SS'        }    {'df'}    {'MS'      }    {'Chi-sq'  }    {'Prob>Chi-sq'}
%    {'Groups'}    {[   76.3889]}    {[ 1]}    {[ 76.3889]}    {[  0.8170]}    {[     0.3661]}
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2);tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2)],[ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2)));ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2)))*2],'off')
%    {'Source'}    {'SS'      }    {'df'}    {'MS'      }    {'Chi-sq'  }    {'Prob>Chi-sq'}
%    {'Groups'}    {[  9.0536]}    {[ 1]}    {[  9.0536]}    {[  0.2353]}    {[     0.6276]}
%errorbar(3, 
%    5.5603   -3.0936   11.2564   -7.7438
%    0.4240    0.5475    0.9490    5.4588
%    17    14    17     2
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2);tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2)],[ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2)));ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2)))*2],'off')
%    {'Source'}    {'SS'        }    {'df'}    {'MS'      }    {'Chi-sq'  }    {'Prob>Chi-sq'}
%    {'Groups'}    {[  132.2222]}    {[ 1]}    {[132.2222]}    {[  1.2593]}    {[     0.2618]}
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2);tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2)],[ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2)));ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2)))*2],'off')
%    {'Source'}    {'SS'      }    {'df'}    {'MS'      }    {'Chi-sq'  }    {'Prob>Chi-sq'}
%    {'Groups'}    {[ 14.5714]}    {[ 1]}    {[ 14.5714]}    {[  0.5714]}    {[     0.4497]}
%errorbar(2
%    2.0045  -24.8944    8.0151  -64.3518
%    0.2119    7.4737    0.6603   31.9515
%    19    10    17     4
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2);tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2)],[ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2)));ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2)))*2],'off')
%    {'Source'}    {'SS'        }    {'df'}    {'MS'      }    {'Chi-sq'  }    {'Prob>Chi-sq'}
%    {'Groups'}    {[  233.1722]}    {[ 1]}    {[233.1722]}    {[  1.8880]}    {[     0.1694]}
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2);tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2)],[ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2)));ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2)))*2],'off')
%    {'Source'}    {'SS'      }    {'df'}    {'MS'      }    {'Chi-sq'  }    {'Prob>Chi-sq'}
%    {'Groups'}    {[ 45.4545]}    {[ 1]}    {[ 45.4545]}    {[  2.0053]}    {[     0.1567]}
%errorbar(1,
%    5.2483  -12.9651   15.3064   -6.8898
%    0.5600    2.2060    3.0339    1.1086
%    14    10    13     5
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2);tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2)],[ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)>=0,2)));ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)>=0,2)))*2],'off')
%    {'Source'}    {'SS'        }    {'df'}    {'MS'      }    {'Chi-sq'  }    {'Prob>Chi-sq'}
%    {'Groups'}    {[   86.7000]}    {[ 1]}    {[ 86.7000]}    {[  1.1187]}    {[     0.2902]}
[p,tbl,stats] = kruskalwallis([tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2);tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2)],[ones(size(tabParamScoreWT(tabParamScoreWT(:,end)==3&abs(tabParamScoreWT(:,2))<outlierLimit&tabParamScoreWT(:,2)<0,2)));ones(size(tabParamScoreKO(tabParamScoreKO(:,end)==3&abs(tabParamScoreKO(:,2))<outlierLimit&tabParamScoreKO(:,2)<0,2)))*2],'off')
%    {'Source'}    {'SS'      }    {'df'}    {'MS'      }    {'Chi-sq'  }    {'Prob>Chi-sq'}
%    {'Groups'}    {[  0.9524]}    {[ 1]}    {[  0.9524]}    {[  0.0272]}    {[     0.8690]}
syms beta
ylabel(texlabel(beta))
%yticks([-400 -300 -200 -100 0 100 200 300 400])
xticks(1:4)
xlabel('Bins of 5 trials')
plot([1 4],[0 0],'--')
plot(1:4,[51.4322 23.6799 5.5603 3.6487], [0.5 0.5 0.5],'LineWidth',2)
plot(1:4,[-210.8015 -65.6529 -3.0936 -2.3759], [0.5 0.5 0.5],'LineWidth',2)
plot(1:4,[139.6058 318.5156 84.8797 3.8058], 'k','LineWidth',2)
plot(1:4,[-122.5171 -248.0037 -234.5668 -4.9023], 'k','LineWidth',2)
print('-bestfit','EHMT1_WT-KO_betaByBinsOf5trials','-dpdf')

%% FIGURE WHEN ELIMINATING OUTLIERS abs(beta)>400
figure
plot([1 4],[0 0],'--')
hold on
errorfill(1:4, [5.2483 2.0045 5.5603 3.6487], [0.5600 0.2119 0.4240 0.3486], [0.5 0.5 0.5])
errorfill(1:4, [-12.9651 -24.8944 -3.0936 -2.3759], [2.2060 7.4737 0.5475 0.4087], [0.5 0.5 0.5])
errorfill(1:4, [15.3064 8.0151 11.2564 3.8058], [3.0339 0.6603 0.9490 0.5303], 'k')
errorfill(1:4, [-6.8898 -64.3518 -7.7438 -4.9023], [1.1086 31.9515 5.4588 1.1640], 'k')
ylabel(texlabel(beta))
%yticks([-400 -300 -200 -100 0 100 200 300 400])
xticks(1:4)
xlabel('Bins of 5 trials')
title('grey circle: WT; black triangle: KO')
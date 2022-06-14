%% Genzel et al. (2019) rat+mouse data

% figure of distrib of likelihoods
figure
subplot(2,4,1)
hist(tabParamScore(:,nbFreeParam+1))
xlabel('likelihood')
if (strcmp(species,'rat'))
    ylabel('number of rats')
else
    ylabel('number of mice')
end
%title(['Fit of Model ' num2str(whichModel) ' on mouse data'])
title('ALL conditions')
hold on
plot([chance_likelihood_restrictive(1) chance_likelihood_restrictive(1)],[0 22],'--k')
%axis([0.8 0.94 0 12])

% figure of distrib of likelihoods for RANDOM condition
subplot(2,4,2)
hist(tabParamScore(tabParamScore(:,end)==1,nbFreeParam+1))
xlabel('likelihood')
title('RANDOM condition')
hold on
plot([chance_likelihood_restrictive(2) chance_likelihood_restrictive(2)],[0 10],'--k')
%axis([0.74 0.96 0 10])

% figure of distrib of likelihoods for STABLE condition
subplot(2,4,3)
hist(tabParamScore(tabParamScore(:,end)==2,nbFreeParam+1))
xlabel('likelihood')
title('STABLE condition')
hold on
plot([chance_likelihood_restrictive(3) chance_likelihood_restrictive(3)],[0 10],'--k')
%axis([0.78 0.95 0 10])

% figure of distrib of likelihoods for OVERLAPPING condition
subplot(2,4,4)
hist(tabParamScore(tabParamScore(:,end)==3,nbFreeParam+1))
xlabel('likelihood')
title('OVERLAPPING condition')
hold on
plot([chance_likelihood_restrictive(4) chance_likelihood_restrictive(4)],[0 10],'--k')
%axis([0.84 0.95 0 10])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure of distrib of parameters
subplot(2,4,5)
theMax = 15; % 25 for rats 1.5 for mice
offset = theMax/3; space = theMax/20;
plot([0 1],[0 0],'--','Color',[0.5 0.5 0.5])
hold on
plot([0.5 0.5],[-theMax theMax],'--','Color',[0.5 0.5 0.5])
plot(tabParamScore(tabParamScore(:,1)>=0.5&tabParamScore(:,2)>=0,1),tabParamScore(tabParamScore(:,1)>=0.5&tabParamScore(:,2)>=0,2),'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
plot(tabParamScore(tabParamScore(:,1)<0.5&tabParamScore(:,2)>=0,1),tabParamScore(tabParamScore(:,1)<0.5&tabParamScore(:,2)>=0,2),'^k','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
plot(tabParamScore(tabParamScore(:,1)>=0.5&tabParamScore(:,2)<0,1),tabParamScore(tabParamScore(:,1)>=0.5&tabParamScore(:,2)<0,2),'ok','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
plot(tabParamScore(tabParamScore(:,1)<0.5&tabParamScore(:,2)<0,1),tabParamScore(tabParamScore(:,1)<0.5&tabParamScore(:,2)<0,2),'^k','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
syms alpha beta
xlabel(texlabel(alpha))
ylabel(texlabel(beta))
axis([0 1 -theMax theMax])
%title(['Fit of Model ' num2str(whichModel) ' on mouse data'])

% adding the outliers
plot([0 1],[theMax-offset+space/2 theMax-offset+space/2],'-k')
plot([0 1],[theMax-offset-space/2 theMax-offset-space/2],'-k')
plot([0 1],[-theMax+offset+space/2 -theMax+offset+space/2],'-k')
plot([0 1],[-theMax+offset-space/2 -theMax+offset-space/2],'-k')
outliers = tabParamScore(abs(tabParamScore(:,2))>theMax-offset+space/2,1:2);
outliersPos = outliers(outliers(:,2)>=0,:);
[boubou, ind] = sort(outliersPos(:,2),1,'descend');
outliersPos = outliersPos(ind,:);
nbOP = size(outliersPos,1);
for iii=1:nbOP
    if (outliersPos(iii,1)>=0.5)
        plot(outliersPos(iii,1),theMax-(offset-space/2)*iii/(nbOP+1),'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
    else
        plot(outliersPos(iii,1),theMax-(offset-space/2)*iii/(nbOP+1),'^k','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
    end
end
outliersNeg = outliers(outliers(:,2)<0,:);
[boubou, ind] = sort(outliersNeg(:,2),1,'descend');
outliersNeg = outliersNeg(ind,:);
nbON = size(outliersNeg,1);
for iii=1:nbON
    if (outliersNeg(iii,1)>=0.5)
        plot(outliersNeg(iii,1),-theMax+(offset+space/2)*iii/(nbON+1),'ok','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
    else
        plot(outliersNeg(iii,1),-theMax+(offset+space/2)*iii/(nbON+1),'^k','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
    end
end

% figure of distrib of parameters RANDOM session
subplot(2,4,6)
plot([0 1],[0 0],'--','Color',[0.5 0.5 0.5])
hold on
plot([0.5 0.5],[-theMax theMax],'--','Color',[0.5 0.5 0.5])
plot(tabParamScore(tabParamScore(:,1)>=0.5&tabParamScore(:,2)>=0&tabParamScore(:,end)==1,1),tabParamScore(tabParamScore(:,1)>=0.5&tabParamScore(:,2)>=0&tabParamScore(:,end)==1,2),'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
plot(tabParamScore(tabParamScore(:,1)<0.5&tabParamScore(:,2)>=0&tabParamScore(:,end)==1,1),tabParamScore(tabParamScore(:,1)<0.5&tabParamScore(:,2)>=0&tabParamScore(:,end)==1,2),'^k','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
plot(tabParamScore(tabParamScore(:,1)>=0.5&tabParamScore(:,2)<0&tabParamScore(:,end)==1,1),tabParamScore(tabParamScore(:,1)>=0.5&tabParamScore(:,2)<0&tabParamScore(:,end)==1,2),'ok','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
plot(tabParamScore(tabParamScore(:,1)<0.5&tabParamScore(:,2)<0&tabParamScore(:,end)==1,1),tabParamScore(tabParamScore(:,1)<0.5&tabParamScore(:,2)<0&tabParamScore(:,end)==1,2),'^k','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
syms alpha beta
xlabel(texlabel(alpha))
%ylabel(texlabel(beta))
%title(['Fit of Model ' num2str(whichModel) ' on mouse data RANDOM session'])
%axis([0 1 -0.8 1.7])
axis([0 1 -theMax theMax])

% adding the outliers
plot([0 1],[theMax-offset+space/2 theMax-offset+space/2],'-k')
plot([0 1],[theMax-offset-space/2 theMax-offset-space/2],'-k')
plot([0 1],[-theMax+offset+space/2 -theMax+offset+space/2],'-k')
plot([0 1],[-theMax+offset-space/2 -theMax+offset-space/2],'-k')
outliers = tabParamScore(tabParamScore(:,end)==1&abs(tabParamScore(:,2))>theMax-offset+space/2,1:2);
outliersPos = outliers(outliers(:,2)>=0,:);
[boubou, ind] = sort(outliersPos(:,2),1,'descend');
outliersPos = outliersPos(ind,:);
nbOP = size(outliersPos,1);
for iii=1:nbOP
    if (outliersPos(iii,1)>=0.5)
        plot(outliersPos(iii,1),theMax-(offset-space/2)*iii/(nbOP+1),'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
    else
        plot(outliersPos(iii,1),theMax-(offset-space/2)*iii/(nbOP+1),'^k','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
    end
end
outliersNeg = outliers(outliers(:,2)<0,:);
[boubou, ind] = sort(outliersNeg(:,2),1,'descend');
outliersNeg = outliersNeg(ind,:);
nbON = size(outliersNeg,1);
for iii=1:nbON
    if (outliersNeg(iii,1)>=0.5)
        plot(outliersNeg(iii,1),-theMax+(offset+space/2)*iii/(nbON+1),'ok','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
    else
        plot(outliersNeg(iii,1),-theMax+(offset+space/2)*iii/(nbON+1),'^k','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
    end
end

% figure of distrib of parameters STABLE session
subplot(2,4,7)
plot([0 1],[0 0],'--','Color',[0.5 0.5 0.5])
hold on
plot([0.5 0.5],[-theMax theMax],'--','Color',[0.5 0.5 0.5])
plot(tabParamScore(tabParamScore(:,1)>=0.5&tabParamScore(:,2)>=0&tabParamScore(:,end)==2,1),tabParamScore(tabParamScore(:,1)>=0.5&tabParamScore(:,2)>=0&tabParamScore(:,end)==2,2),'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
plot(tabParamScore(tabParamScore(:,1)<0.5&tabParamScore(:,2)>=0&tabParamScore(:,end)==2,1),tabParamScore(tabParamScore(:,1)<0.5&tabParamScore(:,2)>=0&tabParamScore(:,end)==2,2),'^k','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
plot(tabParamScore(tabParamScore(:,1)>=0.5&tabParamScore(:,2)<0&tabParamScore(:,end)==2,1),tabParamScore(tabParamScore(:,1)>=0.5&tabParamScore(:,2)<0&tabParamScore(:,end)==2,2),'ok','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
plot(tabParamScore(tabParamScore(:,1)<0.5&tabParamScore(:,2)<0&tabParamScore(:,end)==2,1),tabParamScore(tabParamScore(:,1)<0.5&tabParamScore(:,2)<0&tabParamScore(:,end)==2,2),'^k','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
syms alpha beta
xlabel(texlabel(alpha))
%ylabel(texlabel(beta))
%title(['Fit of Model ' num2str(whichModel) ' on mouse data STABLE session'])
axis([0 1 -theMax theMax])

% adding the outliers
plot([0 1],[theMax-offset+space/2 theMax-offset+space/2],'-k')
plot([0 1],[theMax-offset-space/2 theMax-offset-space/2],'-k')
plot([0 1],[-theMax+offset+space/2 -theMax+offset+space/2],'-k')
plot([0 1],[-theMax+offset-space/2 -theMax+offset-space/2],'-k')
outliers = tabParamScore(tabParamScore(:,end)==2&abs(tabParamScore(:,2))>theMax-offset+space/2,1:2);
outliersPos = outliers(outliers(:,2)>=0,:);
[boubou, ind] = sort(outliersPos(:,2),1,'descend');
outliersPos = outliersPos(ind,:);
nbOP = size(outliersPos,1);
for iii=1:nbOP
    if (outliersPos(iii,1)>=0.5)
        plot(outliersPos(iii,1),theMax-(offset-space/2)*iii/(nbOP+1),'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
    else
        plot(outliersPos(iii,1),theMax-(offset-space/2)*iii/(nbOP+1),'^k','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
    end
end
outliersNeg = outliers(outliers(:,2)<0,:);
[boubou, ind] = sort(outliersNeg(:,2),1,'descend');
outliersNeg = outliersNeg(ind,:);
nbON = size(outliersNeg,1);
for iii=1:nbON
    if (outliersNeg(iii,1)>=0.5)
        plot(outliersNeg(iii,1),-theMax+(offset+space/2)*iii/(nbON+1),'ok','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
    else
        plot(outliersNeg(iii,1),-theMax+(offset+space/2)*iii/(nbON+1),'^k','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
    end
end

% figure of distrib of parameters OVERLAPPING session
subplot(2,4,8)
plot([0 1],[0 0],'--','Color',[0.5 0.5 0.5])
hold on
plot([0.5 0.5],[-theMax theMax],'--','Color',[0.5 0.5 0.5])
plot(tabParamScore(tabParamScore(:,1)>=0.5&tabParamScore(:,2)>=0&tabParamScore(:,end)==3,1),tabParamScore(tabParamScore(:,1)>=0.5&tabParamScore(:,2)>=0&tabParamScore(:,end)==3,2),'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
plot(tabParamScore(tabParamScore(:,1)<0.5&tabParamScore(:,2)>=0&tabParamScore(:,end)==3,1),tabParamScore(tabParamScore(:,1)<0.5&tabParamScore(:,2)>=0&tabParamScore(:,end)==3,2),'^k','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
plot(tabParamScore(tabParamScore(:,1)>=0.5&tabParamScore(:,2)<0&tabParamScore(:,end)==3,1),tabParamScore(tabParamScore(:,1)>=0.5&tabParamScore(:,2)<0&tabParamScore(:,end)==3,2),'ok','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
plot(tabParamScore(tabParamScore(:,1)<0.5&tabParamScore(:,2)<0&tabParamScore(:,end)==3,1),tabParamScore(tabParamScore(:,1)<0.5&tabParamScore(:,2)<0&tabParamScore(:,end)==3,2),'^k','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
syms alpha beta
xlabel(texlabel(alpha))
%ylabel(texlabel(beta))
%title(['Fit of Model ' num2str(whichModel) ' on mouse data OVERLAPPING session'])
%axis([0 1 -0.8 1.7])
axis([0 1 -theMax theMax])

% adding the outliers
plot([0 1],[theMax-offset+space/2 theMax-offset+space/2],'-k')
plot([0 1],[theMax-offset-space/2 theMax-offset-space/2],'-k')
plot([0 1],[-theMax+offset+space/2 -theMax+offset+space/2],'-k')
plot([0 1],[-theMax+offset-space/2 -theMax+offset-space/2],'-k')
outliers = tabParamScore(tabParamScore(:,end)==3&abs(tabParamScore(:,2))>theMax-offset+space/2,1:2);
outliersPos = outliers(outliers(:,2)>=0,:);
[boubou, ind] = sort(outliersPos(:,2),1,'descend');
outliersPos = outliersPos(ind,:);
nbOP = size(outliersPos,1);
for iii=1:nbOP
    if (outliersPos(iii,1)>=0.5)
        plot(outliersPos(iii,1),theMax-(offset-space/2)*iii/(nbOP+1),'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
    else
        plot(outliersPos(iii,1),theMax-(offset-space/2)*iii/(nbOP+1),'^k','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
    end
end
outliersNeg = outliers(outliers(:,2)<0,:);
[boubou, ind] = sort(outliersNeg(:,2),1,'descend');
outliersNeg = outliersNeg(ind,:);
nbON = size(outliersNeg,1);
for iii=1:nbON
    if (outliersNeg(iii,1)>=0.5)
        plot(outliersNeg(iii,1),-theMax+(offset+space/2)*iii/(nbON+1),'ok','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
    else
        plot(outliersNeg(iii,1),-theMax+(offset+space/2)*iii/(nbON+1),'^k','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
    end
end

% print('-bestfit','ModelFig1A-H_simulationM1goodnessOfFitForRats','-dpdf')

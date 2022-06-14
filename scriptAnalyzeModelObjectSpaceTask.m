%% Genzel et al. (2019) rat+mouse data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RAT DATA

whichModel = 1;
nbFreeParam = 2;
subjectList = [1:16 101:116];
nbSubject = length(subjectList);
load(['rat_data_baseline_raw-copie.csv']);
DATA = rat_data_baseline_raw_copie;
DATA(DATA(:,1)>3,:) = []; % get rid of 'oa' conditions
DATA = DATA(:,[1 5:7 8:9 4 3]); % get rid of useless columns
nbTrial = 18;
population = [];
tabParamScore = []; % Param, L, LL, AIC, BIC, Pseudo-R2 scores of model, nbTrial
chance_level = 0.5;
chance_level_restrictive = 0;
for iii=1:nbSubject
    sDATA = DATA(DATA(:,end)==subjectList(iii),1:end-1);
    chance_LL_subject = log(simILonObjectSpaceTask(sDATA, 'optim', whichModel, 0, 0, 0, 0)) * nbTrial;
    chance_level_restrictive = chance_level_restrictive + chance_LL_subject;
    load(['rat/model' num2str(whichModel) '/OptimSpeciesrat_Model1_Init0_Subject' num2str(subjectList(iii)) '_fmsResults.mat'])
    [boubou, ind] = sort(fmsResults(:,4),1,'descend');
    dzack = fmsResults(fmsResults(:,1)>=0&fmsResults(:,1)<=1,:);
    [~, boubou] = max(dzack(:,end));
    population = [population ; dzack(boubou,[1 2 end])];
    % analyzing model fitting performance
    LL = log(dzack(boubou,end)) * nbTrial;
    tabParamScore = [tabParamScore ; [dzack(boubou,[1:end-2 end]) LL (-2*LL+nbFreeParam) (-2*LL+nbFreeParam*log(nbTrial)) (1-LL/chance_LL_subject) nbTrial]];
end
% outliers: 1 6 11 16 104 107 111 113 
LL = sum(tabParamScore(:,nbFreeParam+2)) % log-likelihood totale du modèle sur l'ensemble des sujets
NL = exp(LL/sum(tabParamScore(:,end))) % normalized likelihood, pour avoir une idée qualitative du fit du modèle en terme de proba (donc entre 0 et 1)
AIC = sum(tabParamScore(:,nbFreeParam+3)) % le modèle qui a le plus petit AIC fitte le mieux. Mais AIC ne pénalise pas assez la complexité (nbFreeParam) des modèles
BIC = -2 * sum(LL) + nbSubject * nbFreeParam * log(sum(tabParamScore(:,end))) % le modèle qui a le plus petit BIC fitte le mieux. Mais BIC pénalise trop la complexité (nbFreeParam) des modèles
% donc idéalement, si un modèle a le meilleur AIC et le meilleur BIC, alors il gagne !
pseudoR2 = 1 - LL / (sum(tabParamScore(:,end)) * log(chance_level)) % si > 0, alors le modèle fitte mieux que la chance
chance_LL_restrictive = chance_level_restrictive
chance_likelihood_restrictive = exp(sum(chance_level_restrictive)/sum(tabParamScore(:,end)))
pseudoR2_restrictive = 1 - LL / chance_level_restrictive % si > 0, alors le modèle fitte mieux que la chance

% figure of distrib of likelihoods
figure,hist(population(:,end))
xlabel('likelihood')
ylabel('number of rats')
title('Fit of Model 1 on rat data')
hold on
plot([chance_likelihood_restrictive chance_likelihood_restrictive],[0 12],'--k')
axis([0.8 0.94 0 12])
%print('-bestfit','distribLikelihoodRat_model1','-dpdf')

% figure of distrib of parameters
%figure,plot(population(:,1),population(:,2),'ok')
figure,plot([0 1],[0 0],'--','Color',[0.5 0.5 0.5])
hold on
plot([0.5 0.5],[-2 2],'--','Color',[0.5 0.5 0.5])
plot(population(population(:,1)>=0.5&population(:,2)>=0,1),population(population(:,1)>=0.5&population(:,2)>=0,2),'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
plot(population(population(:,1)<0.5&population(:,2)>=0,1),population(population(:,1)<0.5&population(:,2)>=0,2),'^k','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
plot(population(population(:,1)>=0.5&population(:,2)<0,1),population(population(:,1)>=0.5&population(:,2)<0,2),'ok','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
plot(population(population(:,1)<0.5&population(:,2)<0,1),population(population(:,1)<0.5&population(:,2)<0,2),'^k','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
axis([0 1 -2 2])
syms alpha beta
xlabel(texlabel(alpha))
ylabel(texlabel(beta))
title('Fit of Model 1 on rat data')

% adding the outliers
plot([0 1],[1.32 1.32],'-k')
plot([0 1],[1.28 1.28],'-k')
plot([0 1],[-1.52 -1.52],'-k')
plot([0 1],[-1.48 -1.48],'-k')

plot(0.99995855,1.9,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
text(0.9,1.9,'1360')
plot(0.99964167,1.8,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
text(0.9,1.8,'115')
plot(0.99633932,1.7,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
text(0.9,1.7,'35')
plot(0.99666582,1.6,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
text(0.9,1.6,'13.3')
plot(0.99359828,1.5,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
text(0.9,1.5,'13.1')
plot(0.99775778,-1.6,'ok','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
text(0.9,-1.6,'-11.5')
plot(0.99784084,-1.7,'ok','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
text(0.9,-1.7,'-15.7')
plot(0.99994323,-1.8,'ok','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
text(0.9,-1.8,'-504')

% % OLD (for only 9+4 starting points of gradient descent:)
% plot(0.99995317,1.9,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
% text(0.9,1.9,'1215')
% plot(0.9999739,1.8,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
% text(0.9,1.8,'680')
% plot(0.99338411,1.5,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
% text(0.9,1.5,'12.7')
% plot(0.99602257,1.4,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
% text(0.9,1.4,'11.4')
% plot(0.99994495,-1.8,'ok','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
% text(0.9,-1.8,'-520')
% plot(0.99756359,-1.7,'ok','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
% text(0.9,-1.7,'-14.1')

%print('-bestfit','distribParamRat_model1','-dpdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MOUSE DATA

whichModel = 3;
nbFreeParam = whichModel + 1;
subjectList = [1:8 101:108 201:208 26923:26927 26929:26930];
nbSubject = length(subjectList);
load(['mou_data_baseline_raw-copie.csv']);
DATA = mou_data_baseline_raw_copie;
DATA(DATA(:,1)>3,:) = []; % get rid of 'oa' conditions
DATA = DATA(:,[1 5:7 8:9 4 3]); % get rid of useless columns
%nbTrial = 62 or 63 depending on the subject (include in loop below)
tabParamScore = []; % Param, L, LL, AIC, BIC, Pseudo-R2 scores of model, nbTrial
chance_level = 0.5;
chance_level_restrictive = 0;
for iii=1:nbSubject
    sDATA = DATA(DATA(:,end)==subjectList(iii),1:end-1);
    nbTrial = size(sDATA,1);
    chance_LL_subject = log(simILonObjectSpaceTask(sDATA, 'optim', whichModel, 0, 0, 0, 0)) * nbTrial;
    chance_level_restrictive = chance_level_restrictive + chance_LL_subject;
    load(['mouse/model' num2str(whichModel) '/OptimSpeciesmou_Model' num2str(whichModel) '_Init025_Subject' num2str(subjectList(iii)) '_fmsResults.mat'])
    for iii=1:length(fmsResults)
        if (~isreal(fmsResults(iii,end)))
            fmsResults(iii,end) = 0;
        end
    end
    switch (whichModel)
        case 1
            dzack = fmsResults(fmsResults(:,1)>=0,:);
        case {2,3}
            dzack = fmsResults(fmsResults(:,1)>=0&fmsResults(:,3)>=0,:);
    end
    [~, boubou] = max(dzack(:,end));
    %population = [population ; dzack(boubou,[1 2 end])];
    % analyzing model fitting performance
    LL = log(dzack(boubou,end)) * nbTrial;
    tabParamScore = [tabParamScore ; [dzack(boubou,[1:end-2 end]) LL (-2*LL+nbFreeParam) (-2*LL+nbFreeParam*log(nbTrial)) (1-LL/chance_LL_subject) nbTrial]];
end
LL = sum(tabParamScore(:,nbFreeParam+2)) % log-likelihood totale du modèle sur l'ensemble des sujets
NL = exp(LL/sum(tabParamScore(:,end))) % normalized likelihood, pour avoir une idée qualitative du fit du modèle en terme de proba (donc entre 0 et 1)
AIC = sum(tabParamScore(:,nbFreeParam+3)) % le modèle qui a le plus petit AIC fitte le mieux. Mais AIC ne pénalise pas assez la complexité (nbFreeParam) des modèles
BIC = -2 * sum(LL) + nbSubject * nbFreeParam * log(sum(tabParamScore(:,end))) % le modèle qui a le plus petit BIC fitte le mieux. Mais BIC pénalise trop la complexité (nbFreeParam) des modèles
% donc idéalement, si un modèle a le meilleur AIC et le meilleur BIC, alors il gagne !
pseudoR2 = 1 - LL / (sum(tabParamScore(:,end)) * log(chance_level)) % si > 0, alors le modèle fitte mieux que la chance
chance_LL_restrictive = chance_level_restrictive
chance_likelihood_restrictive = exp(sum(chance_level_restrictive)/sum(tabParamScore(:,end)))
pseudoR2_restrictive = 1 - LL / chance_level_restrictive % si > 0, alors le modèle fitte mieux que la chance

% figure of distrib of likelihoods
figure,hist(tabParamScore(:,nbFreeParam+1))
xlabel('likelihood')
ylabel('number of mice')
title(['Fit of Model ' num2str(whichModel) ' on mouse data'])
hold on
plot([chance_likelihood_restrictive chance_likelihood_restrictive],[0 12],'--k')
axis([0.8 0.94 0 12])

% figure of distrib of parameters
%figure,plot(population(:,1),population(:,2),'ok')
figure,plot([0 1],[0 0],'--','Color',[0.5 0.5 0.5])
hold on
plot([0.5 0.5],[-0.8 1.7],'--','Color',[0.5 0.5 0.5])
plot(tabParamScore(tabParamScore(:,1)>=0.5&tabParamScore(:,2)>=0,1),tabParamScore(tabParamScore(:,1)>=0.5&tabParamScore(:,2)>=0,2),'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
plot(tabParamScore(tabParamScore(:,1)<0.5&tabParamScore(:,2)>=0,1),tabParamScore(tabParamScore(:,1)<0.5&tabParamScore(:,2)>=0,2),'^k','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
plot(tabParamScore(tabParamScore(:,1)>=0.5&tabParamScore(:,2)<0,1),tabParamScore(tabParamScore(:,1)>=0.5&tabParamScore(:,2)<0,2),'ok','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
plot(tabParamScore(tabParamScore(:,1)<0.5&tabParamScore(:,2)<0,1),tabParamScore(tabParamScore(:,1)<0.5&tabParamScore(:,2)<0,2),'^k','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
axis([0 1 -0.8 1.7])
syms alpha beta
xlabel(texlabel(alpha))
ylabel(texlabel(beta))
title(['Fit of Model ' num2str(whichModel) ' on mouse data'])

% adding the outliers
plot([0 1],[1.02 1.02],'-k')
plot([0 1],[0.98 0.98],'-k')
plot([0 1],[-0.52 -0.52],'-k')
plot([0 1],[-0.48 -0.48],'-k')

plot(0.9999892,1.6,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
text(0.9,1.6,'755')
plot(0.9999402,1.5,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
text(0.9,1.5,'171')
plot(0.999089,1.4,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
text(0.9,1.4,'20.8')
plot(0.9991393,1.3,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
text(0.9,1.3,'15.4')
plot(0.9992825,1.2,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
text(0.9,1.2,'9.4')
plot(0.9993826,1.1,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
text(0.9,1.1,'9.1')
plot(0.9963922,-0.6,'ok','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
text(0.9,-0.6,'-10.3')
plot(0.9988648,-0.7,'ok','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
text(0.9,-0.7,'-17.5')
% OLD (with only 9+4 starting points of the gradient descent:)
% plot(1.000000002,1.9,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
% text(0.9,1.9,'705272')
% plot(1.000000236,1.8,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
% text(0.9,1.8,'51097')
% plot(0.9999853963,1.7,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
% text(0.9,1.7,'570')
% plot(0.9999015594,1.6,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
% text(0.9,1.6,'109')
% plot(0.9988625845,1.5,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
% text(0.9,1.5,'17.1')
% plot(0.9989007987,1.4,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
% text(0.9,1.4,'12.4')
% plot(0.9951086749,1.3,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
% text(0.9,1.3,'11.8')
% plot(0.9991503345,1.2,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
% text(0.9,1.2,'8.1')
% plot(0.9992228474,1.1,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
% text(0.9,1.1,'7.4')
% plot(1.00000124,-1.9,'ok','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
% text(0.9,-1.9,'-15372')
% plot(0.9997232045,-1.2,'ok','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
% text(0.9,-1.2,'-8.0')
% plot(0.9962949536,-1.3,'ok','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
% text(0.9,-1.3,'-10.1')
% plot(0.9984450446,-1.4,'ok','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k')
% text(0.9,-1.4,'-13.3')

print('-bestfit','distribParamMouse_model1','-dpdf')

paramPerSubject = [subjectList' tabParamScore(:,[1:(nbFreeParam+1) end-1])];
save('mouseParamPerSubject_Model1init0.csv','paramPerSubject','-ASCII')

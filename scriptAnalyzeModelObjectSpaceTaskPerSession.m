%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MOUSE DATA

whichModel = 1;
nbSession = 3; % 3 considers the 3 conditions separately; 1 groups them together
species = 'R24'; %'mou'; % 'EWT'; % 'C24'

switch whichModel
    case 1
        init = '0';
        nbFreeParam = 2;
    case 2
        init = '0';
        nbFreeParam = 3;
    case 3
        init = '025';
        nbFreeParam = 4;
    case 4
        init = '025';
        nbFreeParam = 3;
end

%% Genzel et al. (2019) mouse data
%subjectList = [1:8 101:108 201:208 26923:26927 26929:26930]; % all data

%% Schut et al. (2020) mouse data
%subjectList = [1:35 201:213 301:307]; % all data
%subjectList = [1:7 11 13:14 16:17 18:19 20 23 26:27 30 32 34 201:203 205:208 212 302 304 307]; % WT
%subjectList = [8:10 12 15 21:22 24:25 28:29 31 33 35 204 209:211 213 301 303 305:306]; % KO
%subjectList = [1:4 11 18:19 201:203 205 20 23 26:27 30 204 209:211 302 304]; % E10
%subjectList = [5:7 13:14 16:17 32 34 206:208 212 307]; % E20
%subjectList = [8:10 12 21 22 24:25 28:29 31 301 303]; % E11
%subjectList = [15 33 35 213 305:306]; % E21

%% Lobato et al. (in prep) mouse data
%subjectList = [1:16 201:221 301:316 401:404 406:408 410:412 414:416]; %all data
switch species
    case 'CON'
        subjectList = [1:2:15 201:2:209 210 213:2:217 301:2:315 402:2:416]; % CON
    case 'C24'
        subjectList = [1:2:15 201:2:209 210 213:2:217 301:2:315 402:2:416]; % C24
    case 'C48'
        subjectList = [201:2:209 210 213:2:217 402:2:408 412:2:416]; % C48
    case 'C72'
        subjectList = [1:2:15 201:2:209 210 213:2:217]; % C72
    case 'RGS'
        subjectList = [2:2:16 202:2:208 211 212:2:218 219:221 302:2:316 401 403 407 411 415]; % RGS
    case 'R24'
        subjectList = [2:2:16 202:2:208 211 212:2:218 219:221 302:2:316 401 403 407 411 415]; % R24
    case 'R48'
        subjectList = [202:2:208 211 212:2:218 219:221 401 403 407 411]; % R48
    case 'R72'
        subjectList = [2:2:16 202:2:208 211 212:2:218 219:221]; % R72
end

nbSubject = length(subjectList);
switch (species)
    case 'rat'
        load(['rat_data_baseline_raw-copie.csv']);
        DATA = rat_data_baseline_raw_copie;
        clear rat_data_baseline_raw_copie;
    case 'mou'
        load(['mou_data_baseline_raw-copie.csv']);
        DATA = mou_data_baseline_raw_copie;
        clear mou_data_baseline_raw_copie;
    case 'EWT'
        load(['EHMT1_WT.mat']);
        DATA = EHMT1_WT;
        %DATA = DATA(:,[1:3 5 4 6:7 9:end]); % permute day/trial column and erase obj number
        clear EHMT1_WT;
    case 'EKO'
        load(['EHMT1_KO.mat']);
        DATA = EHMT1_KO;
        %DATA = DATA(:,[1:3 5 4 6:7 9:end]); % permute day/trial column and erase obj number
        clear EHMT1_KO;
    case 'E10'
        load(['EHMT1_group1gene0.mat']);
        DATA = EHMT1_group1gene0;
        DATA = DATA(:,[1:3 5 4 6:7 9:end]); % permute day/trial column and erase obj number
        clear EHMT1_group1gene0;
    case 'E11'
        load(['EHMT1_group1gene1.mat']);
        DATA = EHMT1_group1gene1;
        DATA = DATA(:,[1:3 5 4 6:7 9:end]); % permute day/trial column and erase obj number
        clear EHMT1_group1gene1;
    case 'E20'
        load(['EHMT1_group2gene0.mat']);
        DATA = EHMT1_group2gene0;
        DATA = DATA(:,[1:3 5 4 6:7 9:end]); % permute day/trial column and erase obj number
        clear EHMT1_group2gene0;
    case 'E21'
        load(['EHMT1_group2gene1.mat']);
        DATA = EHMT1_group2gene1;
        DATA = DATA(:,[1:3 5 4 6:7 9:end]); % permute day/trial column and erase obj number
        clear EHMT1_group2gene1;
    case 'CON' % Control, all tests (new data 2021)
        load(['Controlcorrected_24.mat']);
        DATA = Controlcorrected_24;
        load(['Controlcorrected_72.mat']);
        DATA = [DATA ; Controlcorrected_72];
        load(['Controlcorrected_48.mat']);
        Controlcorrected_48(:,1) = 1;
        DATA = [DATA ; Controlcorrected_48];
    case 'RGS' % RGS14, all tests (new data 2021)
        load(['RGS14corrected_24.mat']);
        DATA = RGS14corrected_24;
        load(['RGS14corrected_72.mat']);
        DATA = [DATA ; RGS14corrected_72];
        load(['RGS14corrected_48.mat']);
        RGS14corrected_48(:,1) = 1;
        DATA = [DATA ; RGS14corrected_48];
    case 'R24' % RGS14, 24h test (new data 2021)
        load(['RGS14corrected_24.mat']);
        DATA = RGS14corrected_24;
    case 'R72' % RGS14, 72h test (new data 2021)
        load(['RGS14corrected_72.mat']);
        DATA = RGS14corrected_72;
    case 'R48' % RGS14, 48h test (new data 2021)
        load(['RGS14corrected_48.mat']);
        DATA = RGS14corrected_48;
        DATA(:,1) = 1;
    case 'C24' % Control, 24h test (new data 2021)
        load(['Controlcorrected_24.mat']);
        DATA = Controlcorrected_24;
    case 'C72' % Control, 72h test (new data 2021)
        load(['Controlcorrected_72.mat']);
        DATA = Controlcorrected_72;
    case 'C48' % Control, 48h test (new data 2021)
        load(['Controlcorrected_48.mat']);
        DATA = Controlcorrected_48;
        DATA(:,1) = 1;
end
DATA(DATA(:,1)>3,:) = []; % get rid of 'oa' conditions
DATA = DATA(:,[1 5:7 8:9 4 3]); % get rid of useless columns    
% DATA should contain the series of trials (1 per line) and the following
    % columns:
    %     1 current condition (1 stable 2 random 3 overlapping)
    %     2 trial number
    %     3:4 locations
    %     5:6 proba
    %     7 day
    %     8 subject
%nbTrial = 62 or 63 depending on the subject (included in loop below)
tabParamScore = []; % Param, L, LL, AIC, BIC, Pseudo-R2 scores of model, nbTrial, subj, session (1,2,3)
chance_level = 0.5;
chance_level_restrictive = zeros(1,4);
for nsub=1:nbSubject
    for nses=1:nbSession
        if (nbSession > 1)
            sDATA = DATA(DATA(:,end)==subjectList(nsub)&DATA(:,1)==nses,1:end-1);
            %load(['WT-KO_removedTrialsWithDouble0s/OptimSpecies' species '_Model' num2str(whichModel) '_Init' init '_Subject' num2str(subjectList(nsub)) '_Session' num2str(nses) '_Trials16-20_fmsResults.mat']) % _Trials16-20
            %load(['mousePerSessionPerBinsOf5trials2019/OptimSpeciesmou_Model' num2str(whichModel) '_Init' init '_Subject' num2str(subjectList(nsub)) '_Session' num2str(nses) '_Trials16-20_fmsResults.mat'])
            load(['ModelOptimData/' species '/OptimSpecies' species '_Model' num2str(whichModel) '_Init' init '_Subject' num2str(subjectList(nsub)) '_Session' num2str(nses) '_fmsResults.mat']) % 2021: Lobato et al
        else
            load(['ModelOptimData/' species '/OptimSpecies' species '_Model' num2str(whichModel) '_Init' init '_Subject' num2str(subjectList(nsub)) '_fmsResults.mat']) % 2021: Lobato et al
        end
        nbTrial = size(sDATA,1);
        if (nbTrial > 0)%load(['mousePerSession/model' num2str(whichModel) '/OptimSpeciesmou_Model' num2str(whichModel) '_Init' init '_Subject' num2str(subjectList(nsub)) '_Session' num2str(nses) '_Trials1-5_fmsResults.mat'])
            chance_LL_subject = log(simILonObjectSpaceTask(sDATA, 'optim', whichModel, 0, 0, 0, 0)) * nbTrial;
            chance_level_restrictive = chance_level_restrictive + chance_LL_subject * [1 (nses==2) (nses==1) (nses==3)];
            switch (whichModel)
                case 1
                    dzack = fmsResults(fmsResults(:,1)>=0&fmsResults(:,1)<=1,:);
                case 2
                    dzack = fmsResults(fmsResults(:,1)>=0&fmsResults(:,1)<=1&fmsResults(:,3)>=0&fmsResults(:,3)<=1&fmsResults(:,end)>=0&fmsResults(:,end)<=1,:);
                case 3
                    for iii=1:length(fmsResults)
                        if (~isreal(fmsResults(iii,end)))
                            fmsResults(iii,end) = 0;
                        end
                    end
                    dzack = fmsResults(fmsResults(:,1)>=0&fmsResults(:,3)>=0&fmsResults(:,3)<=1&fmsResults(:,end)>=0&fmsResults(:,end)<=1,:); % &fmsResults(:,4)>=0
                case 4
                    for iii=1:length(fmsResults)
                        if (~isreal(fmsResults(iii,end)))
                            fmsResults(iii,end) = 0;
                        end
                    end
                    dzack = fmsResults(fmsResults(:,1)>=0&fmsResults(:,1)<=1&fmsResults(:,3)>=0&fmsResults(:,3)<=1&fmsResults(:,end)>=0&fmsResults(:,end)<=1,:); % &fmsResults(:,4)>=0
            end
            [~, boubou] = max(dzack(:,end));
            %population = [population ; dzack(boubou,[1 2 end])];
            % analyzing model fitting performance
            LL = log(dzack(boubou,end)) * nbTrial;
            tabParamScore = [tabParamScore ; [dzack(boubou,[1:end-2 end]) LL (-2*LL+nbFreeParam) (-2*LL+nbFreeParam*log(nbTrial)) (1-LL/chance_LL_subject) nbTrial subjectList(nsub) nses]];
        end
    end
end
LL = [sum(tabParamScore(:,nbFreeParam+2)) sum(tabParamScore(tabParamScore(:,end)==2,nbFreeParam+2)) sum(tabParamScore(tabParamScore(:,end)==1,nbFreeParam+2)) sum(tabParamScore(tabParamScore(:,end)==3,nbFreeParam+2))] % log-likelihood totale du modèle sur l'ensemble des sujets
NL = exp(LL ./ [sum(tabParamScore(:,end-2)) sum(tabParamScore(tabParamScore(:,end)==2,end-2)) sum(tabParamScore(tabParamScore(:,end)==1,end-2)) sum(tabParamScore(tabParamScore(:,end)==3,end-2))]) % normalized likelihood, pour avoir une idée qualitative du fit du modèle en terme de proba (donc entre 0 et 1)
AIC = [sum(tabParamScore(:,nbFreeParam+3)) sum(tabParamScore(tabParamScore(:,end)==2,nbFreeParam+3)) sum(tabParamScore(tabParamScore(:,end)==1,nbFreeParam+3)) sum(tabParamScore(tabParamScore(:,end)==3,nbFreeParam+3))] % le modèle qui a le plus petit AIC fitte le mieux. Mais AIC ne pénalise pas assez la complexité (nbFreeParam) des modèles
BIC = [sum(tabParamScore(:,nbFreeParam+4)) sum(tabParamScore(tabParamScore(:,end)==2,nbFreeParam+4)) sum(tabParamScore(tabParamScore(:,end)==1,nbFreeParam+4)) sum(tabParamScore(tabParamScore(:,end)==3,nbFreeParam+4))] % le modèle qui a le plus petit BIC fitte le mieux. Mais BIC pénalise trop la complexité (nbFreeParam) des modèles
pseudoR2 = (1-LL./[sum(tabParamScore(:,end-2))*log(chance_level) sum(tabParamScore(tabParamScore(:,end)==2,end-2))*log(chance_level) sum(tabParamScore(tabParamScore(:,end)==1,end-2))*log(chance_level) sum(tabParamScore(tabParamScore(:,end)==3,end-2))*log(chance_level)]) % si > 0, alors le modèle fitte mieux que la chance
chance_LL_restrictive = chance_level_restrictive
chance_likelihood_restrictive = exp(chance_level_restrictive ./ [sum(tabParamScore(:,end-2)) sum(tabParamScore(tabParamScore(:,end)==2,end-2)) sum(tabParamScore(tabParamScore(:,end)==1,end-2)) sum(tabParamScore(tabParamScore(:,end)==3,end-2))])
pseudoR2_restrictive = 1 - LL./chance_level_restrictive % si > 0, alors le modèle fitte mieux que la chance

% retransform wrong session numbers
tabParamScore(tabParamScore(:,end)==2,end)=4;
tabParamScore(tabParamScore(:,end)==1,end)=2;
tabParamScore(tabParamScore(:,end)==4,end)=1;

% plot mean+sem of parameters without outliers (abs(beta)<1000)
num2str(mean(tabParamScore(abs(tabParamScore(:,2))<1000,:)))
num2str(std(tabParamScore(abs(tabParamScore(:,2))<1000,:))/sqrt(sum(abs(tabParamScore(:,2))<1000)))

% % save data into csv file
% save('mouseParamPerSubjectPerSession_Model1init0.csv','tabParamScore','-ASCII')

%% Genzel et al. (2019)
% plot Figures 6 and 7 model simulation goodness of fit and distrib params
%scriptPlotAnalysisModelObjectSpaceTask.m

%% Schut et al 2020 (WT/KO)
%scriptObjectSpaceTaskPlotWTKO

%% Lobato et al. (in prep) (RGS14/CONTROL)
scriptObjectSpaceTaskPlotRGS14

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% RAT DATA (Genzel et al., 2019)
% 
% whichModel = 1;
% species = 'rat';
% init = '0';
% nbFreeParam = 2;
% subjectList = [1:16 101:116]; % rats Genzel et al 2019
% nbSubject = length(subjectList);
% nbSession = 3;
% load([species '_data_baseline_raw-copie.csv']);
% DATA = rat_data_baseline_raw_copie;
% DATA(DATA(:,1)>3,:) = []; % get rid of 'oa' conditions
% DATA = DATA(:,[1 5:7 8:9 4 3]); % get rid of useless columns
% nbTrial = 18;
% population = [];
% tabParamScore = []; % Param, L, LL, AIC, BIC, Pseudo-R2 scores of model, nbTrial, subj, session (1,2,3)
% chance_level = 0.5;
% chance_level_restrictive = 0;
% for nsub=1:nbSubject
%     for nses=1:nbSession
%         sDATA = DATA(DATA(:,end)==subjectList(nsub)&DATA(:,1)==nses,1:end-1);
%         chance_LL_subject = log(simILonObjectSpaceTask(sDATA, 'optim', whichModel, 0, 0, 0, 0)) * nbTrial;
%         chance_level_restrictive = chance_level_restrictive + chance_LL_subject * [1 (nses==2) (nses==1) (nses==3)];
%         load(['ratPerSession/model' num2str(whichModel) '/OptimSpeciesrat_Model' num2str(whichModel) '_Init' init '_Subject' num2str(subjectList(nsub)) '_Session' num2str(nses) '_fmsResults.mat'])
%         switch (whichModel)
%             case 1
%                 dzack = fmsResults(fmsResults(:,1)>=0&fmsResults(:,1)<=1,:);
%             case 2
%                 dzack = fmsResults(fmsResults(:,1)>=0&fmsResults(:,1)<=1&fmsResults(:,3)>=0&fmsResults(:,3)<=1&fmsResults(:,end)>=0&fmsResults(:,end)<=1,:);
%             case 3
%                 flag = false;
%                 for iii=1:length(fmsResults)
%                     if (~isreal(fmsResults(iii,end)))
%                         flag = true;
%                         fmsResults(iii,end) = 0;
%                     end
%                 end
% %                 if (flag)
% %                     [nsub nses]
% %                 end
%                 dzack = fmsResults(fmsResults(:,1)>=0&fmsResults(:,1)<=1&fmsResults(:,3)>=0&fmsResults(:,3)<=1&fmsResults(:,end)>=0&fmsResults(:,end)<=1,:); % &fmsResults(:,4)>=0
%             case 4
%                 flag = false;
%                 for iii=1:length(fmsResults)
%                     if (~isreal(fmsResults(iii,end)))
%                         flag = true;
%                         fmsResults(iii,end) = 0;
%                     end
%                 end
% %                 if (flag)
% %                     [nsub nses]
% %                 end
%                 dzack = fmsResults(fmsResults(:,1)>=0&fmsResults(:,1)<=1&fmsResults(:,3)>=0&fmsResults(:,3)<=1&fmsResults(:,end)>=0&fmsResults(:,end)<=1,:);
%         end
%         %[boubou, ind] = sort(fmsResults(:,end),1,'descend');
%         [~, boubou] = max(dzack(:,end));
%         % analyzing model fitting performance
%         LL = log(dzack(boubou,end)) * nbTrial;
%         tabParamScore = [tabParamScore ; [dzack(boubou,[1:end-2 end]) LL (-2*LL+nbFreeParam) (-2*LL+nbFreeParam*log(nbTrial)) (1-LL/chance_LL_subject) nbTrial subjectList(nsub) nses]];
%     end
% end
% % outliers: 1 6 11 16 104 107 111 113 
% LL = [sum(tabParamScore(:,nbFreeParam+2)) sum(tabParamScore(tabParamScore(:,end)==2,nbFreeParam+2)) sum(tabParamScore(tabParamScore(:,end)==1,nbFreeParam+2)) sum(tabParamScore(tabParamScore(:,end)==3,nbFreeParam+2))] % log-likelihood totale du modèle sur l'ensemble des sujets
% NL = exp(LL ./ [sum(tabParamScore(:,end-2)) sum(tabParamScore(tabParamScore(:,end)==2,end-2)) sum(tabParamScore(tabParamScore(:,end)==1,end-2)) sum(tabParamScore(tabParamScore(:,end)==3,end-2))]) % normalized likelihood, pour avoir une idée qualitative du fit du modèle en terme de proba (donc entre 0 et 1)
% AIC = [sum(tabParamScore(:,nbFreeParam+3)) sum(tabParamScore(tabParamScore(:,end)==2,nbFreeParam+3)) sum(tabParamScore(tabParamScore(:,end)==1,nbFreeParam+3)) sum(tabParamScore(tabParamScore(:,end)==3,nbFreeParam+3))] % le modèle qui a le plus petit AIC fitte le mieux. Mais AIC ne pénalise pas assez la complexité (nbFreeParam) des modèles
% BIC = [sum(tabParamScore(:,nbFreeParam+4)) sum(tabParamScore(tabParamScore(:,end)==2,nbFreeParam+4)) sum(tabParamScore(tabParamScore(:,end)==1,nbFreeParam+4)) sum(tabParamScore(tabParamScore(:,end)==3,nbFreeParam+4))] % le modèle qui a le plus petit BIC fitte le mieux. Mais BIC pénalise trop la complexité (nbFreeParam) des modèles
% pseudoR2 = (1-LL./[sum(tabParamScore(:,end-2))*log(chance_level) sum(tabParamScore(tabParamScore(:,end)==2,end-2))*log(chance_level) sum(tabParamScore(tabParamScore(:,end)==1,end-2))*log(chance_level) sum(tabParamScore(tabParamScore(:,end)==3,end-2))*log(chance_level)]) % si > 0, alors le modèle fitte mieux que la chance
% chance_LL_restrictive = chance_level_restrictive
% chance_likelihood_restrictive = exp(chance_level_restrictive ./ [sum(tabParamScore(:,end-2)) sum(tabParamScore(tabParamScore(:,end)==2,end-2)) sum(tabParamScore(tabParamScore(:,end)==1,end-2)) sum(tabParamScore(tabParamScore(:,end)==3,end-2))])
% pseudoR2_restrictive = 1 - LL./chance_level_restrictive % si > 0, alors le modèle fitte mieux que la chance
% 
% % retransform wrong session numbers
% tabParamScore(tabParamScore(:,end)==2,end)=4;
% tabParamScore(tabParamScore(:,end)==1,end)=2;
% tabParamScore(tabParamScore(:,end)==4,end)=1;
% 
% % plot mean+sem of parameters without outliers (abs(beta)<1000)
% num2str(mean(tabParamScore(abs(tabParamScore(:,2))<1000,:)))
% num2str(std(tabParamScore(abs(tabParamScore(:,2))<1000,:))/sqrt(sum(abs(tabParamScore(:,2))<1000)))
% 
% % % save data into csv file
% % save('ratParamPerSubjectPerSession_Model1init0.csv','tabParamScore','-ASCII')
% 
% %% Genzel et al. (2019)
% % plot Figures 6 and 7 model simulation goodness of fit and distrib params
% %scriptPlotAnalysisModelObjectSpaceTask.m
% 
% 

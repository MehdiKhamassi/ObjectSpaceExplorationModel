function [results, DATA] = launchSimILonObjectSpaceTask(species, whichModel, nbTrials, alpha, beta, gamma, init)
    % INPUT
    %       species = rat / mou
    %       nbTrials = nb trials per condition x 6
    %       alpha = learning rate (if < 0.5, overlapping>stable, else
    %                               overlapping<stable)
    %       beta = exploration rate (determines amplitude of D.I.)
    %               0.2 seems good for Fig 2C-D
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% INIT FIRST VERSION
%     % rat experiment taken from Fig 2A
%     DATA = [1 1 1 2 0 0;1 2 2 4 0 0; 1 3 3 4 0 0; 1 4 1 3 0 0; 1 5 4 2 0 0; 1 6 1 4 0 0; 2 1 1 2 0 0; 2 2 1 2 0 0; 2 3 1 2 0 0; 2 4 1 2 0 0; 2 5 1 2 0 0; 2 6 1 4 0 0; 3 1 1 2 0 0; 3 2 1 3 0 0; 3 3 1 2 0 0; 3 4 1 3 0 0; 3 5 1 4 0 0; 3 6 1 4 0 0];
%     DATA2= [1 1 3 4 0 0;1 2 4 2 0 0; 1 3 1 3 0 0; 1 4 4 1 0 0; 1 5 2 3 0 0; 1 6 1 4 0 0; 2 1 1 2 0 0; 2 2 1 2 0 0; 2 3 1 2 0 0; 2 4 1 2 0 0; 2 5 1 2 0 0; 2 6 1 4 0 0; 3 1 1 2 0 0; 3 2 1 3 0 0; 3 3 1 2 0 0; 3 4 1 3 0 0; 3 5 1 4 0 0; 3 6 1 4 0 0];
%     DATA3= [1 1 1 3 0 0;1 2 1 4 0 0; 1 3 3 2 0 0; 1 4 4 1 0 0; 1 5 2 4 0 0; 1 6 1 4 0 0; 2 1 1 2 0 0; 2 2 1 2 0 0; 2 3 1 2 0 0; 2 4 1 2 0 0; 2 5 1 2 0 0; 2 6 1 4 0 0; 3 1 1 2 0 0; 3 2 1 3 0 0; 3 3 1 2 0 0; 3 4 1 3 0 0; 3 5 1 4 0 0; 3 6 1 4 0 0];
%     DATA = [DATA ; DATA2];
%     DATA = [DATA ; DATA3];
%     DATA = [DATA ; DATA];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% INIT SECOND VERSION
%     % automatically generate object locations in DATA based on Fig.2A and 3A
%     if (strcmp(species,'rat'))
%         nbSample = 5; % (rat version)
%     else
%         nbSample = 20; % (mouse version)
%     end
%     DATA = zeros(nbTrials,7);
% %     Mehdi : DATA columns
% %     1 current condition (1 stable 2 random 3 overlapping)
% %     2 trial number
% %     3:4 locations
% %     5:6 proba
% %     7 day
%     testO1 = randi(4); testO2 = randi(3); if (testO2 >= testO1) testO2 = testO2 + 1; end
%     vecteur = 1:4; vecteur(vecteur==testO1) = []; vecteur(vecteur==testO2) = [];
%     for iii=1:size(DATA,1)
%         switch(mod(iii,3)+1)
%             case 1 % random
%                 for jjj=1:nbSample
%                     DATA((iii-1)*(nbSample+1)+jjj,:) = [mod(iii,3)+1 jjj 0 0 0 0 1];
%                     O1 = randi(4); O2 = randi(3); if (O2 >= O1) O2 = O2 + 1; end
%                     DATA((iii-1)*(nbSample+1)+jjj,3:4) = [O1 O2];
%                     % test trial
%                     DATA((iii-1)*(nbSample+1)+nbSample+1,:) = [mod(iii,3)+1 99 testO1 testO2 0 0 2];
%                 end
%             case 2 % stable
%                 testO1 = randi(4); testO2 = randi(3); if (testO2 >= testO1) testO2 = testO2 + 1; end
%                 vecteur = 1:4; vecteur(vecteur==testO1) = []; vecteur(vecteur==testO2) = [];
%                 O2 = randi(2); O2 = vecteur(O2);
%                 DATA((iii-1)*(nbSample+1)+1,:) = [mod(iii,3)+1 1 testO1 O2 0 0 1];
%                 for jjj=2:nbSample
%                     DATA((iii-1)*(nbSample+1)+jjj,:) = [mod(iii,3)+1 jjj testO1 O2 0 0 1];
%                 end
%                 % test trial
%                 DATA((iii-1)*(nbSample+1)+nbSample+1,:) = [mod(iii,3)+1 99 testO1 testO2 0 0 2];
%             case 3 % overlapping
%                 for jjj=1:(nbSample-1)
%                     O2 = randi(2); O2 = vecteur(O2);
%                     DATA((iii-1)*(nbSample+1)+jjj,:) = [mod(iii,3)+1 jjj testO1 O2 0 0 1];
%                 end
%                 % last sample trial
%                 DATA((iii-1)*(nbSample+1)+nbSample,:) = [mod(iii,3)+1 nbSample testO1 testO2 0 0 1];
%                 % test trial
%                 DATA((iii-1)*(nbSample+1)+nbSample+1,:) = [mod(iii,3)+1 99 testO1 testO2 0 0 2];
%         end
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% INIT THIRD VERSION
    % simulate the model on exactly the series of trials experienced by the
    % animals
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% loading behavioral data for all subjects
    switch (species)
        case 'rat'
            load(['rat_data_baseline_raw-copie.csv']);
            DATA = rat_data_baseline_raw_copie;
            clear rat_data_baseline_raw_copie;
        case 'mou'
            load(['mou_data_baseline_raw-copie.csv']);
            DATA = mou_data_baseline_raw_copie;
            clear mou_data_baseline_raw_copie;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% reorganizing columns of DATA for the model simulation
    % DATA should contain the series of trials (1 per line) and the following
    % columns:
    %     1 current condition (1 stable 2 random 3 overlapping)
    %     2 trial number
    %     3:4 locations
    %     5:6 proba
    %     7 day
    %     8 subject
    debut = 1; fin = 30000; % we want all subjects
    DATA(DATA(:,1)>3,:) = []; % get rid of 'oa' conditions
    DATA(DATA(:,3)<debut,:) = []; % get rid of non-considered subjects
    DATA(DATA(:,3)>fin,:) = []; % get rid of non-considered subjects
    DATA = DATA(:,[1 5:7 8:9 4 3]); % get rid of useless columns
    if strcmp(species,'rat')
        DATA(DATA(:,2)==6,7) = 2; % test trial is 24h later
    end
    %nbSubject = size(unique(DATA(:,end)),1)
    nbSubject = max(DATA(:,end));
    excludedSubjects = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% reorganizing lines of DATA into separate blocks per subject
    sDATA = [];
    nbSample = 1000;
    for nsub = debut:min(nbSubject,fin)
        if ((sum(excludedSubjects==nsub) == 0)&&(sum(DATA(:,end)==nsub)>0))
            % preventing 2 consecutive sessions to have the same number
            % once the subject number is removed, and the program can thus
            % think that it is the same continuing session without
            % resetting
            if ((strcmp(species,'rat'))&(nsub>100))
                sDATA = [sDATA ; DATA(DATA(:,end)==nsub&DATA(:,1)==2,:)];
                sDATA = [sDATA ; DATA(DATA(:,end)==nsub&DATA(:,1)==1,:)];
                sDATA = [sDATA ; DATA(DATA(:,end)==nsub&DATA(:,1)==3,:)];
            else if ((strcmp(species,'mou'))&(nsub==6))
                sDATA = [sDATA ; DATA(DATA(:,end)==nsub&DATA(:,1)==3,:)];
                sDATA = [sDATA ; DATA(DATA(:,end)==nsub&DATA(:,1)==1,:)];
                sDATA = [sDATA ; DATA(DATA(:,end)==nsub&DATA(:,1)==2,:)];
            else if ((strcmp(species,'mou'))&(nsub==101))
                sDATA = [sDATA ; DATA(DATA(:,end)==nsub&DATA(:,1)==1,:)];
                sDATA = [sDATA ; DATA(DATA(:,end)==nsub&DATA(:,1)==3,:)];
                sDATA = [sDATA ; DATA(DATA(:,end)==nsub&DATA(:,1)==2,:)];
            else if ((strcmp(species,'mou'))&(nsub==106))
                sDATA = [sDATA ; DATA(DATA(:,end)==nsub&DATA(:,1)==3,:)];
                sDATA = [sDATA ; DATA(DATA(:,end)==nsub&DATA(:,1)==1,:)];
                sDATA = [sDATA ; DATA(DATA(:,end)==nsub&DATA(:,1)==2,:)];
            else if ((strcmp(species,'mou'))&(nsub==203))
                sDATA = [sDATA ; DATA(DATA(:,end)==nsub&DATA(:,1)==2,:)];
                sDATA = [sDATA ; DATA(DATA(:,end)==nsub&DATA(:,1)==3,:)];
                sDATA = [sDATA ; DATA(DATA(:,end)==nsub&DATA(:,1)==1,:)];
            else if ((strcmp(species,'mou'))&(nsub==206))
                sDATA = [sDATA ; DATA(DATA(:,end)==nsub&DATA(:,1)==1,:)];
                sDATA = [sDATA ; DATA(DATA(:,end)==nsub&DATA(:,1)==3,:)];
                sDATA = [sDATA ; DATA(DATA(:,end)==nsub&DATA(:,1)==2,:)];
            else if ((strcmp(species,'mou'))&(nsub==208))
                sDATA = [sDATA ; DATA(DATA(:,end)==nsub&DATA(:,1)==3,:)];
                sDATA = [sDATA ; DATA(DATA(:,end)==nsub&DATA(:,1)==2,:)];
                sDATA = [sDATA ; DATA(DATA(:,end)==nsub&DATA(:,1)==1,:)];
            else if ((strcmp(species,'mou'))&(nsub==26929))
                sDATA = [sDATA ; DATA(DATA(:,end)==nsub&DATA(:,1)==3,:)];
                sDATA = [sDATA ; DATA(DATA(:,end)==nsub&DATA(:,1)==1,:)];
                sDATA = [sDATA ; DATA(DATA(:,end)==nsub&DATA(:,1)==2,:)];
            else
                sDATA = [sDATA ; DATA(DATA(:,end)==nsub,:)];
            end; end; end; end; end; end; end; end
            nbSample = min(nbSample,max(sDATA(:,2))-1);
        end
    end
    DATA = sDATA;
    clear sDATA;
    %num2str(DATA(DATA(:,2)==1,:)) % log
    %nbSample
    DATA = DATA(:,1:end-1); % removing the subject column
    DATA(DATA(:,2)==nbSample+1,2) = 99; % labeling test trial as 99
    DATA(DATA(:,1)==1,1) = 4; % wrong labeling of stable condition
    DATA(DATA(:,1)==2,1) = 1; % wrong labeling of random condition
    DATA(DATA(:,1)==4,1) = 2; % wrong labeling of stable condition
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TASK simulation
    [~, DATA] = simILonObjectSpaceTask(DATA, 'simul', whichModel, alpha, beta, gamma, init);
    % note that the gamma day-to-day forgetting rate is not used yet
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% FIGURES
    % figure indexes
    if (strcmp(species,'mou'))
        idx1 = 3:4;
        idx2 = 7:8;
        idx3 = 11:12;
    else
        idx1 = 1:2;
        idx2 = 5:6;
        idx3 = 9:10;
    end
    %figure
    %subplot(2,1,1)
    % plot Discrimination Index
    subplot(3,4,idx1)
    vecteurRANDOM = [];
    minSize = 1000; % we take the min nb of trials per condition (to equal nb of trials between conditions)
    for iii=1:nbSample
        minSize = min(minSize,size(DATA(DATA(:,1)==1&DATA(:,2)==iii,5)));
    end
    minSize = min(minSize,size(DATA(DATA(:,1)==1&DATA(:,2)==99,5)));
    for iii=1:nbSample
        %[iii size(DATA(DATA(:,1)==1&DATA(:,2)==iii,6)) size(DATA(DATA(:,1)==1&DATA(:,2)==iii,5))]
        DI = DATA(DATA(:,1)==1&DATA(:,2)==iii,6)-DATA(DATA(:,1)==1&DATA(:,2)==iii,5);
        vecteurRANDOM = [vecteurRANDOM DI(1:minSize)];
    end
    DI = DATA(DATA(:,1)==1&DATA(:,2)==99,6)-DATA(DATA(:,1)==1&DATA(:,2)==99,5);
    vecteurRANDOM = [vecteurRANDOM DI(1:minSize)];
    errorbar([1:nbSample+1]-0.1,mean(vecteurRANDOM,1),std(vecteurRANDOM,1),'.','Color', [0.5 0.5 0.5])
    hold on
    plot([1:nbSample+1]-0.1,mean(vecteurRANDOM,1), ':o', 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], 'Color', [0.5 0.5 0.5])
    
    vecteurSTABLE = [];
    minSize = 1000; % we take the min nb of trials per condition (to equal nb of trials between conditions)
    for iii=1:nbSample
        minSize = min(minSize,size(DATA(DATA(:,1)==2&DATA(:,2)==iii,5)));
    end
    minSize = min(minSize,size(DATA(DATA(:,1)==2&DATA(:,2)==99,5)));
    for iii=1:nbSample
        DI = DATA(DATA(:,1)==2&DATA(:,2)==iii,6)-DATA(DATA(:,1)==2&DATA(:,2)==iii,5);
        vecteurSTABLE = [vecteurSTABLE DI(1:minSize)];
    end
    DI = DATA(DATA(:,1)==2&DATA(:,2)==99,6)-DATA(DATA(:,1)==2&DATA(:,2)==99,5);
    vecteurSTABLE = [vecteurSTABLE DI(1:minSize)];
    errorbar([1:nbSample+1]-0,mean(vecteurSTABLE,1),std(vecteurRANDOM,1),'.k')
    plot([1:nbSample+1]-0,mean(vecteurSTABLE,1), '--s', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], 'Color', [0 0 0])
    
    vecteurOVERLA = [];
    minSize = 1000; % we take the min nb of trials per condition (to equal nb of trials between conditions)
    for iii=1:nbSample
        minSize = min(minSize,size(DATA(DATA(:,1)==3&DATA(:,2)==iii,5)));
    end
    minSize = min(minSize,size(DATA(DATA(:,1)==3&DATA(:,2)==99,5)));
    for iii=1:nbSample
        DI = DATA(DATA(:,1)==3&DATA(:,2)==iii,6)-DATA(DATA(:,1)==3&DATA(:,2)==iii,5);
        vecteurOVERLA = [vecteurOVERLA DI(1:minSize)];
    end
    DI = DATA(DATA(:,1)==3&DATA(:,2)==99,6)-DATA(DATA(:,1)==3&DATA(:,2)==99,5);
    vecteurOVERLA = [vecteurOVERLA DI(1:minSize)];
    errorbar([1:nbSample+1]+0.1,mean(vecteurOVERLA,1),std(vecteurRANDOM,1),'.k')
    plot([1:nbSample+1]+0.1,mean(vecteurOVERLA,1), '-^', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], 'Color', [0 0 0])
    
    theMax = 0;
    for iii = 1:nbSample+1
        if (abs(mean(vecteurRANDOM(:,iii)) + std(vecteurRANDOM(:,iii))) > theMax)
            theMax = abs(mean(vecteurRANDOM(:,iii)) + std(vecteurRANDOM(:,iii)));
        end
        if (abs(mean(vecteurSTABLE(:,iii)) + std(vecteurSTABLE(:,iii))) > theMax)
            theMax = abs(mean(vecteurSTABLE(:,iii)) + std(vecteurSTABLE(:,iii)));
        end
        if (abs(mean(vecteurOVERLA(:,iii)) + std(vecteurOVERLA(:,iii))) > theMax)
            theMax = abs(mean(vecteurOVERLA(:,iii)) + std(vecteurOVERLA(:,iii)));
        end
    end
    theMax = 1.5 * theMax;
    xlabel('trials')
    %xticks([1:(nbSample+1)])
    %xticklabels([num2str(1:nbSample) 'Test'])
    axis([0.5 (nbSample+1+0.5) min(-0.01,-theMax) max(0.01,theMax)]) %1.5*(max(mean([vecteurRANDOM;vecteurSTABLE;vecteurOVERLA],1))+max(std([vecteurRANDOM;vecteurSTABLE;vecteurOVERLA],1)))])
    if (idx1(1) == 1)
        ylabel('D.I.')
    end
    title(species)
    
    % plot Discrimination Index at Test
    %subplot(2,1,2)
    subplot(3,4,idx2)
    errorbar(1,mean(vecteurRANDOM(:,nbSample+1),1),std(vecteurRANDOM(:,nbSample+1),1),'.k')
    hold on
    errorbar(5,mean(vecteurSTABLE(:,nbSample+1),1),std(vecteurSTABLE(:,nbSample+1),1),'.k')
    errorbar(9,mean(vecteurOVERLA(:,nbSample+1),1),std(vecteurOVERLA(:,nbSample+1),1),'.k')
    bar(1,mean(vecteurRANDOM(:,nbSample+1),1), 'k')
    bar(5,mean(vecteurSTABLE(:,nbSample+1),1), 'k')
    bar(9,mean(vecteurOVERLA(:,nbSample+1),1), 'k')
    [h, p] = ttest(vecteurRANDOM(:,nbSample+1));
    if (p < 0.05)
        text(1,0.9*theMax,'*')
    end
    [h, p] = ttest(vecteurSTABLE(:,nbSample+1));
    if (p < 0.05)
        text(5,0.9*theMax,'*')
    end
    [h, p] = ttest(vecteurOVERLA(:,nbSample+1));
    if (p < 0.05)
        text(9,0.9*theMax,'*')
    end
    if (strcmp(species,'mou'))
        errorbar(0,mean(vecteurRANDOM(:,nbSample),1),std(vecteurRANDOM(:,nbSample),1),'.k')
        errorbar(4,mean(vecteurSTABLE(:,nbSample),1),std(vecteurSTABLE(:,nbSample),1),'.k')
        errorbar(8,mean(vecteurOVERLA(:,nbSample),1),std(vecteurOVERLA(:,nbSample),1),'.k')
        bar(0,mean(vecteurRANDOM(:,nbSample),1), 'w')
        bar(4,mean(vecteurSTABLE(:,nbSample),1), 'w')
        bar(8,mean(vecteurOVERLA(:,nbSample),1), 'w')
        [h, p] = ttest(vecteurRANDOM(:,nbSample));
        if (p < 0.05)
            text(0,0.9*theMax,'*')
        end
        [h, p] = ttest(vecteurSTABLE(:,nbSample));
        if (p < 0.05)
            text(4,0.9*theMax,'*')
        end
        [h, p] = ttest(vecteurOVERLA(:,nbSample));
        if (p < 0.05)
            text(8,0.9*theMax,'*')
        end
    end
    axis([-2 11 min(-0.01,-theMax) max(0.01,theMax)])
    xticks('')
    xlabel('Random   Stable   Overlapping')
    if (idx2(1) == 5)
        ylabel('D.I.')
    end
    
    %% evolution of entropy
    subplot(3,4,idx3)
    % random condition
    vecteurENTO1 = [];
    vecteurENTO2 = [];
    minSize = 1000; % we take the min nb of trials per condition (to equal nb of trials between conditions)
    for iii=1:nbSample
        minSize = min(minSize,size(DATA(DATA(:,1)==1&DATA(:,2)==iii,5)));
    end
    minSize = min(minSize,size(DATA(DATA(:,1)==1&DATA(:,2)==99,5)));
    for iii=1:nbSample
        ENTO1 = DATA(DATA(:,1)==1&DATA(:,2)==iii,8);
        ENTO2 = DATA(DATA(:,1)==1&DATA(:,2)==iii,9);
        vecteurENTO1 = [vecteurENTO1 ENTO1(1:minSize)];
        vecteurENTO2 = [vecteurENTO2 ENTO2(1:minSize)];
    end
    ENTO1 = DATA(DATA(:,1)==1&DATA(:,2)==99,8);
    ENTO2 = DATA(DATA(:,1)==1&DATA(:,2)==99,9);
    vecteurENTO1 = [vecteurENTO1 ENTO1(1:minSize)];
    vecteurENTO2 = [vecteurENTO2 ENTO2(1:minSize)];
    errorbar([1:nbSample+1]-0.1,mean(vecteurENTO2,1),std(vecteurENTO2,1),'.','Color', [0.5 0.5 0.5])
    hold on
    plot([1:nbSample+1]-0.1,mean(vecteurENTO2,1), ':o', 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], 'Color', [0.5 0.5 0.5])
    % stable condition
    vecteurENTO1 = [];
    vecteurENTO2 = [];
    minSize = 1000; % we take the min nb of trials per condition (to equal nb of trials between conditions)
    for iii=1:nbSample
        minSize = min(minSize,size(DATA(DATA(:,1)==2&DATA(:,2)==iii,5)));
    end
    minSize = min(minSize,size(DATA(DATA(:,1)==2&DATA(:,2)==99,5)));
    for iii=1:nbSample
        ENTO1 = DATA(DATA(:,1)==2&DATA(:,2)==iii,8);
        ENTO2 = DATA(DATA(:,1)==2&DATA(:,2)==iii,9);
        vecteurENTO1 = [vecteurENTO1 ENTO1(1:minSize)];
        vecteurENTO2 = [vecteurENTO2 ENTO2(1:minSize)];
    end
    ENTO1 = DATA(DATA(:,1)==2&DATA(:,2)==99,8);
    ENTO2 = DATA(DATA(:,1)==2&DATA(:,2)==99,9);
    vecteurENTO1 = [vecteurENTO1 ENTO1(1:minSize)];
    vecteurENTO2 = [vecteurENTO2 ENTO2(1:minSize)];
    errorbar([1:nbSample+1]-0,mean(vecteurENTO2,1),std(vecteurENTO2,1),'.k')
    plot([1:nbSample+1]-0,mean(vecteurENTO2,1), '--s', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], 'Color', [0 0 0])
    % overlapping condition
    vecteurENTO1 = [];
    vecteurENTO2 = [];
    minSize = 1000; % we take the min nb of trials per condition (to equal nb of trials between conditions)
    for iii=1:nbSample
        minSize = min(minSize,size(DATA(DATA(:,1)==3&DATA(:,2)==iii,5)));
    end
    minSize = min(minSize,size(DATA(DATA(:,1)==3&DATA(:,2)==99,5)));
    for iii=1:nbSample
        ENTO1 = DATA(DATA(:,1)==3&DATA(:,2)==iii,8);
        ENTO2 = DATA(DATA(:,1)==3&DATA(:,2)==iii,9);
        vecteurENTO1 = [vecteurENTO1 ENTO1(1:minSize)];
        vecteurENTO2 = [vecteurENTO2 ENTO2(1:minSize)];
    end
    ENTO1 = DATA(DATA(:,1)==3&DATA(:,2)==99,8);
    ENTO2 = DATA(DATA(:,1)==3&DATA(:,2)==99,9);
    vecteurENTO1 = [vecteurENTO1 ENTO1(1:minSize)];
    vecteurENTO2 = [vecteurENTO2 ENTO2(1:minSize)];
    errorbar([1:nbSample+1]+0.1,mean(vecteurENTO2,1),std(vecteurENTO2,1),'.k')
    plot([1:nbSample+1]+0.1,mean(vecteurENTO2,1), '-^', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], 'Color', [0 0 0])
    if (idx3(1) == 9)
        ylabel('entropy O2')
    end
    xlabel('trials')
    axis([0.5 (nbSample+1+0.5) -0.5 2]) %1.5*(max(mean([vecteurRANDOM;vecteurSTABLE;vecteurOVERLA],1))+max(std([vecteurRANDOM;vecteurSTABLE;vecteurOVERLA],1)))])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% EXTRACTION of main results
    results = zeros(1,8);
    results(1) = mean(vecteurSTABLE(:,nbSample+1),1); % stable D.I. at test
    results(2) = abs(mean(vecteurOVERLA(:,nbSample+1),1) - mean(vecteurSTABLE(:,nbSample+1),1)); % overlap-stable D.I. at test
    results(3) = abs(mean(vecteurOVERLA(:,nbSample),1) - mean(vecteurSTABLE(:,nbSample),1)); % overlap-stable D.I. at last sample
    results(4) = mean(vecteurSTABLE(:,nbSample),1); % stable D.I. at last sample
    results(5) = mean(vecteurOVERLA(:,nbSample+1),1); % overlapping D.I. at test
    results(6) = mean(vecteurOVERLA(:,nbSample),1); % overlapping D.I. at last sample
    results(7) = mean(vecteurRANDOM(:,nbSample+1),1); % random D.I. at test
    results(8) = mean(vecteurRANDOM(:,nbSample),1); % random D.I. at last sample
    
    % clearing memory when calling function with
    % scriptPlotSimILonObjectSpaceTask.m
    %clear DATA;
    clear sDATA;
    clear vecteurRANDOM;
    clear vecteurSTABLE;
    clear vecteurOVERLA;
end
function [likelihood, DATA] = simILonObjectSpaceTask(DATA, mode, whichModel, alpha, beta, gamma, init)

%     Mehdi : DATA columns
%     1 current condition (1 stable 2 random 3 overlapping)
%     2 trial number
%     3:4 locations
%     5:6 proba
%     7 day
%     8:9 entropies

%     mode = 'simul' / 'optim'

    % fixed parameters
    nC = 3; % 3 conditions
    nO = 2; % 2 objects
    nL = 4; % 4 locations
    
    %alpha = 0.6; % learning rate
    %beta = 0.3; % exploration rate
    %gamma = T.B.D. % day-to-day forgetting rate
    
    % input parameters
    nbTrials = size(DATA,1);
    likelihood = zeros(nbTrials,1); % loglikelihoods of data/reward given the model
    
    tabINIT = ones(nO,nL) * init; % USED to be zeros(); CONSIDER ones() or ones()*param here!!!
    tabOL = [];
    for ccc=1:nC
        tabOL = [tabOL ; tabINIT];
    end
    probaPerObject = zeros(nbTrials,nO);
    
    % launching the task
    for iii=1:nbTrials
        switch(DATA(iii,1)) % current condition
            case {1,6,7} % random condition
                % reset values for 2 new objects at first trial
                if (DATA(iii,2) == 1)
                    tabOL(1:2,:) = tabINIT;
                else
                    if (((whichModel == 2)||(whichModel == 3))&&(iii>1)&&(DATA(iii,7)~=DATA(iii-1,7)))
                        % day-to-day forgetting
                        tabOL(1,:) = tabOL(1,:) + (1 - gamma) * (tabINIT(1,:) - tabOL(1,:));
                        tabOL(2,:) = tabOL(2,:) + (1 - gamma) * (tabINIT(2,:) - tabOL(2,:));
                    end
                end
                % learning object locations
                O1location = zeros(1, nL); O1location(DATA(iii,3)) = 1;
                tabOL(1,:) = (1 - alpha) * tabOL(1,:) + alpha * O1location;
                O2location = zeros(1, nL); O2location(DATA(iii,4)) = 1;
                tabOL(2,:) = (1 - alpha) * tabOL(2,:) + alpha * O2location;
                % exploring objects
                entO1 = entropyProba(tabOL(1,:) / sum(tabOL(1,:)));
                entO2 = entropyProba(tabOL(2,:) / sum(tabOL(2,:)));
                [~, probaPerObject(iii,:)] = valueBasedDecision([entO1 entO2],'softmax',beta,0);
%                 % debug
%                 tabOL(1:2,:)
%                 [DATA(iii,1:2) entO1 entO2 probaPerObject(iii,:)]
            case 2 % stable condition
                % reset for 2 new objects
                if (DATA(iii,2) == 1)
                    tabOL(3:4,:) = tabINIT;
                else
                    if (((whichModel == 2)||(whichModel == 3))&&(iii>1)&&(DATA(iii,7)~=DATA(iii-1,7)))
                        % day-to-day forgetting
                        tabOL(3,:) = tabOL(3,:) + (1 - gamma) * (tabINIT(1,:) - tabOL(3,:));
                        tabOL(4,:) = tabOL(4,:) + (1 - gamma) * (tabINIT(2,:) - tabOL(4,:));
                    end
                end
                % learning object locations
                O1location = zeros(1, nL); O1location(DATA(iii,3)) = 1;
                tabOL(3,:) = (1 - alpha) * tabOL(3,:) + alpha * O1location;
                O2location = zeros(1, nL); O2location(DATA(iii,4)) = 1;
                tabOL(4,:) = (1 - alpha) * tabOL(4,:) + alpha * O2location;
                % exploring objects
                entO1 = entropyProba(tabOL(3,:) / sum(tabOL(3,:)));
                entO2 = entropyProba(tabOL(4,:) / sum(tabOL(4,:)));
                [~, probaPerObject(iii,:)] = valueBasedDecision([entO1 entO2],'softmax',beta,0);
%                 % debug
%                 tabOL(3:4,:)
%                 [DATA(iii,1:2) entO1 entO2 probaPerObject(iii,:)]
            case {3,4,5} % overlapping condition
                % reset for 2 new objects
                if (DATA(iii,2) == 1)
                    tabOL(5:6,:) = tabINIT;
                else
                    if (((whichModel == 2)||(whichModel == 3))&&(iii>1)&&(DATA(iii,7)~=DATA(iii-1,7)))
                        % day-to-day forgetting
                        tabOL(5,:) = tabOL(5,:) + (1 - gamma) * (tabINIT(1,:) - tabOL(5,:));
                        tabOL(6,:) = tabOL(6,:) + (1 - gamma) * (tabINIT(2,:) - tabOL(6,:));
                    end
                end
                % learning object locations
                O1location = zeros(1, nL); O1location(DATA(iii,3)) = 1;
                tabOL(5,:) = (1 - alpha) * tabOL(5,:) + alpha * O1location;
                O2location = zeros(1, nL); O2location(DATA(iii,4)) = 1;
                tabOL(6,:) = (1 - alpha) * tabOL(6,:) + alpha * O2location;
                % exploring objects
                entO1 = entropyProba(tabOL(5,:) / sum(tabOL(5,:)));
                entO2 = entropyProba(tabOL(6,:) / sum(tabOL(6,:)));
                [~, probaPerObject(iii,:)] = valueBasedDecision([entO1 entO2],'softmax',beta,0);
%                 % debug
%                 tabOL(5:6,:)
%                 [DATA(iii,1:2) entO1 entO2 probaPerObject(iii,:)]
        end % end of switch(condition)
        if (mode == 'simul')
            % storing probas and entropies
            DATA(iii,[5:6 8:9]) = [probaPerObject(iii,:) entO1 entO2];
        else % mode == 'optim'
            % comparing with subject's choice
            probaSubject = DATA(iii,5:6) / max(10e-6,sum(DATA(iii,5:6))); % avoid division by 0
            likelihood(iii) = log(1 - abs(probaSubject(1) - probaPerObject(iii,1))); % storing log likelihood
        end
    end % end of the task
    likelihood = exp(sum(likelihood)/nbTrials); % LL on all trials
%     %likelihood = exp(sum(likelihood(6:end))/(nbTrials-5)); % LL on trials 6:21
%     %likelihood = exp(sum(likelihood(11:end))/(nbTrials-10)); % LL on trials 11:21
%     lastTrial = 20;
%     if (nbTrials < lastTrial)
%         lastTrial = nbTrials;
%     end;
%     likelihood = exp(sum(likelihood(16:lastTrial))/(lastTrial-16+1)); % LL on 5 trials (1 week)
%     %likelihood = exp(sum(likelihood(11:lastTrial))/(lastTrial-16+1)); % LL on 5 trials (1 week)
end
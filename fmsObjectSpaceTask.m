function LL = fmsObjectSpaceTask( vectFreeParam, init, whichModel, DATA )
%fmsObjectSpaceTask performs an optimization of parameters minimizing the
%negative log likelihood of the model fitting data in the rat/mouse object
%space task of Lisa Genzel and Francesco Battaglia.
%It uses a gradient descent method (fminsearch or fmincon)
%initialized at a series of starting points

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fixed parameters
    nbTrials = size(DATA,1);
    idxLik = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % optimized parameters
    alpha = vectFreeParam(1);
    beta = vectFreeParam(2);
    switch (whichModel)
        case 4 % no day-to-day forgetting but adding an Init parameter
            gamma = 0;
            init = vectFreeParam(3);
        case 3 % hypothesizing a day-to-day forgetting and adding an Init parameter
            gamma = vectFreeParam(3);
            init = vectFreeParam(4);
        case 2 % hypothesizing a day-to-day forgetting
            gamma = vectFreeParam(3);
            %init = 0.25;
        otherwise % no day-to-day forgetting
            gamma = 0;
            %init = 0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % launching the simulation of the object space task with the required parameters
    likelihood = simILonObjectSpaceTask(DATA, 'optim', whichModel, alpha, beta, gamma, init);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % computing LL of the requested model summed over all trials of the session
    % this corresponds to an aggregate score over check_cue, object choice and RTs
    LL = - 1 * (log(likelihood(:,idxLik))*nbTrials); % on all trials
%     lastTrial = 20;
%     if (nbTrials < lastTrial)
%         lastTrial = nbTrials;
%     end;
%     LL = - 1 * (log(likelihood(:,idxLik))*(lastTrial-16+1)); % on 5 trials (one week)
end


%function parallelOptiObjectSpaceTask

    species = 'rat';
    methode = 'fms';
    whichModel = 4;
    init = 0.25;
    
    %% grid search + gradient descent
    %subjectList = [1:16 101:116]; % rat
    %subjectList = [1:8 101:108 201:208 26923:26927 26929:26930]; % mouse
    subjectList = [1:35 201:213 301:307]; % 2019
    nWorkers = length(subjectList);
    
    myCluster=parcluster('local'); 
    myCluster.NumWorkers= nWorkers; 
    poolobj = parpool(myCluster,nWorkers); 

    parfor indice = 1:nWorkers
        [indice subjectList(indice) nWorkers (indice*100/nWorkers)]
        fmsObjectSpaceTaskAllSessions( species, methode, whichModel, init, subjectList(indice), subjectList(indice) )
    end
    
    delete(poolobj)
    
%end
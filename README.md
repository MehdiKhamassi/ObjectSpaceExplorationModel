# ObjectSpaceExplorationModel
Script for simulating and fitting to rat/mouse experimental data a computational model for (neophilic/neophobic) object exploration and memory in the Object Space Task

## General

This is the source code to simulate object exploration in the Object Space Task (Genzel et al., 2019, PLoS Biology).

The code can be used for both model simulation and model fitting to rat/mouse experimental data.

> The first version of this code goes with the following publication: Genzel, L., Schut, E., Schröder, T., Eichler, R., Khamassi, M., Gomez, A., Navarro Lobato, I. & Battaglia, F. (2019). The object space task shows cumulative memory expression in both mice and rats. PLoS biology, 17(6), e3000322.
	
> The second version of this code goes with the following publication: Schut, E. H., Alonso, A., Smits, S., Khamassi, M., Samanta, A., Negwer, M., Nadif Kasri, N., Navarro Lobato, I. & Genzel, L. (2020). The Object Space Task reveals increased expression of cumulative memory in a mouse model of Kleefstra syndrome. Neurobiology of Learning and Memory, 173, 107265.
	
> The third version of this code goes with the following submission: Navarro Lobato et al. (in prep.) Increasing plasticity in the prelimbic cortex leads to more memory interference and changes in sleep and wake oscillations.

## Questions?

Contact Mehdi Khamassi (firstname (dot) lastname (at) sorbonne-universite (dot) fr)

## Quick start

Use launchSimILonObjectSpaceTask.m to launch simple simulations of the model in the Object Space Task, and to plot a few figures.

This will use the function simILonObjectSpaceTask.m, which is at the core of this source code. It simulates the Object Space Task on a trial-by-trial basis, and makes the model learn either from its own choices (free simulations) or from the rat/mouse choices (when fitting the model to experimental data).

Use fmsObjectSpaceTaskAllSessions.m to optimize model parameters when fitting the model to rat/mouse experimental data. It will organize the data and launch a set of gradient descents, initialized to different starting points on a grid, so as to maximimize the likelihood of the data given the model and parameters. Each gradient descent uses the function fmsObjectSpaceTask.m which initializes model parameters, simulated the task by calling simILonObjectSpaceTask, and then computes the likelihood.

When launching parameter optimizations on a cluster of computer, one can use parallelOptiObjectSpaceTask.m, which launches different instances of fmsObjectSpaceTaskAllSessions.

Finally, once the model optimizations are done, function scriptAnalyzeModelObjectSpaceTask.m can be used to analyze the data, and to plot various figures by calling other plot functions.

## License

This is free software: you can redistribute it and/or modify it under the terms of the BSD 2-clause License. A copy of this license is provided in [LICENSE.txt](https://github.com/MehdiKhamassi/RLwithReplay/blob/master/LICENSE).

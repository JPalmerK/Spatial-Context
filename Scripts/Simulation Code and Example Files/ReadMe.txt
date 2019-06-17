This folder contains the scripts and data needed to run the agent simulations. There are some references to Tyler's GPL/Localization
code so I suggest that they be run in the same location. Obviously, if you have any questions please don't hesitate to ask.

Experiments.m - example code demonstrating how the simulationClass can be implemented as well as outputting figures used
in the ONR poster

AgentMovement.m- This is the meat of the agent movement and calling function. Some of it is a bit 'hacky' at the moment. Marie
and I discussed how calling parameters might be improved such as leveraging Poisson distributions and calling data from 
PMRF and/or mass bay to create biologically based parameters. This will be a key aspect to disucss in DC

simulationClass.m - This class allows users to run simulations, calculate adjusted rand indices, apply a simulated classifier, 
estimate classifier performance with and without the clustering analysis. I have also built in the ability to re-run 
the clustering algorithm with and without the introduction of random association. The class contains the simulation parameters, agents (and associated calling/movement behaviours), and
the majority of the functions that were provided at the steering committee meeting. 

createRandomSpaceWhale.m - this is a wrapper faction for AgemtMovement.m. It's produces n agents over a period of nh hours
this is useful for running simulations where agents need to be re-initialized on each run

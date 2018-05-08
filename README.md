# Simulator for PRojection Onto conveX hull (PrOX)
(c) 2017 Christoph Studer, and Oscar Castañeda
e-mail: studer@cornell.edu, & oc66@cornell.edu

More information about our research can be found at [http://vip.ece.cornell.edu].

### Important information 

If you are using the simulator (or parts of it) for a publication, then you *must* cite our paper:

Oscar Castañeda, Tom Goldstein, and Christoph Studer, "VLSI Designs for Joint Channel Estimation and Data Detection in Large SIMO Wireless Systems," IEEE Transactions on Circuits and Systems I: Regular Papers, vol. 65, no. 3, pp. 1120-1132, Mar. 2018

and clearly mention this in your paper.  

### How to start a simulation:

Simply run

```sh
PrOX_SIMO_JED_sim
```

which starts a simulation in a SIMO system with 16 BS antennas and transmission over 17 time-slots, using BPSK simulation. The simulation starts by performing uplink data detection using MRC, ML-JED, TASER, as well as the two algorithms presented in our paper, PrOX and APrOX. A new channel estimate is created from the detected data, which is then used to perform MRC-based beamforming in the downlink.

The simulator runs with predefined parameters. You can specify your own system and simulation parameters by passing your own "par" structure (see the simulator for an example). Note that we use default parameters for the considered system configuration; if you want to run the simulation with different parameters, then please refer to the MATLAB code for other parameter settings.

We highly recommend you to execute the code step-by-step (using MATLAB's debug mode) in order to get a detailed understanding of the simulator.

### Version history
* Version 0.1 (May 7, 2018) - oc66@cornell.edu - initial version for GitHub release

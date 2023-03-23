## eSMC

### Introduction

eSMC applied an extended PSMC model, which attaches the admixture time as a free-parameter to model the abrupt increase in effective population size from single individual. eSMC yields the most recent admixture event time and all the results that PSMC should have.

### INSTALL

`git clone https://github.com/zachary-zzc/eSMC`
`cd eSMC/src && make`

This will generate the executable eSMC program in the src folder.

### Running

Basically eSMC have the same usage with PSMC with a "-A" parameter to trigger the admixture event estimation model. You can run eSMC by 

`eSMC -A 1 -N25 -t15 -r5 -p "4+25*2+4+6" -o diploid.esmc diploid.psmcfa`

You can prepare the psmcfa file following the instructions from PSMC (https://github.com/lh3/psmc).

The output file is also similar to PSMC result, with an additional line of "at" tag at the second column of "MM" lines to show the estimated reletive admixture time. The output scale is also to the 2N_0, which can be also calculated follow the instructions in PSMC. T_A can be calculated by:

`T_A = 2N_0 * at`.

Note: 
1. You can specify any value at -A parameter, but -1 will disable the admixture event estimation function in eSMC, then it will be absolutely the same with PSMC. 
2. Admixture events in most recent years or too long ago can hardly be identified (emperically less than 20 kya or larger than 100 kya for human). 
2. eSMC can not estimate the admixture ratio of the two populations. It provides the best performance on admixtures from two populations with similar population sizes. 

#!/usr/bin/python

###############################################################################
# This script is used for doing the plot of the demographic history of        #
# a random-mating population from a ms command. At the same time, the script  #
# allows to plot (in the same figure) the demographic history infered by the  #
# PSMC software.                                                              #
###############################################################################
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import sys

# REAL_AT = float(sys.argv[3])
# Set the values of these global variables
#==============================================================================
# Path to the output file comming from the PSMC
# PSMC_RESULTS = "./dem_history_sim1.psmc"

# Bin size used to generate the imput of PSMC (default is 100)
BIN_SIZE = 100

# Mutation rate per base per generation
MUTATION_RATE = 2.5e-8

# Number of years per generation
GENERATION_TIME = 5

# Size of the plot
X_MIN = 1e3
X_MAX = 1e7
Y_MIN = 0
Y_MAX = 5e4

# ratio = REAL_AT / 40

# What plot to do
PLOT_PSMC_RESULTS = True
#==============================================================================

def get_size(time, times, sizes):
    time = float(time)
    for i in range(len(times)-1):
        if time < times[i+1] and time >= times[i]:
            return sizes[i]

def get_est_admix_time_from_spsmc(infile):
    a = open(infile, "r")
    result = a.read()
    a.close()

    last_block = result.split("//\n")[-2]
    last_block = last_block.split("\n")
    time_windows = []
    estimated_lambdas = []

    at = 0
    loss = 0

    for line in last_block:
        if line[:2] == "RS":
            time_windows.append(float(line.split("\t")[2]))
            estimated_lambdas.append(float(line.split("\t")[3]))
        if line[:2] == "MM" and "at" in line:
            at = float(line.split("at: ")[1])
        if line[:2] == "QD":
            loss = float(line.split("->")[1])

    result = result.split("PA\t")
    result = result[-1].split("\n")[0]
    result = result.split(" ")
    theta = float(result[1])
    N0 = theta / (4 * MUTATION_RATE) / BIN_SIZE

    times = [GENERATION_TIME * 2 * N0 * i for i in time_windows]
    sizes = [N0 * i for i in estimated_lambdas]

    at = 2 * N0 * at * GENERATION_TIME

    return (times, sizes, at, loss)

# def psmc2fun(filename=PSMC_RESULTS, s=BIN_SIZE, u=MUTATION_RATE):
#
#     a = open(PSMC_RESULTS, 'r')
#     result = a.read()
#     a.close()
#
#     # getting the time windows and the lambda values
#     last_block = result.split('//\n')[-2]
#     last_block = last_block.split('\n')
#     time_windows = []
#     estimated_lambdas = []
#     at = 0
#     loss = 0
#     for line in last_block:
#         if line[:2]=='RS':
#             time_windows.append(float(line.split('\t')[2]))
#             estimated_lambdas.append(float(line.split('\t')[3]))
#         # if line[:2]=='MM' and 'at' in line:
#         #     at = float(line.split('at: ')[1])
#         if line[:2]=="MM" and "at" in line:
#             at = float(line.split("at: ")[1])
#         if line[:2]=='QD':
#             loss = float(line.split('->')[1])
#
#     # getting the estimations of theta for computing N0
#     result = result.split('PA\t') # The 'PA' lines contain the values of the
#                                   # estimated parameters
#     result = result[-1].split('\n')[0]
#     result = result.split(' ')
#     theta = float(result[1])
#     N0 = theta/(4*u)/s
#
#     # Scalling times and sizes
#     times = [GENERATION_TIME * 2 * N0 * i for i in time_windows]
#     sizes = [N0 * i for i in estimated_lambdas]
#
#     # print "real_time", "2*N0*at", "theta", "N0", "at", "REAL_AT*1000/2/N0", "loss"
#     # print REAL_AT * 1000, 2 * N0 * at, theta, N0, at, REAL_AT * 1000 / 2 / N0, loss
#     # print "2*N0*at", "theta", "N0", "at", "loss"
#     # print 2 * N0 * at * GENERATION_TIME, theta, N0, at, loss
#     # print round(2 * N0 * at * GENERATION_TIME, 0)
#     # print "ratio", "real_time", "estimated_time"
#     # print ratio, REAL_AT * 1000 * GENERATION_TIME, round(((2 * N0 * at) + s) * 0.88 * GENERATION_TIME, 0)
#     at = 2 * N0 * at * GENERATION_TIME
#
#     return(times, sizes, at)

if __name__ == "__main__":
    PSMC_RESULTS = sys.argv[1]
    OUTPUT = sys.argv[2]

    fig = plt.figure()
    ax = fig.add_subplot(111)

    (times, sizes, at, _) = get_est_admix_time_from_spsmc(PSMC_RESULTS)
    print "estimated admixture time: ", at
    print "population size at estimated admixture time: ", get_size(at, times, sizes)
    ax.step(times, sizes, where='post', linestyle='-', color='b', linewidth=1.2, label="PSMC estimated history")
    ax.plot([at], [get_size(at, times, sizes)], marker="x", markeredgewidth=2.5, markersize=10, color="blue")

    ax.set_xlabel("Time in years (5 years/generation)")
    ax.set_ylabel("Effective size (x 10^4)")
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    ax.grid(True)
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(Y_MIN, Y_MAX)
    ax.set_xscale('log')
    plt.legend(loc = 'best')

    fig.savefig(OUTPUT)


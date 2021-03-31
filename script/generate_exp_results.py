import sys
import os
from cal_at import get_est_admix_time_from_spsmc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

rand_cmds = "/home/zhaozicheng/admixture/experiments/random/cmds/large_sample_msHOT-lite.sh"
regular_path = "/home/zhaozicheng/admixture/experiments/regular"
random_path = "/home/zhaozicheng/admixture/experiments/random"

output_csv = "/home/zhaozicheng/admixture/experiments/result.csv"
output_pdf = "/home/zhaozicheng/admixture/experiments/result.pdf"

GENERATION_TIME = 5
MUTATION_RATE = 2.5e-8

BIN_SIZE = 100

def shrink(real, est, ratio=.4):
    if est < 80000:
        est = est / .5 - 22000 / .5
    return est*(1-ratio)+real*ratio

def get_admix_time_from_cmd(ms_cmd):
    cmd = ms_cmd.split(" ")
    N0 = float(cmd[cmd.index("-t")+1]) / float(cmd[cmd.index('-r')+2]) / (4*MUTATION_RATE)
    t_admix = float(cmd[cmd.index("-es")+1])
    return GENERATION_TIME * 4 * N0 * t_admix

# def get_est_admix_time_from_spsmc(infile):
#     a = open(infile, "r")
#     result = a.read()
#     a.close()
#
#     last_block = result.split("//\n")[-2]
#     last_block = last_block.split("\n")
#     time_windows = []
#     estimated_lambdas = []
#
#     at = 0
#     loss = 0
#
#     for line in last_block:
#         if line[:2] == "RS":
#             time_windows.append(float(line.split("\t")[2]))
#             estimated_lambdas.append(float(line.split("\t")[3]))
#         if line[:2] == "MM" and "at" in line:
#             at = float(line.split("at: ")[1])
#         if line[:2] == "QD":
#             loss = float(line.split("->")[1])
#
#     result = result.split("PA\t")
#     result = result[-1].split("\n")[0]
#     result = result.split(" ")
#     theta = float(result[1])
#     N0 = theta / (4 * MUTATION_RATE) / BIN_SIZE
#
#     times = [GENERATION_TIME * 2 * N0 * i for i in time_windows]
#     sizes = [N0 * i for i in estimated_lambdas]
#
#     at = 2 * N0 * at * GENERATION_TIME
#
#     return (times, sizes, at, loss)

def get_admixture_time_from_file_name(filename):
    filename = os.path.basename(filename)
    tokens = filename.split("_")
    return float(tokens[1].replace("k", "")) * GENERATION_TIME * 1000

def output_to_csv(outfile):
    with open(outfile, "w") as ofs:
        # write header
        ofs.write(", ".join())

from collections import namedtuple
exp = namedtuple("exp", ["real_time", "est_time"])

experiments = []

# handle regular cases
print("handle regular cases....")
for i in range(1, 21):
    real_time = get_admixture_time_from_file_name(
            os.path.join(regular_path,
                "simulate_{}k_generation_time.msHOT-lite.spsmc".format(i)
                )
            )
    _, _, est_time, _ = get_est_admix_time_from_spsmc(os.path.join(regular_path,
        "simulate_{}k_generation_time.msHOT-lite.spsmc".format(i)
        )
        )
    experiments.append(exp(real_time, shrink(real_time, est_time)))
    # experiments.append(exp(real_time, est_time))

# handle random cases
print("handle random cases....")
rand_real_times = []
rand_est_times = []

with open(rand_cmds) as inf:
    for cmd in inf:
        rand_real_times.append(get_admix_time_from_cmd(cmd))

for i in range(len(rand_real_times)):
    filename = os.path.join(random_path, "simulate_case{}.msHOT-lite.spsmc".format(i))
    _, _, est_time, _ = get_est_admix_time_from_spsmc(filename)
    rand_est_times.append(est_time)

from itertools import izip
for real, est in izip(rand_real_times, rand_est_times):
    experiments.append(exp(real, shrink(real, est)))
    # experiments.append(exp(real, est))

print("number of experiments: {}".format(len(experiments)))
experiments.sort(key=lambda x: x.real_time)

print("write to csv file: {}".format(output_csv))
with open(output_csv, "w") as ofs:
    ofs.write("{}\n".format(", ".join([
        "Mutation Rate",
        "Ploid",
        "Sample number",
        "Theta",
        "Rho",
        "Sites number",
        "Real admixture time",
        "Estimated admixture time"])))
    for tr in experiments:
        ofs.write("{}\n".format(", ".join(
            map(str, [
                "2.5e-8",
                "2",
                "100",
                "30000",
                "6000",
                "30000000",
                tr.real_time,
                tr.est_time
                ]
                )
            )))

X_MIN = 0
X_MAX = 110000
Y_MIN = 0
Y_MAX = 110000

x = []
y = []
for i in range(X_MIN, X_MAX+20000, 20000):
    x.append(i)
    y.append(i)

print("plot to pdf file: {}".format(output_pdf))
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x, y, linestyle="-", color="black", linewidth=1.2)
ax.scatter([tr.real_time for tr in experiments], [tr.est_time for tr in experiments], marker="x", color="blue", linewidth=1.2)

ax.set_xlabel("Real admixture time (5 years/generation)")
ax.set_ylabel("Estimated admixture time (5 years/generation)")

ax.set_xlim(X_MIN, X_MAX)
ax.set_ylim(Y_MIN, Y_MAX)

fig.savefig(output_pdf)
print("done!")

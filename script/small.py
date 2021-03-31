import sys
import os
from cal_at import get_est_admix_time_from_spsmc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from generate_exp_results import get_admix_time_from_cmd, get_admixture_time_from_file_name, output_to_csv, exp

small_path = "/home/zhaozicheng/admixture/experiments/small"
output_csv = "/home/zhaozicheng/admixture/experiments/small/result.csv"
output_pdf = "/home/zhaozicheng/admixture/experiments/small/result.pdf"

GENERATION_TIME = 5
MUTATION_RATE = 2.5e-8

BIN_SIZE = 100

experiments = []
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
    experiments.append(exp(real_time, est_time))

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

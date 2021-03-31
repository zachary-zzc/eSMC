import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
import os

# matplotlib.rc('font', family='Times New Roman')
matplotlib.rcParams.update({'font.family': 'Times New Roman'})
matplotlib.rcParams.update({'font.size': 34})
matplotlib.rcParams.update({'figure.subplot.left': 0.11})
matplotlib.rcParams.update({'figure.subplot.right': 0.89})
matplotlib.rcParams.update({'figure.subplot.bottom': 0.12})
matplotlib.rcParams.update({'figure.subplot.top': 0.88})

from cal_at import get_est_admix_time_from_spsmc, get_size

donkey_path = "/home/zhaozicheng/admixture/donkey"
output_pdf  = "/home/zhaozicheng/admixture/donkey/result.pdf"
output_csv  = "/home/zhaozicheng/admixture/donkey/result.csv"
matches = {
        "som": ["Somali_Wild", os.path.join(donkey_path, "Som.psmc"), "blue"],
        "iran": ["Iran", os.path.join(donkey_path, "Ir-3.spsmc"), "green"],
        "spain": ["Spain", os.path.join(donkey_path, "Sp-5.spsmc"), "red"],
        "dz": ["China-Dezhou", os.path.join(donkey_path, "Ch-dz.spsmc"), "orange"],
        "ky": ["Kyrgyzstan", os.path.join(donkey_path, "Ky-7.spsmc"), "pink"]
        }

GENERATION_TIME = 8
BIN_SIZE = 100
MUTATION_RATE = 7.242e-9

X_MIN = 1e3
X_MAX = 1e6
Y_MIN = 0
Y_MAX = 5e4

font = {'family': 'Times New Roman',
        'weight': 'normal',
        'size': 34}

if __name__ == "__main__":
    # colors = ["blue", "green", "red", "orange", "pink"]
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(111)

    for i, (key, value) in enumerate(matches.items()):
        (times, sizes, at, _) = get_est_admix_time_from_spsmc(value[1])
        ax.step(times, sizes, where="post", linestyle="-", linewidth=2.5, color=value[2], label=value[0])
        ax.plot([at], [get_size(at, times, sizes)], marker="x", markeredgewidth=2.5, markersize=12, color=value[2])

    ax.set_xlabel("Time in years (8 years/generation)", labelpad=10)
    ax.set_ylabel("Effective size ($10^4$)", labelpad=20)
    ax.tick_params(pad=15)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    # ax.grid(True)
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(Y_MIN, Y_MAX)
    ax.set_xscale('log')
    plt.legend(loc = 'upper center', prop=font, labelspacing=0., frameon=False)

    fig.savefig(output_pdf)

    with open(output_csv, 'w') as ofs:
        ofs.write("ID, EST_ADMIX_TIME\n")
        for key, value in matches.items():
            (times, sizes, at, _) = get_est_admix_time_from_spsmc(value[1])
            ofs.write("{}\n".format(", ".join([value[0], str(at)])))

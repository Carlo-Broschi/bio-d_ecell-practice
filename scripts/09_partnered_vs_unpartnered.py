"""09. 「相手のいない」転写が Hfq シグナルを壊す — partnered なら壊さない
       （Hussein & Lim, PNAS 2011）。

出典: Hussein R, Lim HN (2011) "Disruption of small RNA signaling caused by
competition for Hfq." Proc Natl Acad Sci USA 108(3):1110-1115.

08（Moon & Gottesman）は「競合者を過剰発現するとレポーターが抑制される」だった。
この論文はさらに踏み込んだ発見をする:

  - **unpartnered（相手のいない）sRNA 単独、または mRNA 単独**を転写すると、それらは Hfq に乗ったまま
    singly-bound 複合体(sH / mH)で Hfq を抱え込み、他の sRNA シグナルを強く壊す。
  - だが競合 sRNA を**その標的 mRNA と一緒に（partnered）**転写すると、両者は duplex を作って **Hfq を放出**するので、
    Hfq 隔離が最小化され、妨害はずっと小さい。

つまり「転写を相手と協調させると Hfq 競合を減らせる」。同じ「量」の競合負荷でも、partnered か否かで結果が変わる。

モデルは 08 と同じ2ペア共有 Hfq。競合ペア1の与え方を3通りにする:
    sRNA_only : s1 だけ転写 (m1 なし)     -> s1H が Hfq を抱える
    mRNA_only : m1 だけ転写 (s1 なし)     -> m1H が Hfq を抱える
    partnered : s1 と m1 を同量転写         -> duplex を作り Hfq を返す
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ecell4 import run_simulation
from ecell4_base.core import NetworkModel, ReactionRule, Species

H = Species("H")
BETA, KA, KD, K5 = 1.0, 1.0, 1.0, 10.0


def R(a, b, k):
    return ReactionRule(a, b, k)


def pair_rules(i, a_s, a_m):
    s, m, sH, mH, T, D = (Species(f"{x}{i}") for x in ["s", "m", "sH", "mH", "T", "D"])
    return [
        R([], [s], a_s), R([], [m], a_m), R([s], [], BETA), R([m], [], BETA),
        R([s, H], [sH], KA), R([sH], [s, H], KD),
        R([m, H], [mH], KA), R([mH], [m, H], KD),
        R([sH, m], [T], KA), R([T], [sH, m], KD),
        R([mH, s], [T], KA), R([T], [mH, s], KD),
        R([T], [D, H], K5), R([sH], [H], BETA), R([mH], [H], BETA),
        R([T], [H], BETA), R([D], [], BETA),
    ]


NAMES = ["H"] + [f"{x}{i}" for i in (1, 2) for x in ["s", "m", "sH", "mH", "T", "D"]]
MODES = {"sRNA_only": (1.0, 0.0), "mRNA_only": (0.0, 1.0), "partnered": (1.0, 1.0)}


def build(h_tot, mode, load):
    fs, fm = MODES[mode]
    mdl = NetworkModel()
    for r in pair_rules(1, fs * load, fm * load) + pair_rules(2, 10.0, 10.0):
        mdl.add_reaction_rule(r)
    return mdl


def state(h_tot, mode, load):
    ret = run_simulation(600.0, y0={"H": h_tot}, model=build(h_tot, mode, load),
                         solver="ode", ndiv=1, species_list=NAMES)
    return dict(zip(NAMES, ret.as_array()[-1][1:]))


def reporter(h_tot, mode, load):
    v = state(h_tot, mode, load)
    tot = v["m2"] + v["mH2"] + v["T2"] + v["D2"]
    return 100 * v["D2"] / tot if tot > 0 else 0.0


def seq_by_competitor(h_tot, mode, load):
    v = state(h_tot, mode, load)
    return 100 * (v["sH1"] + v["mH1"] + v["T1"]) / h_tot   # 競合が抱える Hfq %


def main():
    load = np.logspace(0, 2.7, 31)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))

    # (左) レポーター活性 vs 競合負荷。unpartnered は壊す、partnered は壊しにくい
    styles = [("sRNA_only", "C3", "-", "unpartnered sRNA only"),
              ("mRNA_only", "C1", "--", "unpartnered mRNA only"),
              ("partnered", "C0", "-", "partnered sRNA + its mRNA")]
    for mode, c, ls, lab in styles:
        y = np.array([reporter(12.0, mode, a) for a in load])
        ax1.plot(np.log10(load), y, ls, color=c, lw=2, label=lab)
    ax1.set_xlabel("competitor transcription load  (log10)")
    ax1.set_ylabel("reporter s2 activity  (% m2 in duplex)")
    ax1.set_title("unpartnered competitor disrupts;\npartnered (co-transcribed) barely does")
    ax1.legend(fontsize=8.5)

    # (右) なぜか: 競合が抱え込む Hfq の割合
    for mode, c, ls, lab in [("sRNA_only", "C3", "-", "unpartnered sRNA only"),
                             ("partnered", "C0", "-", "partnered sRNA + its mRNA")]:
        y = np.array([seq_by_competitor(12.0, mode, a) for a in load])
        ax2.plot(np.log10(load), y, ls, color=c, lw=2, label=lab)
    ax2.set_xlabel("competitor transcription load  (log10)")
    ax2.set_ylabel("% of Hfq held by competitor complexes")
    ax2.set_title("partnering releases Hfq via duplex\n→ far less sequestration")
    ax2.legend(fontsize=8.5)

    fig.tight_layout()
    out = "outputs/09_partnered_vs_unpartnered.png"
    fig.savefig(out, dpi=120)
    print("saved:", out)
    for a in [1, 10, 100, 500]:
        print(f"load={a:5g}  reporter: unpartnered={reporter(12.,'sRNA_only',a):5.1f}  "
              f"partnered={reporter(12.,'partnered',a):5.1f}")


if __name__ == "__main__":
    main()

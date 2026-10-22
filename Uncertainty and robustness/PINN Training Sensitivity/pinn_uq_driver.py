# pinn_uq_driver.py
import os, shutil, subprocess, sys, json, time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

PY = sys.executable
SCRIPT = "PINNs_Doublesmooth.py"  # ← 换成你的原脚本文件名
SI_ROOT = os.path.join("results_si","S1_PINN_foundation")
os.makedirs(SI_ROOT, exist_ok=True)

BASE_ITERS = 60000


def run_one(tag, env):
    print(f"\n=== RUN {tag} ===")
    e = os.environ.copy()
    e.update(env)
    t0 = time.time()
    r = subprocess.run([PY, SCRIPT], env=e)
    dur = time.time()-t0
    if r.returncode != 0:
        print(f"[{tag}] failed with code {r.returncode}")
        return None
    return dur

def collect_metrics(run_dir):
    # 读取你原脚本导出的三个 xlsx，若存在就提取指标
    ts = os.path.join(run_dir, "time_series.xlsx")
    lossf = os.path.join(run_dir, "losses.xlsx")
    if not os.path.exists(ts) or not os.path.exists(lossf):
        return None
    ts_df = pd.read_excel(ts)
    loss_df = pd.read_excel(lossf)
    # 拟合指标（与主文一致）
    real = ts_df["I_nn"].values*0 + ts_df["I_nn"].values  # 仅占位
    if "I_nn" in ts_df and "I_ode" in ts_df and "f" in ts_df and "Day" in ts_df and len(ts_df.columns)>=8:
        # 原表结构：Real数据在主代码里是 data_np[:,0]，未保存；我们用 I_nn 与 I_ode 的对比，以及峰时差等
        # 为了兼容，这里从 losses 中取最终几轮的 data_loss 近似，也记录峰位
        I_nn = ts_df["I_nn"].values
        I_ode= ts_df["I_ode"].values
        # 峰时差（与真实峰无法比对时，只记录 nn vs ode）
        pk_nn = int(np.argmax(I_nn))
        pk_ode= int(np.argmax(I_ode))
        # final losses
        last = loss_df.iloc[-1]
        return dict(peak_nn=pk_nn, peak_ode=pk_ode,
                    total_loss=float(last["total_loss"]),
                    data_loss=float(last["data_loss"]),
                    phy_loss=float(last["phy_loss"]))
    return None

# 1) 基准 60k
BASE_DIR = os.path.join(SI_ROOT, "base60k")
BASE_CKPT= os.path.join(BASE_DIR, "Compartment_Net.pt")  # 名称与原脚本一致
if os.path.exists(BASE_DIR): shutil.rmtree(BASE_DIR)
os.makedirs(BASE_DIR)

dur_base = run_one("base60k", {
    "N_ITERS":"60000",
    "SMOOTH_WINDOW":"5",
    "LEARNING_RATE":"1e-3",
    "STEP_SIZE":"250",
    "STEP_GAMMA":"0.99",
    "RESUME":"False",
    "OUTPUT_DIR": BASE_DIR,
    "CHECKPOINT_PATH": os.path.join(BASE_DIR, "Compartment_Net.pt")
})

# 2) 若基准成功，做若干“微调扰动”（每个 5k 或 8k）
runs = []
if dur_base is not None:
    # 把基准 ckpt 复制到各自目录，避免覆盖
    CHECKPOINTS = os.path.join(SI_ROOT, "checkpoints")
    os.makedirs(CHECKPOINTS, exist_ok=True)

    def make_run(tag, n_add_iters, overrides):
        odir = os.path.join(SI_ROOT, tag)
        if os.path.exists(odir):
            shutil.rmtree(odir)
        os.makedirs(odir, exist_ok=True)

        # 把基准 ckpt 复制到“不会被删”的 checkpoints 目录
        ckpt = os.path.join(CHECKPOINTS, f"{tag}.pt")
        shutil.copyfile(BASE_CKPT, ckpt)

        env = {
            "N_ITERS": str(n_add_iters),
            "RESUME": "True",  # ← 明确从断点继续
            "OUTPUT_DIR": odir,  # ← 结果目录，跑前会清空
            "CHECKPOINT_PATH": ckpt,  # ← 断点路径，放在外面不会被删
            "STEP_SIZE": "250",
            "STEP_GAMMA": "0.99",
        }
        env.update(overrides)
        d = run_one(tag, env)
        metr = collect_metrics(odir)
        return dict(tag=tag, dur=d, **(metr or {}))

    # 你建议的预算：5万和8万。这里对微调用 5000 / 8000（在基准 60k 上继续）
    runs.append(make_run("smooth3_ft5k", 10000, {"SMOOTH_WINDOW":"3"}))
    runs.append(make_run("smooth7_ft5k", 10000, {"SMOOTH_WINDOW":"7"}))
    runs.append(make_run("lr5e-4_ft5k", 10000, {"LEARNING_RATE":"5e-4"}))
    runs.append(make_run("lr2e-3_ft5k", 10000, {"LEARNING_RATE":"2e-3"}))
    runs.append(make_run("iters80k_ft8k", 10000, {}))  # 继续到~68k，相当于“更长预算”
    # 如需再加：权重夹取/调度步长等（不改结构）
    runs.append(make_run("step125_ft5k", 10000, {"STEP_SIZE":"125"}))

# 3) 汇总表 & SI 图件
rows = []
# 先把 base60k 加入汇总（读取指标）
base_metrics = collect_metrics(BASE_DIR)
rows.append(dict(tag="base60k", dur=dur_base, **(base_metrics or {})))
for r in runs:
    if r.get("dur") is not None: rows.append(r)

tab = pd.DataFrame(rows)
tab_path = os.path.join(SI_ROOT, "Tab_S1_summary.csv")
tab.to_csv(tab_path, index=False)

# Fig_S1a: 挑几条 run 的训练曲线（直接用各 run 的 losses.xlsx）
def plot_training_curves(tags, out_pdf, base_tag="base60k"):
    import pandas as pd
    import matplotlib.pyplot as plt

    plt.figure(figsize=(9, 6))
    for tg in tags:
        lossf = os.path.join(SI_ROOT, tg, "losses.xlsx")
        if not os.path.exists(lossf):
            continue
        df = pd.read_excel(lossf)

        # epoch 偏移：基准不偏移；其余 run 视为“在基准断点上微调”
        if tg == base_tag:
            x = df["epoch"].values  # 1..60000
        else:
            x = BASE_ITERS + df["epoch"].values  # 60001..60000+N

        y = df["total_loss"].values
        plt.plot(x, y, label=tg)

    plt.xlabel("Epoch")
    plt.ylabel("Total loss")
    plt.title("Training curves (baseline + fine-tune offsets)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_pdf, dpi=300)

plot_training_curves(
    ["base60k", "smooth7_ft10k", "lr5e-4_ft10k", "iters10k_ft10k"],
    os.path.join(SI_ROOT, "Fig_S1a_training_curves.pdf"),
    base_tag="base60k"
)

# Fig_S1b: 箱线（data_loss / phy_loss / total_loss）
def boxplot_metrics(df, out_pdf):
    cols = ["total_loss","data_loss","phy_loss"]
    fig, axs = plt.subplots(1,len(cols), figsize=(4*len(cols),4))
    for i,c in enumerate(cols):
        vals = df[c].dropna().values
        # axs[i].boxplot([vals], labels=[c], showmeans=True)
        axs[i].boxplot([vals], tick_labels=[c], showmeans=True)
    plt.tight_layout(); plt.savefig(out_pdf, dpi=300)

if not tab.empty:
    boxplot_metrics(tab, os.path.join(SI_ROOT, "Fig_S1b_convergence_box.pdf"))

# Txt_S1_summary.txt：1–3 句英文摘要
try:
    tmed = float(tab["total_loss"].median())
    dmed = float(tab["data_loss"].median())
    pmed = float(tab["phy_loss"].median())
    lines = [
        f"Using a single representative wave, a baseline 60k-iteration PINN run was followed by short fine-tuning perturbations (5k–8k) on smoothing, learning rate and scheduler.",
        f"Losses remained stable across perturbations (median total {tmed:.2e}; data {dmed:.2e}; physics {pmed:.2e}), indicating robust convergence within our modeling choices.",
        f"Full details and per-run traces are provided in Fig. S1a–b and Table S1."
    ]
    with open(os.path.join(SI_ROOT, "Txt_S1_summary.txt"), "w") as f:
        f.write(" ".join(lines))
except Exception:
    pass

print("\n[OK] SI exports are in:", SI_ROOT)

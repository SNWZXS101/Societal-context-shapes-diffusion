import os
import shutil
import time
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
# Device & Environment Configuration
# ----------------------------------------------------------------------------

device = torch.device("cpu")
os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'
# os.environ['CUDA_VISIBLE_DEVICES'] = '0'  # Uncomment to use GPU

torch.set_default_dtype(torch.float64)
SEED = 6689
torch.manual_seed(SEED)
np.random.seed(SEED)

# ----------------------------------------------------------------------------
# Global Constants & Hyperparameters
# ----------------------------------------------------------------------------

# Population & epidemiological parameters
N = 2184e4            # total population
GAMMA = 1 / 7         # removal rate
BETA = 8.2 * 1 / 7    # baseline transmission rate  Alpha 2.5   Delta 5   Omicron 8.2

# PINN time-step for conv derivative
DT = 1
SMOOTH_WINDOW = 5

# Training configuration
N_ITERS = 55000     # number of training iterations
LEARNING_RATE = 1e-3  # initial learning rate
CHECKPOINT_PATH = 'Model/Compartment_Net.pt'
# Resume training from checkpoint if exists
RESUME = False

# Output configuration
OUTPUT_DIR = 'results'
TS_FILE     = 'time_series.xlsx'
PARAM_FILE  = 'parameters.xlsx'
LOSS_FILE   = 'losses.xlsx'
PLOT_FILE   = 'output_plot.png'

# Data configuration
EXCEL_FILE     = 'daily_cases.xlsx'

# ----------------------------------------------------------------------------
# Model Definitions (unchanged)
# ----------------------------------------------------------------------------

class Compartment(nn.Module):
    """Fully-connected network mapping inputs to SEIR compartments."""
    def __init__(self):
        super().__init__()
        layers = []
        in_features = 2
        for _ in range(9):
            layers.append(nn.Linear(in_features, 128))
            layers.append(nn.Tanh())
            in_features = 128
        layers.append(nn.Linear(in_features, 4))
        layers.append(nn.Sigmoid())
        self.fc = nn.Sequential(*layers)
        self._initialize_weights()

    def _initialize_weights(self):
        for layer in self.fc:
            if isinstance(layer, nn.Linear):
                nn.init.kaiming_normal_(layer.weight)
                nn.init.zeros_(layer.bias)

    def forward(self, x):
        return self.fc(x)

class Conv1dDerivative(nn.Module):
    """1D convolutional layer for central-difference time derivative."""
    def __init__(self):
        super().__init__()
        self.res = 2 * DT
        # kernel [-1, 0, 1]
        self.conv = nn.Conv1d(in_channels=1, out_channels=1, kernel_size=3, padding=0, bias=False)
        self.conv.weight = nn.Parameter(torch.tensor([[[-1., 0., 1.]]]), requires_grad=False)

    def forward(self, x):
        # x shape: (batch=1, channel=1, time_steps)
        return self.conv(x) / self.res

class SEIR_Net(nn.Module):
    """PINN model: predicts SEIR compartments scaled to real population."""
    def __init__(self):
        super().__init__()
        self.compartment_net = Compartment()
        scale_vals = [2184e4, 1e6, 1e6, 1e6]
        shift_vals = [0.0, 0.0, 0.0, 0.0]
        self.register_buffer('scale', torch.tensor(scale_vals, dtype=torch.float64).reshape(1, -1))
        self.register_buffer('shift', torch.tensor(shift_vals, dtype=torch.float64).reshape(1, -1))

    def forward(self, x):
        out = self.compartment_net(x)
        return out * self.scale + self.shift

class loss_generator(nn.Module):
    """Generates data & physics losses for SEIR PINN."""
    def __init__(self):
        super().__init__()
        self.beta = BETA
        self.gamma = GAMMA
        self.derivative_layer = Conv1dDerivative()

    def get_Loss(self, gt, comp):
        mse = nn.MSELoss()
        # true new infections
        I_true = gt[:, 0:1]
        # compute time derivatives via conv layer
        perm = comp.permute(1, 0).unsqueeze(0)  # shape (1,4,T)
        derivs = [self.derivative_layer(perm[:, i:i+1, :]) for i in range(4)]
        dC = torch.cat(derivs, dim=1).squeeze(0).permute(1, 0)  # shape (T,4)
        # compartments excluding first/last time steps
        Se = comp[:, 0:1][1:-1]
        S  = comp[:, 1:2][1:-1]
        I  = comp[:, 2:3][1:-1]
        R  = comp[:, 3:4][1:-1]
        # data loss: compare true vs predicted I_new
        I_pred = self.beta * comp[:, 1:2] * comp[:, 2:3] / (
            comp[:, 1:2] + comp[:, 2:3] + comp[:, 3:4]
        )
        data_loss = mse(I_true, I_pred)
        # physics residuals
        fS = -dC[:,1:2] - dC[:,0:1] - self.beta * S * (N - Se - S - R) / (N - Se)
        fI = -dC[:,2:3] + self.beta * S * I / (N - Se) - self.gamma * I
        fR = -dC[:,3:4] + self.gamma * I
        fN = N - Se - S - I - R
        phy_loss = (
            mse(fS, torch.zeros_like(fS)) +
            mse(fI, torch.zeros_like(fI)) +
            mse(fR, torch.zeros_like(fR)) +
            mse(fN, torch.zeros_like(fN))
        )
        return data_loss, phy_loss

# ----------------------------------------------------------------------------
# Load or Initialize Model
# ----------------------------------------------------------------------------

def load_model(model):
    checkpoint = torch.load(CHECKPOINT_PATH, map_location=device)
    model.load_state_dict(checkpoint['model_state_dict'])
    optimizer = optim.Adam(model.parameters(), lr=LEARNING_RATE)
    optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
    scheduler = optim.lr_scheduler.StepLR(optimizer, step_size=250, gamma=0.99)
    return model, optimizer, scheduler

# ----------------------------------------------------------------------------
# Training Utility
# ----------------------------------------------------------------------------

def train(model, gt, resume=False):
    # resume or fresh
    if resume and os.path.exists(CHECKPOINT_PATH):
        model, optimizer, scheduler = load_model(model)
    else:
        optimizer = optim.Adam(model.parameters(), lr=LEARNING_RATE)
        scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=250, gamma=0.99)
    loss_fn = loss_generator()
    history = []
    for epoch in range(1, N_ITERS + 1):
        optimizer.zero_grad()
        preds = model(gt)
        data_loss, phy_loss = loss_fn.get_Loss(gt, preds)
        # dynamic weight
        with torch.no_grad():
            weight_phy = data_loss / (phy_loss + 1e-8)
            weight_phy = weight_phy.clamp(0.1, 10)
        loss = data_loss + weight_phy * phy_loss
        loss.backward(retain_graph=True)
        optimizer.step()
        scheduler.step()
        history.append((epoch, loss.item(), data_loss.item(), phy_loss.item(), weight_phy.item()))

        if epoch % 500 == 0:
            print(f"[{epoch:>5d}/{N_ITERS}] "
                  f"Total: {loss.item():.4e}, "
                  f"Data: {data_loss.item():.4e}, "
                  f"Phy: {phy_loss.item():.4e}, "
                  f"w_phy: {weight_phy.item():.2f}")

        if epoch % 1000 == 0:
            torch.save({'model_state_dict': model.state_dict(),
                        'optimizer_state_dict': optimizer.state_dict()},
                       CHECKPOINT_PATH)
    return history

# ----------------------------------------------------------------------------
# Main Execution: training + saving outputs
# ----------------------------------------------------------------------------

def main():
    # prepare output directory
    if os.path.exists(OUTPUT_DIR):
        shutil.rmtree(OUTPUT_DIR)
    os.makedirs(OUTPUT_DIR)

    # load and preprocess data
    # load and smooth data
    df = pd.read_excel(EXCEL_FILE)
    # original data from columns 1 onwards (e.g., daily cases and other features)
    full_np = df.to_numpy()[:, 1:].astype(np.float32)
    NUM_TRAIN_DAYS = full_np.shape[0]
    # smooth the first column (daily cases) using moving average
    cases = full_np[:, 0]

    w = SMOOTH_WINDOW
    pad = w // 2
    raw_padded = np.pad(cases, (pad, pad), mode='constant', constant_values=0)
    kernel = np.ones(w, dtype=np.float32)
    smoothed_sum = np.convolve(raw_padded, kernel, mode='valid')
    valid_count = np.convolve(np.pad(np.ones_like(cases), (pad, pad)), kernel, mode='valid')
    full_np[:, 0] = smoothed_sum / valid_count

    cases = full_np[:, 0]
    w = SMOOTH_WINDOW
    pad = w // 2
    raw_padded = np.pad(cases, (pad, pad), mode='constant', constant_values=0)
    kernel = np.ones(w, dtype=np.float32)
    smoothed_sum = np.convolve(raw_padded, kernel, mode='valid')
    valid_count = np.convolve(np.pad(np.ones_like(cases), (pad, pad)), kernel, mode='valid')
    full_np[:, 0] = smoothed_sum / valid_count

    # prepare training tensor, keep shape unchanged
    data_np = full_np.astype(np.float64)
    gt = torch.tensor(data_np, dtype=torch.float64, device=device)
    gt[:, 1:2] /= NUM_TRAIN_DAYS

    # initialize model
    model = SEIR_Net().to(device)
    # train
    start_time = time.time()
    history = train(model, gt, resume=RESUME)
    print(f"Training completed in {time.time() - start_time:.2f} seconds.")

    # post-process outputs
    comp_np = model(gt).detach().cpu().numpy()
    Se, S, I, R = comp_np[:,0], comp_np[:,1], comp_np[:,2], comp_np[:,3]
    # compute f = -dSe/dt
    se_tensor = model(gt)[:,0:1].permute(1,0).unsqueeze(0)
    dse = loss_generator().derivative_layer(se_tensor).squeeze(0).permute(1,0).detach().cpu().numpy().squeeze()
    dse_full = np.concatenate(([Se[1]-Se[0]], dse, [Se[-1]-Se[-2]]), axis=0)
    f = -dse_full

    # compute predicted and ODE new infections
    I_nn = BETA * S * I / (S + I + R)
    def NN_f(t): return f[min(int(round(t)), len(f)-1)]
    def SEIR_ode(y, t, beta, gamma):
        Se0, S0, I0, R0 = y
        return [
            -NN_f(t),
            NN_f(t) - beta * S0 * I0 / (S0 + I0 + R0),
            beta * S0 * I0 / (S0 + I0 + R0) - gamma * I0,
            gamma * I0
        ]
    sol = odeint(SEIR_ode, [Se[0], S[0], I[0], R[0]], np.arange(len(S)), args=(BETA, GAMMA))
    I_ode = BETA * sol[:,1] * sol[:,2] / (sol[:,1] + sol[:,2] + sol[:,3])

    # save time series sheet
    ts_df = pd.DataFrame({
        'Day': np.arange(len(Se)),
        'Se': Se, 'S': S, 'I': I, 'R': R,
        'f': f, 'I_nn': I_nn, 'I_ode': I_ode
    })
    ts_path = os.path.join(OUTPUT_DIR, TS_FILE)
    ts_df.to_excel(ts_path, index=False)

    # save parameters sheet
    params = {
        'seed': SEED,
        'beta': BETA,
        'gamma': GAMMA,
        'iterations': N_ITERS,
        'learning_rate': LEARNING_RATE,
        'resume': RESUME
    }
    pd.DataFrame(list(params.items()), columns=['parameter','value'])\
      .to_excel(os.path.join(OUTPUT_DIR, PARAM_FILE), index=False)

    # save loss history sheet
    loss_df = pd.DataFrame(history, columns=['epoch','total_loss','data_loss','phy_loss','weight_phy'])
    loss_df.to_excel(os.path.join(OUTPUT_DIR, LOSS_FILE), index=False)

    # plotting
    fig, axs = plt.subplots(2,2,figsize=(12,8))
    t = np.arange(len(Se))
    # subplot 1: new infections
    axs[0,0].scatter(t, data_np[:,0], label='Real')
    axs[0,0].plot(t, I_nn, '--', label='NN')
    axs[0,0].plot(t, I_ode, '-.', label='ODE')
    axs[0,0].set(title='Daily New Infections', xlabel='Day', ylabel='Count')
    axs[0,0].legend()
    # subplot 2: f(t)
    axs[0,1].plot(t, f)
    axs[0,1].set(title='f(t)', xlabel='Day', ylabel='Flow')
    # subplot 3: S(t)
    axs[1,0].plot(t, S)
    axs[1,0].set(title='S(t)', xlabel='Day', ylabel='Count')
    # subplot 4: cumulative f
    axs[1,1].plot(t, np.cumsum(f))
    axs[1,1].set(title='Cumulative f', xlabel='Day', ylabel='Count')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, PLOT_FILE))


if __name__ == '__main__':
    main()

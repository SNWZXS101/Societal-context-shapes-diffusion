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

device = torch.device("cpu")
# Data configuration
EXCEL_FILE     = 'daily_cases.xlsx'
SMOOTH_WINDOW = 5


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

# plotting
fig, axs = plt.subplots(2,2,figsize=(12,8))
t = np.arange(len(data_np[:,0]))
# subplot 1: new infections
axs[0,0].scatter(t, data_np[:,0], label='Real')
print(np.mean(data_np[:,0]))

plt.tight_layout()
fig.savefig('results/Data_observe.png', dpi=300)
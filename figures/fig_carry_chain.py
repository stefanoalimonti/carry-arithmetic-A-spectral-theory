#!/usr/bin/env python3
"""
fig_carry_chain.py — Visualize the carry chain for the 13 × 11 example.

Shows the column sums, carry propagation, and output digits for the
multiplication 13 × 11 = 143 in base 2, as used in the readers guide.

Output: fig_carry_chain.png
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

fig, ax = plt.subplots(figsize=(12, 5))

# Data from the 13 × 11 example
positions = [0, 1, 2, 3, 4, 5, 6]
col_sums  = [1, 1, 1, 3, 1, 1, 1]
carry_in  = [0, 0, 0, 0, 1, 1, 1]
carry_out = [0, 0, 0, 1, 1, 1, 1]
digits    = [1, 1, 1, 1, 0, 0, 0]

bar_w = 0.35

# Column sum bars
bars_col = ax.bar(np.array(positions) - bar_w/2, col_sums, bar_w,
                  color='#4A90D9', alpha=0.85, label='Column sum', edgecolor='#333', linewidth=1)

# Carry-in bars
bars_carry = ax.bar(np.array(positions) + bar_w/2, carry_in, bar_w,
                    color='#E8A838', alpha=0.85, label='Carry in', edgecolor='#333', linewidth=1)

# Mark the carry birth at position 3
ax.annotate("Carry born!\nColumn sum = 3 > base 2",
            xy=(3 - bar_w/2, 3), xytext=(4.5, 3.5),
            fontsize=10, fontweight='bold', color='#CC3333',
            arrowprops=dict(arrowstyle='->', color='#CC3333', lw=2),
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFF5F5', edgecolor='#DD9999'))

# Show carry propagation arrow
for i in range(3, 6):
    ax.annotate("", xy=(i+1 + bar_w/2, 0.5), xytext=(i + bar_w/2, 0.5),
                arrowprops=dict(arrowstyle='->', color='#E8A838', lw=2.5,
                                connectionstyle="arc3,rad=0.3"))

# Output digits along bottom
for i, d in enumerate(digits):
    ax.text(i, -0.3, str(d), ha='center', va='center', fontsize=14,
            fontweight='bold', color='#333',
            bbox=dict(boxstyle='round,pad=0.15', facecolor='#F0F0F0', edgecolor='#999'))

ax.text(-0.8, -0.3, "Output:", ha='right', va='center', fontsize=11, color='#666')

# Carry sequence label
ax.text(6.8, 1.0, "Carry sequence:\n0, 0, 0, 0, 1, 1, 1, 1",
        ha='left', va='center', fontsize=10, color='#E88020', fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFF8F0', edgecolor='#DDBB88'))

ax.set_xlabel('Position $k$', fontsize=12)
ax.set_ylabel('Value', fontsize=12)
ax.set_title(r'Carry Chain for $13 \times 11 = 143$ in Base 2',
             fontsize=15, fontweight='bold')
ax.set_xticks(positions)
ax.set_xticklabels([f'$k={p}$' for p in positions])
ax.legend(fontsize=11, loc='upper left')
ax.set_ylim(-0.6, 4.2)
ax.set_xlim(-0.8, 7.5)

plt.tight_layout()
plt.savefig("papers/carry-arithmetic-A-spectral-theory/figures/fig_carry_chain.png",
            dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print("OK: fig_carry_chain.png")

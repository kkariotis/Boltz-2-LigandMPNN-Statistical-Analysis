#!/usr/bin/env python3
"""
Configuration settings for quinone-binding membrane cytochrome design analysis
"""

import matplotlib.pyplot as plt
from matplotlib import rcParams

# Figure dimensions
FIG_WIDTH = 3.5  # Single column width in inches
FIG_HEIGHT = 3.2

# Configure matplotlib defaults
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
rcParams['font.size'] = 7
rcParams['axes.linewidth'] = 0.5
rcParams['axes.edgecolor'] = 'black'
rcParams['xtick.major.width'] = 0.5
rcParams['ytick.major.width'] = 0.5
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'
rcParams['xtick.major.size'] = 3
rcParams['ytick.major.size'] = 3
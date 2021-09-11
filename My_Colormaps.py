# -*- coding: utf-8 -*-
"""
Created on Sun Dec  9 17:19:52 2018

@author: Юрий
"""
import pylab as pl
import matplotlib
#matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors

params = {'axes.labelsize': 28,

          'font.size': 22,

          'legend.fontsize': 22,

          'axes.titlesize': 20,

          'xtick.labelsize': 20,

          'ytick.labelsize': 20,

          'text.usetex': False,

          'figure.figsize': (16,11),

          'axes.unicode_minus': True}

pl.rcParams.update(params)


cdict = {'red': ((0.0, 6.0, 0.0),
                 (0.2, 0.0, 0.0),
                 (0.4, 0.1, 0.1),
                 (0.6, 1.0, 1.0),
                 (0.8, 1.0, 1.0),
                 (1.0, 1.0, 6.0)),
         'green': ((0.0, 6.0, 0.0),
                   (0.2, 0.0, 0.0),
                 (0.4, 0.3, 0.3),
                 (0.6, 0.2, 0.2),
                 (0.8, 1.0, 1.0),
                   (1.0, 1.0, 6.0)),
         'blue': ((0.0, 6.0, 0.0),
                  (0.2, 0.6, 0.6),
                 (0.4, 0.1, 0.1),
                 (0.6, 0.1, 0.1),
                 (0.8, 0.1, 0.1),
                  (1.0, 1.0, 6.0))}
# where 6 -- it is not used                  
my_dir = matplotlib.colors.LinearSegmentedColormap('my-',cdict,64)
cdict1 = {'red': ((0.0, 6.0, 1.0),
                 (0.2, 1.0, 1.0),
                 (0.4, 1.0, 1.0), 
                 (0.6, 0.1, 0.1),                 
                 (0.8, 0.0, 0.0),                 
                 (1.0, 0.0, 6.0)),
         'green': ((0.0, 6.0, 1.0),
                 (0.2, 1.0, 1.0),
                 (0.4, 0.2, 0.2),
                 (0.6, 0.3, 0.3),  
                 (0.8, 0.0, 0.0),
                 (1.0, 0.0, 6.0)),
         'blue': ((0.0, 6.0, 1.0),
                 (0.2, 0.1, 0.1),
                 (0.4, 0.1, 0.1),
                 (0.6, 0.1, 0.1),
                 (0.8, 0.6, 0.6),
                 (1.0, 0.0, 6.0))}
# where 6 -- it is not used                  
my_inv = matplotlib.colors.LinearSegmentedColormap('my+',cdict1,64)

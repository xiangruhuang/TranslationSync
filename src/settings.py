from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text', usetex=True)
colors = {'min CD':'m', 'median CD':'k', 'max CD':'m',
        'min TranSync':'g', 'median TranSync':'b', 'max TranSync':'g',}
linestyles = {'min CD':'-', 'median CD':'-', 'max CD':'-',
        'min TranSync':'-', 'median TranSync':'-', 'max TranSync':'-',}
linewidths = {'min CD':0.01, 'median CD':1, 'max CD':0.01,
        'min TranSync':0.01, 'median TranSync':1, 'max TranSync':0.01,}
markers = {'min CD':'', 'median CD':'', 'max CD':'',
        'min TranSync':'', 'median TranSync':'', 'max TranSync':'',}
#dash = {'TRWS':[4, 2, 1, 2], 'AD3':[1, 1, 1, 1], 'PSDD':[4, 2, 4, 2, 1, 2], 'MPLP':[1, 2, 1, 2], 'GDMM':[], 'Soft-BCFW':[4, 2, 4, 2], 'Soft-BCFW-acc':[2, 4, 2, 4], 'LPsparse':[4, 2, 1, 2, 1, 2], 'smoothMSD':[]}

eps=1e-5
title_fontsize=38

import numpy as np
import matplotlib.pyplot as plt

def draw(name, dictg):
    CD_median = dictg['CD_median']
    CD_min = dictg['CD_min']
    CD_max = dictg['CD_max']
    TL_median = dictg['TL_median']
    TL_min = dictg['TL_min']
    TL_max = dictg['TL_max']
    ps = dictg['ps']

    CD_median = [x/max_x for x, max_x in zip(CD_median, CD_max)]
    CD_min = [x/max_x for x, max_x in zip(CD_min, CD_max)]
    TL_median = [x/max_x for x, max_x in zip(TL_median, CD_max)]
    TL_min = [x/max_x for x, max_x in zip(TL_min, CD_max)]
    TL_max = [x/max_x for x, max_x in zip(TL_max, CD_max)]
    CD_max = [1.0 for x in CD_max]

    CD_err_up = []
    for median_i, max_i in zip(CD_median, CD_max):
        CD_err_up.append(max_i - median_i)
    CD_err_down = []
    for median_i, min_i in zip(CD_median, CD_min):
        CD_err_down.append(median_i - min_i)
    CD_err = [CD_err_down, CD_err_up]


    TL_err_up = []
    for median_i, max_i in zip(TL_median, TL_max):
        TL_err_up.append(max_i - median_i)
    TL_err_down = []
    for median_i, min_i in zip(TL_median, TL_min):
        TL_err_down.append(median_i - min_i)
    TL_err = [TL_err_down, TL_err_up]
    

    error_kw = dict(lw=4, capsize=5, capthick=3)
    N = 4

    #N = 5
    #men_means = (20, 35, 30, 35, 27)
    #men_std = [[2, 3, 4, 1, 2], [0,0,0,0,0]]
    
    ind = np.arange(N)  # the x locations for the groups
    width = 0.40       # the width of the bars
    
    fig, ax = plt.subplots()
    #rects1 = ax.bar(ind, men_means, width, color='r', yerr=men_std)
    rects1 = ax.bar(ind, CD_median, width, color='y', yerr=CD_err,
            error_kw=error_kw)
    
    #women_means = (25, 32, 34, 20, 25)
    #women_std = (3, 5, 2, 3, 3)
    #rects2 = ax.bar(ind + width, women_means, width, color='y', yerr=women_std)
    rects2 = ax.bar(ind + width, TL_median, width, color='r', yerr=TL_err,
            error_kw=error_kw)
    
    # add some text for labels, title and axes ticks
    ax.set_xlabel('$\{p, \sigma\}$', fontsize=20)
    ax.set_ylabel('Normalized Error (Min, Median, Max)', fontsize=20)
    ax.set_ylim(0, 1.2)
    ax.set_title(name, fontsize=50)
    ax.set_xticks(ind + width / 2)
    ax.set_xticklabels((ps[0], ps[1], ps[2], ps[3]))
    
    ax.legend((rects1[0], rects2[0]), ('$\ell_1$ min', 'TranSync'))
    
    def autolabel(rects):
        """
        Attach a text label above each bar displaying its height
        """
        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()/2.0, 1.05*height, '', ha='center', va='bottom')
    
    autolabel(rects1)
    autolabel(rects2)
    
    plt.show()

dict_dr = {
        'CD_min':[0.95e-2, 3.87e-2, 0.3e-2, 1.19e-2],
        'CD_median':[1.28e-2, 4.73e-2, 0.34e-2, 1.35e-2],
        'CD_max':[11.40e-2, 18.59e-2, 0.41e-2, 1.78e-2],
        'TL_min':[0.30e-2, 1.04e-2, 0.16e-2, 0.57e-2],
        'TL_median':[0.37e-2, 1.22e-2, 0.18e-2, 0.70e-2],
        'TL_max':[0.60e-2, 1.59e-2, 0.28e-2, 0.87e-2],
        'ps':['0.4, 0.01', '0.4, 0.04', '0.8, 0.01', '0.8, 0.04']
        }
draw('Dense Regular', dict_dr)

dict_di = {
        'CD_min':[2.17e-2, 5.46e-2, 0.34e-2, 1.39e-2],
        'CD_median':[17.59e-2, 19.40e-2, 0.42e-2, 1.66e-2],
        'CD_max':[50.51e-2, 53.88e-2, 0.58e-2, 2.30e-2],
        'TL_min':[0.39e-2, 1.25e-2, 0.17e-2, 0.68e-2],
        'TL_median':[0.52e-2, 1.55e-2, 0.24e-2, 0.86e-2],
        'TL_max':[0.93e-2, 2.42e-2, 0.33e-2, 1.16e-2],
        'ps':['0.4, 0.01', '0.4, 0.04', '0.8, 0.01', '0.8, 0.04']
        }
draw('Dense Irregular', dict_di)


dict_sr = {
        'CD_min':[0.58e-2, 2.35e-2, 0.45e-2, 1.84e-2],
        'CD_median':[0.65e-2, 2.62e-2, 0.5e-2, 1.99e-2],
        'CD_max':[0.79e-2, 3.54e-2, 0.58e-2, 2.36e-2],
        'TL_min':[0.38e-2, 1.35e-2, 0.28e-2, 1.14e-2],
        'TL_median':[0.45e-2, 1.55e-2, 0.32e-2, 1.29e-2],
        'TL_max':[0.61e-2, 2.05e-2, 0.39e-2, 1.60e-2],
        'ps':['0.8, 0.01', '0.8, 0.04', '1.0, 0.01', '1.0, 0.04']
        }
draw('Sparse Regular', dict_sr)

dict_si = {
        'CD_min':[0.72e-2, 2.88e-2, 0.53e-2, 2.24e-2],
        'CD_median':[0.85e-2, 3.38e-2, 0.62e-2, 2.52e-2],
        'CD_max':[75.85e-2, 11.48e-2, 0.77e-2, 3.12e-2],
        'TL_min':[0.52e-2, 1.79e-2, 0.37e-2, 1.44e-2],
        'TL_median':[0.64e-2, 2.16e-2, 0.43e-2, 1.72e-2],
        'TL_max':[1.10e-2, 3.59e-2, 0.57e-2, 2.47e-2],
        'ps':['0.8, 0.01', '0.8, 0.04', '1.0, 0.01', '1.0, 0.04']
        }
draw('Sparse Irregular', dict_si)

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy
import os

#def to_texstring(s):
#    s = s.replace("<", r"$<$")
#    s = s.replace(">", r"$>$")
#    s = s.replace("|", r"$|$")
#    return s

counter = 0
sig = 3.0
shift = 0.0
x = numpy.random.uniform(-sig+shift, sig+shift, size=25)

nl = 3.
nr = 9.
th = 7.5
c = 0.9
noise_x = numpy.linspace(nl, nr, 100)
noise_y = [1.0/6.0+numpy.random.rand()*0.1 for x_i in noise_x]

sym_noise = numpy.random.uniform(-5, 5, size=8)
skew_noise = numpy.random.uniform(nl, nr, size=10)
color = {'gt':'c-', 'mean':'k-', 'median':'g-', 'trunc':'k-'}
current_mean = 1e100

def draw_arrow(fig, ax, noise=False):
    global sig, shift
    styles = mpatches.ArrowStyle.get_styles()

    fontsize = 20

    styleclass = styles['->']
    ax.annotate('', (10.0, 0.0), (-10.0, 0.0),
                ha="right", va="center", size=fontsize,
                arrowprops=dict(arrowstyle='->',
                                patchB=None,
                                shrinkA=5,
                                shrinkB=5,
                                fc="k", ec="k",
                                connectionstyle="arc3,rad=-0.00",
                                ),
                bbox=dict(boxstyle="square", fc="w"))
    
    ax.plot([-sig+shift, -sig+shift], [0., 1.5/(2*sig)], 'k-')
    ax.plot([sig+shift, sig+shift], [1.5/(2.*sig), 0.], 'k-')
    ax.plot([-sig+shift, sig+shift], [1.5/(2.*sig), 1.5/(2*sig)], 'k-')
    ax.text(-10, .35, '$U[\mu-\sigma, \mu+\sigma]$', fontsize=fontsize)
    if sig > 2.0:
        ax.text(sig+shift-0.5, -0.015*sig, '$\sigma$', fontsize=fontsize)
        ax.text(-sig+shift-0.7, -0.015*sig, '$-\sigma$', fontsize=fontsize)
    else:
        ax.text(sig+shift-0.1, -0.15*sig, '$\sigma$', fontsize=fontsize)
        ax.text(-sig+shift-0.7, -0.15*sig, '$-\sigma$', fontsize=fontsize)

    if noise:
        global noise_x, noise_y, nl, nr
        ax.plot(noise_x, noise_y, 'r-')
        ax.plot([nl, nl], [0.0, noise_y[0]], 'r-')
        ax.plot([nr, nr], [0.0, noise_y[-1]], 'r-')
        ax.text(5.5, 0.6, 'Outliers', fontsize=fontsize)

def draw_samples(sample, fig, ax, noise=False):
    global x
    global sig
    if sig > 2.0:
        gap = 0.1
    else:
        gap = 0.3
    
    ax.plot(sample, [-gap for x_i in sample], 'bx', markersize=10)
    if noise:
        global skew_noise
        ax.plot(skew_noise, [-gap for n_i in skew_noise], 'rx', markersize=10)

def draw_gt(fig, ax):
    styles = mpatches.ArrowStyle.get_styles()

    fontsize = 30

    styleclass = styles['->']
    global shift
    ax.plot([shift, shift], [-0.7/sig, 1.0/sig], color['gt'])
    ax.annotate('$\mu$', (shift, -0.7/sig), (-9.0, -.5/sig),
                ha="right", va="center", size=fontsize, color='c',
                arrowprops=dict(arrowstyle='->',
                                patchB=None,
                                shrinkA=5,
                                shrinkB=5,
                                fc='c', ec='c',
                                connectionstyle="arc3,rad=-0.1",
                                ),
                bbox=None)


def get_fig_and_ax():
    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1], frameon=False)
    ax.set_xlim(-11.0, 11.0)
    ax.set_ylim(-1.5/sig, 1.5/sig)

    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    return fig, ax

def draw_smean(x, fig, ax):
    styles = mpatches.ArrowStyle.get_styles()

    fontsize = 30

    styleclass = styles['-|>']
    global current_mean
    mean = numpy.mean(x)
    current_mean = mean
    ax.plot([mean]*2, [-0.7/sig, 1.0/sig], color['mean'])
    ax.annotate('mean', (mean, 1.0/sig), (8.0, 1.3/sig),
                ha="right", va="center", size=fontsize,
                arrowprops=dict(arrowstyle='-|>',
                                patchB=None,
                                shrinkA=5,
                                shrinkB=5,
                                fc="k", ec="k",
                                connectionstyle="arc3,rad=0.2",
                                ),
                bbox=None)

def draw_smedian(x, fig, ax):
    styles = mpatches.ArrowStyle.get_styles()

    fontsize = 30

    styleclass = styles['-|>']
    median = numpy.median(x)
    ax.plot([median]*2, [-0.7/sig, 1.0/sig], color['median'])
    ax.annotate('median', (median, 0.8/sig), (-6.0, 1.3/sig),
                ha="right", va="center", size=fontsize, color='g', 
                arrowprops=dict(arrowstyle='-|>',
                                patchB=None,
                                shrinkA=5,
                                shrinkB=5,
                                fc="g", ec="g",
                                connectionstyle="arc3,rad=0.3",
                                ),
                bbox=None)

def save_figure():
    global counter
    plt.savefig('scene%d.eps' % counter)
    counter += 1
    plt.draw()
    plt.show(block=False)
    raw_input('Press <ENTER> to continue')

def draw_idea_text(fig, ax, itr = 0, stage=1):
    if itr == 0:
        cst = ''
    else:
        cst = 'c^%d' % itr
    ax.text(-10.0, -0.3, '1. Delete Sample $t_j$ if $|t_j - $mean$| > \epsilon '+cst+'$',
        fontsize=25)
    if stage > 1:
        ax.text(-10.0, -0.4, '2. Recompute mean and Shrink threshold', fontsize=25)
    global sig, th, c
    ax.plot([current_mean - th*(c ** itr)] * 2, [-0.7/sig, 1.0/sig], 'k--')
    ax.plot([current_mean + th*(c ** itr)] * 2, [-0.7/sig, 1.0/sig], 'k--')

fig, ax = get_fig_and_ax()
draw_arrow(fig, ax)
save_figure()
draw_gt(fig, ax)
save_figure()
draw_samples(x, fig, ax)
save_figure()
draw_smean(x, fig, ax)
save_figure()
plt.close(fig)

"""Scene #2"""

sig = 3.0
shift=-4.0
x = x + shift

fig, ax = get_fig_and_ax()
draw_arrow(fig, ax, noise=True)
save_figure()

draw_samples(x, fig, ax, noise=True)
save_figure()
draw_gt(fig, ax)
all_sample = numpy.concatenate((x, skew_noise))
truncated_sample = all_sample
draw_smean(all_sample, fig, ax)
save_figure()
draw_smedian(all_sample, fig, ax)
save_figure()
draw_idea_text(fig, ax)
save_figure()
plt.close(fig)

def truncate(x, skew_noise, th, c, k):
    all_sample = numpy.concatenate((x, skew_noise))
    mean = numpy.mean(all_sample)
    deleted = []
    #truncated = []
    th_k = th * (c ** k)
    for i, y_i in enumerate(x):
        if abs(y_i - mean) > th_k:
            deleted.append(i)
            #truncated.append(y_i)
    #truncated = numpy.asarray(truncated)
    deleted = deleted[::-1]

    for idx in deleted:
        x = numpy.delete(x, obj=idx, axis=0)

    deleted = []
    for i, y_i in enumerate(skew_noise):
        if abs(y_i - mean) > th_k:
            deleted.append(i)
    deleted = deleted[::-1]

    for idx in deleted:
        skew_noise = numpy.delete(skew_noise, obj=idx, axis=0)

    return x, skew_noise

    
"""delete samples"""
fig, ax = get_fig_and_ax()
draw_arrow(fig, ax, noise=True)

draw_gt(fig, ax)
draw_smean(numpy.concatenate((x, skew_noise)), fig, ax)
draw_smedian(all_sample, fig, ax)
draw_idea_text(fig, ax, 0)
x, skew_noise = truncate(x, skew_noise, th, c, 0)
draw_samples(x, fig, ax, noise=True)
save_figure()
plt.close(fig)

"""Need to update mean and threshold"""
fig, ax = get_fig_and_ax()
draw_arrow(fig, ax, noise=True)

draw_gt(fig, ax)
draw_smean(numpy.concatenate((x, skew_noise)), fig, ax)
draw_smedian(all_sample, fig, ax)
draw_idea_text(fig, ax, 1, stage=2)
draw_samples(x, fig, ax, noise=True)
save_figure()
plt.close(fig)

for k in range(1, 5):
    x, skew_noise = truncate(x, skew_noise, th, c, k)
    fig, ax = get_fig_and_ax()
    draw_arrow(fig, ax, noise=True)

    draw_gt(fig, ax)
    draw_smean(numpy.concatenate((x, skew_noise)), fig, ax)
    draw_smedian(all_sample, fig, ax)
    draw_idea_text(fig, ax, k+1, stage=2)
    draw_samples(x, fig, ax, noise=True)
    save_figure()
    plt.close(fig)


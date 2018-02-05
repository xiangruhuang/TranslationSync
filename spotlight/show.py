import matplotlib.pyplot as plt
import numpy

def truncate(all_sample, skew_noise, th, c, k):
    mean = numpy.mean(all_sample)
    deleted = []
    truncated = []
    th_k = th * (c ** k)
    for i, y_i in enumerate(all_sample):
        if abs(y_i - mean) > th_k:
            deleted.append(i)
            truncated.append(y_i)
    truncated = numpy.asarray(truncated)
    deleted = deleted[::-1]

    for idx in deleted:
        all_sample = numpy.delete(all_sample, obj=idx, axis=0)

    deleted = []
    for i, y_i in enumerate(skew_noise):
        if abs(y_i - mean) > th_k:
            deleted.append(i)
    deleted = deleted[::-1]

    for idx in deleted:
        skew_noise = numpy.delete(skew_noise, obj=idx, axis=0)

    return all_sample, skew_noise, truncated


delta = 1.0

x = numpy.random.uniform(-delta, delta, size=20)
sym_noise = numpy.random.uniform(-5, 5, size=6)
skew_noise = numpy.random.uniform(5, 8, size=6)

color = {'gt':'c-', 'mean':'k--+', 'median':'g-+', 'trunc':'k-+'}
loc = 'upper left'
ft = 20

"""Scene #0"""

h0 = plt.figure()
plt.xlim(-10, 10)
plt.ylim(-1.0, 1.0)
#plt.tick_params(axis='x', which='both', bottom='off', top='off',
#labelbottom='off')

plt.plot([-10, -1], [0.0, 0.0], 'b')
plt.plot([1, 10], [0.0, 0.0], 'b')
plt.plot([-1, -1], [0.0, 0.5], 'b')
plt.plot([1, 1], [0.0, 0.5], 'b')
plt.plot([-1, 1], [0.5, 0.5], 'b')
#plt.plot(x, [0.0 for x_i in x], 'x', label='U[$\mu - \sigma, \mu + \sigma$]')
plt.text(2, 0.7, '$U[-\sigma, \sigma]$', fontsize=20)
plt.savefig('scene0_1.eps', bbox_inches='tight')

plt.plot([0.0]*2, [-1.0, 1.0], color['gt'])
plt.savefig('scene0_2.eps', bbox_inches='tight')

"""Scene #1"""

h1 = plt.figure()
plt.xlim(-10, 10)
plt.ylim(-1.0, 1.0)
plt.axis('off')

plt.plot([0.0]*2, [-1.0, 1.0], color['gt'])
plt.plot(x, [0.0 for x_i in x], 'x', label='U[$\mu - \sigma, \mu + \sigma$]')
#plt.legend(loc=loc, fontsize=ft)
plt.savefig('scene1_1.eps', bbox_inches='tight')

plt.plot([numpy.mean(x)]*2, [-0.8, 0.8], color['mean'], label='mean')
#plt.legend(loc=loc, fontsize=ft)
plt.savefig('scene1_2.eps', bbox_inches='tight')

"""Scene #2: Adding noise in"""

h2 = plt.figure()
plt.xlim(-10, 10)
plt.ylim(-1.0, 1.0)
plt.axis('off')

plt.plot([0.0]*2, [-1.0, 1.0], color['gt'])
plt.plot(x, [0.0 for x_i in x], 'x', label='U[$\mu - \sigma, \mu + \sigma$]')

all_sample = numpy.concatenate((x, sym_noise))

plt.plot([numpy.mean(all_sample)]*2, [-0.8, 0.8], color['mean'], label='mean')
plt.plot(sym_noise, [0.0 for n_i in sym_noise], 'rx', label='noise')
#plt.legend(loc=loc, fontsize=ft)
plt.savefig('scene2_1.eps', bbox_inches='tight')


"""Scene #3: Skewed Noise"""

h3 = plt.figure()
plt.xlim(-10, 10)
plt.ylim(-1.0, 1.0)
plt.axis('off')

plt.plot([0.0]*2, [-1.0, 1.0], color['gt'])
plt.plot(x, [0.0 for x_i in x], 'x', label='U[$\mu - \sigma, \mu + \sigma$]')

all_sample = numpy.concatenate((x, skew_noise))

plt.plot([numpy.mean(all_sample)]*2, [-0.8, 0.8], color['mean'], label='mean')
plt.plot(skew_noise, [0.0 for n_i in skew_noise], 'rx', label='noise')
#plt.legend(loc=loc, fontsize=ft)
plt.savefig('scene3_1.eps', bbox_inches='tight')

plt.plot([numpy.median(all_sample)]*2, [-0.8, 0.8], color['median'], label='median')
#plt.legend(loc=loc, fontsize=ft)
plt.savefig('scene3_2.eps', bbox_inches='tight')


"""Scene #4-9: Remove Noise 1-6"""

trunc_sample = numpy.copy(all_sample)

th = 6
c = 0.9
    
h = plt.figure()
plt.xlim(-10, 10)
plt.ylim(-1.0, 1.0)
plt.axis('off')
plt.text(2, 0.7, 'Iteratively Delete', fontsize=20)
plt.text(2, 0.5, '$\{j | \quad |x^{(k)} - t_j| < \epsilon c^k\}$', fontsize=20)

plt.plot([0.0]*2, [-1.0, 1.0], color['gt'])
plt.plot(x, [0.0 for x_i in x], 'x', label='U[$\mu - \sigma, \mu + \sigma$]')

#trunc_sample, skew_noise, deleted = truncate(trunc_sample, skew_noise, th, c, scene)
#for d in deleted:
#    plt.arrow(d, 0.7, 0, -0.6, shape='full', lw=0, length_includes_head=True,
#            head_width=.15)
plt.plot([numpy.mean(all_sample)]*2, [-0.8, 0.8], color['mean'], label='mean')
#plt.plot([numpy.mean(trunc_sample)]*2, [-0.8, 0.8], color['trunc'], label='truncated')
plt.plot(skew_noise, [0.0 for n_i in skew_noise], 'rx', label='noise')
plt.plot([numpy.median(all_sample)]*2, [-0.8, 0.8], color['median'], label='median')
#plt.scatter(deleted, [0.0 for d in deleted], s=100, edgecolors='b',
#        facecolors='none')
#plt.legend(loc=loc, fontsize=ft)
plt.savefig('scene%d_1.eps' % (4), bbox_inches='tight')

for scene in range(6):
    h = plt.figure()
    plt.xlim(-10, 10)
    plt.ylim(-1.0, 1.0)
    plt.axis('off')
    plt.text(2, 0.7, 'Iteratively Delete', fontsize=20)
    plt.text(2, 0.5, '$\{j | \quad |x^{(k)} - t_j| < \epsilon c^k\}$', fontsize=20)

    plt.plot([0.0]*2, [-1.0, 1.0], color['gt'])
    plt.plot(x, [0.0 for x_i in x], 'x', label='U[$\mu - \sigma, \mu + \sigma$]')

    trunc_sample, skew_noise, deleted = truncate(trunc_sample, skew_noise, th, c, scene)


    #for d in deleted:
    #    plt.arrow(d, 0.7, 0, -0.6, shape='full', lw=0, length_includes_head=True,
    #            head_width=.15)
    plt.plot([numpy.mean(all_sample)]*2, [-0.8, 0.8], color['mean'], label='mean')
    plt.plot([numpy.mean(trunc_sample)]*2, [-0.8, 0.8], color['trunc'], label='truncated')
    plt.plot(skew_noise, [0.0 for n_i in skew_noise], 'rx', label='noise')
    plt.plot([numpy.median(all_sample)]*2, [-0.8, 0.8], color['median'], label='median')
    #plt.scatter(deleted, [0.0 for d in deleted], s=100, edgecolors='b',
    #        facecolors='none')
    #plt.legend(loc=loc, fontsize=ft)
    plt.savefig('scene%d_1.eps' % (scene+5), bbox_inches='tight')
    if len(skew_noise) == 0:
        break


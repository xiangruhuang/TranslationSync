import numpy
import os

def process(folder):
    min_TL2 = []
    median_TL2 = []
    max_TL2 = []

    min_CD = []
    median_CD = []
    max_CD = []

    tmean_TL2 = []
    tmean_CD = []

    zp_TL2 = []
    zp_CD = []
    ratios = range(100)
    for ratio in ratios:
        ml_TL2 = []
        time_TL2 = []
        ml_CD = []
        time_CD = []
        with open(folder+'/ratio'+str(ratio)+'_summary', 'r') as fin:
            lines = fin.readlines()
            assert len(lines) == 100
            for line in lines:
                vals = [float(token) for token in line.strip().split(' ')]
                ml_TL2.append(vals[0])
                time_TL2.append(vals[1])
                ml_CD.append(vals[2])
                time_CD.append(vals[3]) 
        if len(ml_TL2) > 0 and len(ml_CD) > 0:
            min_TL2.append(min(ml_TL2))
            median_TL2.append(numpy.median(ml_TL2))
            max_TL2.append(max(ml_TL2))
            min_CD.append(min(ml_CD))
            median_CD.append(numpy.median(ml_CD))
            max_CD.append(max(ml_CD))
            zp_TL2.append(len([e for e in ml_TL2 if abs(e) < 1e-2])/100.0)
            zp_CD.append(len([e for e in ml_CD if abs(e) < 1e-2])/100.0)
            tmean_CD.append(numpy.mean(time_CD))
            tmean_TL2.append(numpy.mean(time_TL2))

    ratios = [1.0-ratio/100.0 for ratio in ratios]
    return min_TL2, median_TL2, max_TL2, min_CD, median_CD, max_CD, tmean_TL2, \
tmean_CD, zp_TL2, zp_CD, ratios

def process_in_details(folder):
    min_TL2 = []
    median_TL2 = []
    max_TL2 = []

    min_CD = []
    median_CD = []
    max_CD = []

    tmean_TL2 = []
    tmean_CD = []

    zp_TL2 = []
    zp_CD = []
    ratios = range(100)
    for ratio in ratios:
        ml_TL2 = []
        time_TL2 = []
        ml_CD = []
        time_CD = []
        for eid in range(100):
            name_TL2 = folder+'/ratio'+str(ratio)+'_'+str(eid)+'.TL2'
            name_CD = folder+'/ratio'+str(ratio)+'_'+str(eid)+'.CD'
            if os.path.isfile(name_TL2) and os.path.isfile(name_CD):
                with open(name_TL2) as fin:
                    lines = fin.readlines()
                    min_loss=1e100
                    time=0.0
                    for line in lines:
                        ml = 0.0
                        t = 0.0
                        for token in line.strip().split(', '):
                            if token.startswith('min_loss'):
                                ml = float(token.split('=')[1])
                            if token.startswith('elapsed'):
                                t = float(token.split('=')[1])
                        if min_loss > ml:
                            min_loss = ml
                            time = t
                        if min_loss < 1e-5:
                            break
                    ml_TL2.append(min_loss)
                    time_TL2.append(time)

                with open(name_CD) as fin:
                    lines = fin.readlines()
                    min_loss=1e100
                    time=0.0
                    for line in lines:
                        ml = 0.0
                        t = 0.0
                        for token in line.strip().split(', '):
                            if token.startswith('min_loss'):
                                ml = float(token.split('=')[1])
                            if token.startswith('elapsed'):
                                t = float(token.split('=')[1])
                        if min_loss > ml:
                            min_loss = ml
                            time = t
                        if min_loss < 1e-5:
                            break
                    ml_CD.append(min_loss)
                    time_CD.append(time)
        if len(ml_TL2) > 0 and len(ml_CD) > 0:
            min_TL2.append(min(ml_TL2))
            median_TL2.append(numpy.median(ml_TL2))
            max_TL2.append(max(ml_TL2))
            min_CD.append(min(ml_CD))
            median_CD.append(numpy.median(ml_CD))
            max_CD.append(max(ml_CD)) 
            zp_TL2.append(len([e for e in ml_TL2 if abs(e) < 1e-2])/100.0)
            zp_CD.append(len([e for e in ml_CD if abs(e) < 1e-2])/100.0)
            tmean_CD.append(numpy.mean(time_CD))
            tmean_TL2.append(numpy.mean(time_TL2))
    
    ratios = [1.0-ratio/100.0 for ratio in ratios]
    return min_TL2, median_TL2, max_TL2, min_CD, median_CD, max_CD, tmean_TL2, \
tmean_CD, zp_TL2, zp_CD, ratios

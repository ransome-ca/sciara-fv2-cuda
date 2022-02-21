#!/bin/env python3

from cProfile import label
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import re


BW_HBM =  224.32
BW_L1  = 2119.98
BW_L2  =  755.02

PF_FP32 = 4978.0
PF_FP64 =  155.6
PE_FP64 =  153.2



parser = argparse.ArgumentParser(description='Plot roofline')
parser.add_argument('-i', '--input', type=str, required=True, help='Input directory')
parser.add_argument('-o', '--output', type=str, required=True, help='Output image')
parser.add_argument('-s', '--show', action='store_true', help='Show plot')
parser.add_argument('-t', '--type', type=str, required=False, default='cuda-multi-gpu', help='Type of cuda version')
args = parser.parse_args()

if not os.path.isdir(args.input):
    raise RuntimeError('Input directory does not exist')



db = {}


print('Loading data')

for version in os.listdir(args.input):

    print(' + Processing version {}'.format(version))

    db[version] = {}
    
    for runtime in os.listdir(os.path.join(args.input, version)):

        if not runtime.endswith('8_runtime.csv'):
            continue


        print('   - Processing runtime {}'.format(runtime))

        with open(os.path.join(args.input, version, runtime), 'r') as fp:
            
            csv = fp.readlines()
            csv = csv[3:]
            
            for row in csv:
                    
                if not row.startswith('"GPU activities"'):
                    continue

                row = re.findall(r'".+?"|[\w\.]+', row)

                if row[7].find('(') == -1:
                    continue

                k = row[7]
                k = k.split('(')
                k = k[0]
                k = k.replace('"', '')
                k = k.strip()

                db[version][k] = {}
                db[version][k]['time'] = {}
                db[version][k]['time']['min'] = float(row[5]) / 1000.0
                db[version][k]['time']['max'] = float(row[6]) / 1000.0
                db[version][k]['time']['avg'] = float(row[4]) / 1000.0


    for k in os.listdir(os.path.join(args.input, version)):

        if not os.path.isdir(os.path.join(args.input, version, k)):
            continue


        print('   - Processing kernel {}'.format(k))

        with open(os.path.join(args.input, version, k, '%s_%s_metrics.csv' % (version, k)), 'r') as fp:
            
            csv = fp.readlines()
            csv = csv[5:]


            db[version][k]['metrics'] = {}

            for row in csv:

                if not row.startswith('"NVIDIA'):
                    continue

                row = re.findall(r'".+?"|[\w\.]+', row)

                row[3] = row[3].replace('"', '')
                row[3] = row[3].strip()

                if row[3] in db[version][k]['metrics']:
                    db[version][k]['metrics'][row[3]] = float(row[7]) + db[version][k]['metrics'][row[3]]
                else:
                    db[version][k]['metrics'][row[3]] = float(row[7])



print('Processing data')

for version in db:

    print(' + Processing version {}'.format(version))

    for k in db[version]:

        if not 'metrics' in db[version][k]:
            continue

        if not 'time' in db[version][k]:
            continue

        if not 'avg' in db[version][k]['time']:
            continue


        print('   - Processing kernel {}'.format(k))



        flp = 0.0
        flp += db[version][k]['metrics']['flop_count_dp'] if 'flop_count_dp' in db[version][k]['metrics'] else 0.0
        flp += db[version][k]['metrics']['flop_count_sp'] if 'flop_count_sp' in db[version][k]['metrics'] else 0.0
        flp += db[version][k]['metrics']['flop_count_hp'] if 'flop_count_hp' in db[version][k]['metrics'] else 0.0

        db[version][k]['perf'] = float(flp / db[version][k]['time']['avg']) if db[version][k]['time']['avg'] > 0.0 else 0.0

        print('     > Performance: {} GFLOP/s'.format(round(db[version][k]['perf'] / 1e9, 2)))



        db[version][k]['ai'] = {}

        ops = 0.0
        ops += db[version][k]['metrics']['dram_read_transactions'] if 'dram_read_transactions' in db[version][k]['metrics'] else 0.0
        ops += db[version][k]['metrics']['dram_write_transactions'] if 'dram_read_transactions' in db[version][k]['metrics'] else 0.0

        db[version][k]['ai']['HBM'] = float(flp / (ops * 32)) if ops > 0.0 else 0.0

        print('     > HBM: {} FLOP/byte'.format(round(db[version][k]['ai']['HBM'], 3)))


        ops = 0.0
        ops += db[version][k]['metrics']['l2_read_transactions'] if 'l2_read_transactions' in db[version][k]['metrics'] else 0.0
        ops += db[version][k]['metrics']['l2_write_transactions'] if 'l2_write_transactions' in db[version][k]['metrics'] else 0.0

        db[version][k]['ai']['L2'] = float(flp / (ops * 32)) if ops > 0.0 else 0.0

        print('     > L2: {} FLOP/byte'.format(round(db[version][k]['ai']['L2'], 3)))


        ops = 0.0
        ops += db[version][k]['metrics']['gld_transactions'] if 'gld_transactions' in db[version][k]['metrics'] else 0.0
        ops += db[version][k]['metrics']['gst_transactions'] if 'gst_transactions' in db[version][k]['metrics'] else 0.0
        ops += db[version][k]['metrics']['atomic_transactions'] if 'atomic_transactions' in db[version][k]['metrics'] else 0.0
        ops += db[version][k]['metrics']['local_load_transactions'] if 'local_load_transactions' in db[version][k]['metrics'] else 0.0
        ops += db[version][k]['metrics']['local_store_transactions'] if 'local_store_transactions' in db[version][k]['metrics'] else 0.0
        ops += db[version][k]['metrics']['shared_load_transactions'] if 'shared_load_transactions' in db[version][k]['metrics'] else 0.0
        ops += db[version][k]['metrics']['shared_store_transactions'] if 'shared_store_transactions' in db[version][k]['metrics'] else 0.0

        db[version][k]['ai']['L1'] = float(flp / (ops * 32)) if ops > 0.0 else 0.0

        print('     > L1: {} FLOP/byte'.format(round(db[version][k]['ai']['L1'], 3)))



print('Plotting data')


colors = ['#0b84a5','#9dd866','#ca472f','y','m','c']
styles = ['o','D','s','o','D',">","<","*","h","H","+","1","2","3","4","8","p","d","|","_",".",","]



fig = plt.figure(figsize=(10.67, 6.6))
plt.clf()


xmin = 2.5e-2
xmax = 1e3
ymin = 1e1
ymax = 1e4

ax = fig.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Arithimetic intensity (FLOP/byte)')
ax.set_ylabel('Performance (GFLOP/s)')
ax.set_title(str(args.type).replace('-', ' ').upper())
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)


if not args.type in db:
    raise Exception('Invalid version')


x = np.linspace(xmin, xmax, 1000000)

x_fp32_hbm = x[np.argwhere(x * BW_HBM < PF_FP32).flatten()]
x_fp32_l1  = x[np.argwhere(x * BW_L1 < PF_FP32).flatten()]
x_fp32_l2  = x[np.argwhere(x * BW_L2 < PF_FP32).flatten()]

x_fp64_hbm = x[np.argwhere(x * BW_HBM < PF_FP64).flatten()]
x_fp64_l1  = x[np.argwhere(x * BW_L1 < PF_FP64).flatten()]
x_fp64_l2  = x[np.argwhere(x * BW_L2 < PF_FP64).flatten()]


ax.plot(x_fp32_hbm, x_fp32_hbm * BW_HBM, color='gray', linestyle='--')
ax.plot(x_fp32_l1, x_fp32_l1 * BW_L1, color='gray', linestyle='--')
ax.plot(x_fp32_l2, x_fp32_l2 * BW_L2, color='gray', linestyle='--')

ax.plot(x_fp64_hbm, x_fp64_hbm * BW_HBM, color='black', linestyle='-')
ax.plot(x_fp64_l1, x_fp64_l1 * BW_L1, color='black', linestyle='-')
ax.plot(x_fp64_l2, x_fp64_l2 * BW_L2, color='black', linestyle='-')

ax.text(6.2, ymax - (ymax * 0.95), 'HBM: %.2f GB/s' % BW_HBM, verticalalignment='bottom', horizontalalignment='right', color='gray', rotation=45)
ax.text(0.6, ymax - (ymax * 0.95), 'L1: %.2f GB/s' % BW_L1, verticalalignment='bottom', horizontalalignment='right', color='gray', rotation=45)
ax.text(1.6, ymax - (ymax * 0.95), 'L2: %.2f GB/s' % BW_L2, verticalalignment='bottom', horizontalalignment='right', color='gray', rotation=45)




bound_hbm = PF_FP64 / BW_HBM
bound_l1  = PF_FP64 / BW_L1
bound_l2  = PF_FP64 / BW_L2

ax.vlines(x=bound_hbm, ymin=ymin, ymax=PF_FP64, linestyle='--', color=colors[0], alpha=0.2)
ax.vlines(x=bound_l1, ymin=ymin, ymax=PF_FP64, linestyle='--', color=colors[1], alpha=0.2)
ax.vlines(x=bound_l2, ymin=ymin, ymax=PF_FP64, linestyle='--', color=colors[2], alpha=0.2)



ax.hlines(y=PF_FP32, xmin=x_fp32_l1[-1] if len(x_fp32_l1) else 0.0, xmax=xmax, color='gray', linestyle='--')
ax.text(xmax - (xmax * 0.1), PF_FP32 + (PF_FP32 * 0.05), 'FP32: %.2f GFLOP/s' % PF_FP32, verticalalignment='bottom', horizontalalignment='right', color='gray')

ax.hlines(y=PF_FP64, xmin=x_fp64_l1[-1] if len(x_fp64_l1) else 0.0, xmax=xmax, color='black', linestyle='-')
ax.text(xmax - (xmax * 0.1), PF_FP64 + (PF_FP64 * 0.05), 'FP64: %.2f GFLOP/s' % PF_FP64, verticalalignment='bottom', horizontalalignment='right')






patch_handles = list()
marker_handles = list()

for k in db[args.type]:

    if not 'ai' in db[args.type][k]:
        continue

    if not 'perf' in db[args.type][k]:
        continue

    if(db[args.type][k]['perf'] / 1e9 < ymin):
        continue

    ax.plot(db[args.type][k]['ai']['HBM'], db[args.type][k]['perf'] / 1e9, color=colors[0], marker=styles[list(db[args.type].keys()).index(k)])
    ax.plot(db[args.type][k]['ai']['L1'],  db[args.type][k]['perf'] / 1e9, color=colors[1], marker=styles[list(db[args.type].keys()).index(k)])
    ax.plot(db[args.type][k]['ai']['L2'],  db[args.type][k]['perf'] / 1e9, color=colors[2], marker=styles[list(db[args.type].keys()).index(k)])

    marker_handles.append(ax.plot([], [], color='gray', marker=styles[list(db[args.type].keys()).index(k)], label=k)[0])



patch_handles.append(mpatches.Patch(color=colors[1], label='L1'))
patch_handles.append(mpatches.Patch(color=colors[2], label='L2'))
patch_handles.append(mpatches.Patch(color=colors[0], label='HBM'))


leg1 = plt.legend(handles=marker_handles, loc=4, ncol=1)
leg2 = plt.legend(handles=patch_handles, loc='lower right', bbox_to_anchor=(0.545, 0), scatterpoints=1)

ax.add_artist(leg1)


if args.show:
    plt.show()

plt.savefig(args.output)
import numpy as np
import os
import re 
import subprocess
import sys
import multiprocessing


def parse_output(out):
    total_paths = []
    for m in re.finditer(r'Total_paths: (\d+[.]?\d+[E]?\d+)', out):
        total_paths.append(int(eval(m.group(1))))
    core_sizes = []
    for m in re.finditer(r'Core_size: (\d+)', out):
        core_sizes.append(int(m.group(1)))
    m = re.search(r'H_score: ([-]?0.\d+)', out)
    h_score = float(m.group(1))
    return total_paths, core_sizes, h_score

NUM_EXPERIMENTS = 20
#          S=M=T,     S=2M=2T,    T=2S=2M,   M=2S=2T,  S=T=2M
# ratios = [(1, 1, 1), (2, 1, 1), (1, 1, 2), (1, 2, 1), (2, 2, 1)]
# degrees = [(1, 2, 3), 1, 1, 1, 1]
ratios = [(1, 1, 1)]
degrees = [(1, 2, 3)]
alphas = [a / 10.0 for a in xrange(-20, 22, 2)]

for ratio, degree in zip(ratios, degrees):
    for d in [degree] if isinstance(degree, (int, long)) else degree:
        # for v in xrange(200, 400, 200):
        for v in xrange(200, 2000, 200):
            s, m, t = [int(round(float(v * r)/sum(ratio))) for r in ratio]
            # for a in [-2.0]:
            for a in alphas:
                pool = multiprocessing.Pool(4)
                arguments = ['python', 'main.py', '-v', str(v), '-s', str(s), '-t', str(t), '-d', str(d), '-a', str(a), '-n']
                results = pool.map(subprocess.check_output, (arguments + [str(n)] for n in xrange(NUM_EXPERIMENTS)))
                for n in xrange(NUM_EXPERIMENTS):
                    results_mine = parse_output(results[n])
                    print '%d\t%d\t%d\t%d\t%d\t%f\t%d\t'%(v, s, t, (v - s - t), d, a, n),
                    print '%d\t%d\t%f'%(results_mine[1][0], results_mine[1][1], results_mine[2])
                print
                sys.stdout.flush()
                # for n in xrange(NUM_EXPERIMENTS):
                    # arguments = ['python', 'main.py', '-v', str(v), '-s', str(s), '-t', str(t), '-d', str(d), '-a', str(a), '-n', str(n)]
                    # arguments = ['python', 'main.py', '-v', str(v), '-s', str(s), '-t', str(t), '-d', str(d), '-a', str(a), '-n', str(n), '-o', 'Original/paper']
                    # print ' '.join(arguments)
                    # out_mine = subprocess.check_output(arguments)
                    # results_mine = parse_output(out_mine)

                    # os.chdir('Original')
                    # lines = subprocess.check_output(['java', 'HourglassAnalysis']).split('\n')
                    # os.chdir('..')

                    # lines = filter(lambda s: s.startswith('Original_network') or s.startswith('Flat_network') or s.startswith('Total_paths') or s.startswith('Core_size') or s.startswith('H_score'), lines)
                    # for i in xrange(len(lines)):
                        # if lines[i].startswith('Core_size'):
                            # lines[i] = lines[i][0: lines[i].find('Core_set')-1]
                    # out_theirs = '\n'.join(lines)
                    # results_theirs = parse_output(out_theirs)

                    # if (a <= 1.2) and (results_mine[0][0] == results_mine[0][1]) and (results_theirs[0][0] == results_theirs[0][1]):
                        # assert results_mine[0][0] == results_theirs[0][0], '%d != %d'%(results_mine[0][0], results_theirs[0][0])

                    # print '%d\t%d\t%d\t%d\t%d\t%f\t%d\t'%(v, s, t, (v - s - t), d, a, n),
                    # print '%d\t%d\t%f'%(results_mine[1][0], results_mine[1][1], results_mine[2])

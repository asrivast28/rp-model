import sys
# import multiprocessing

import hourglass


def parse_output(out):
    import re
    total_paths = []
    for m in re.finditer(r'Total_paths: (\d+[.]?\d+[E]?\d+)', out):
        total_paths.append(int(eval(m.group(1))))
    P, P_f = total_paths
    core_sizes = []
    for m in re.finditer(r'Core_size: (\d+)', out):
        core_sizes.append(int(m.group(1)))
    C, C_f = core_sizes
    m = re.search(r'H_score: ([-]?0.\d+)', out)
    H = float(m.group(1))
    return P, C, P_f, C_f, H

NUM_EXPERIMENTS = 20
#          S=M=T,     S=M=T,     S=M=T,     S=2M=2T,    T=2S=2M,   M=2S=2T,  S=T=2M
ratios = [(1, 1, 1), (1, 1, 1), (1, 1, 1), (2, 1, 1), (1, 1, 2), (1, 2, 1), (2, 2, 1)]
degrees = [1, 2, 3, 1, 1, 1, 1]
if len(sys.argv) > 1:
    index = int(sys.argv[1])
    ratios = [ratios[index]]
    degrees = [degrees[index]]
# ratios = [(1, 1, 1)]
# degrees = [(1, 2, 3)]
alphas = [a / 10.0 for a in xrange(-20, 22, 2)]

# pool = multiprocessing.Pool(1)
for ratio, degree in zip(ratios, degrees):
    # for v in xrange(200, 400, 200):
    for v in xrange(200, 2200, 200):
        s, m, t = [int(round(float(v * r)/sum(ratio))) for r in ratio]
        t = (v - s - m)
        # for a in [-2.0]:
        for a in alphas:
            arguments = ['-v', str(v), '-s', str(s), '-t', str(t), '-d', str(degree), '-a', str(a), '-n']
            # results = pool.map(hourglass.properties, (arguments + [str(n)] for n in xrange(NUM_EXPERIMENTS)))
            for n in xrange(NUM_EXPERIMENTS):
                # P, C, P_f, C_f, H = results[n]
                P, C, P_f, C_f, H = hourglass.properties(arguments + [str(n)])
                print '%d\t%d\t%d\t%d\t%d\t%f\t%d\t'%(v, s, t, (v - s - t), degree, a, n),
                print '%d\t%d\t%f'%(C, C_f, H)
            print
            sys.stdout.flush()
            # for n in xrange(NUM_EXPERIMENTS):
                # arguments = ['-v', str(v), '-s', str(s), '-t', str(t), '-d', str(d), '-a', str(a), '-n', str(n), '-o', 'Original/paper']
                # result_mine = hourglass.properties(arguments)

                # import os
                # import subprocess
                # os.chdir('Original')
                # lines = subprocess.check_output(['java', 'HourglassAnalysis']).split('\n')
                # os.chdir('..')

                # lines = filter(lambda s: s.startswith('Original_network') or s.startswith('Flat_network') or s.startswith('Total_paths') or s.startswith('Core_size') or s.startswith('H_score'), lines)
                # for i in xrange(len(lines)):
                    # if lines[i].startswith('Core_size'):
                        # lines[i] = lines[i][0: lines[i].find('Core_set')-1]
                # out_theirs = '\n'.join(lines)
                # result_theirs = parse_output(out_theirs)

                # print '%d\t%d\t%d\t%d\t%d\t%f\t%d\t'%(v, s, t, (v - s - t), d, a, n),
                # print '%d\t%d\t%d\t%d\t%f\t%f'%(result_mine[1], result_theirs[1], result_mine[3], result_theirs[3], result_mine[4], result_theirs[4])

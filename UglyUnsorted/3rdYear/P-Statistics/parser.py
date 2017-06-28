import sys


def process_file(filename):
    fl = (open(filename, 'r'))
    f = fl.read().split('\n-')
    for block in f:
        if (block == ''):
            continue

        lines = block.split('\n')

        sum_cancer = 0.0
        count_cancer = 0
        sum_control = 0.0
        count_control = 0

        for l in lines:

            if 'cancer' in l:
                d = [int(s) for s in l.split() if s.isdigit()]
                sum_cancer += float(d[0])/float(d[1])
                count_cancer += 1
            elif 'control' in l:
                d = [int(s) for s in l.split() if s.isdigit()]
                sum_control += float(d[0])/float(d[1])
                count_control += 1

        print(sum_cancer/float(count_cancer),
              sum_control/float(count_control))

for arg in sys.argv:
    if not arg == __file__:
        process_file(arg)

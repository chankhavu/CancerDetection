#!/bin/python

import os
import re
import random
import subprocess


def collect_data(source_folder, regex_match, destination_folder,
                 destination, suffix, ratio):
    dest_train = open(os.path.join(destination_folder,
                                   destination+'_'+suffix), 'w')
    dest_test = open(os.path.join(destination_folder,
                                  'test_'+destination+'_'+suffix), 'w')

    pattern = re.compile(regex_match)
    files = [file for file in os.listdir(source_folder) if
             pattern.match(file)]
    random.shuffle(files)
    counter = 0

    print ('  training items', '('+destination+'):',
           int(float(len(files))*float(ratio)))
    print ('  testing itmes', '('+destination+'):',
           int(float(len(files)) * (float(1)-float(ratio))))

    for file in files:
        counter += 1
        f = open(os.path.join(source_folder, file), 'r')
        values = [x for x in f.read().split()]
        f.close()

        if float(counter) < float(len(files)) * float(ratio):
            dest_train.write(' '.join(values) + '\n')
        else:
            dest_test.write(' '.join(values) + '\n')

    dest_train.close()
    dest_test.close()


def prepare_data(dirs, ratio, color):
    for d in dirs:
        sf = d[0]
        df = d[1]
        name = d[2]

        collect_data(sf, '^[0-9]+_red$', os.path.join(df, color), name,
                     ratio[1], ratio[0])


def cross_validate(dirs, ratios, color, tmp_dir, K):

    result = open('%s_summary2_%d' % (color, K), 'a+')

    for r in ratios:
        print ('processing', color, 'filter data with koefficient', K,
               'and ratio', r[1])

        result.write('\n----- Ratio = %s, K = %d -----\n' % (r[1], K))

        for repeat in range(10):
            prepare_data(dirs, r, color)

            for d in dirs:
                args = [
                    os.path.join('bin', 'knn'),
                    str(K),
                ]

                for _d in dirs:
                    args = args + [os.path.join(tmp_dir, color, _d[2]+'_'+r[1]), ]

                args = args + [os.path.join(tmp_dir, color,
                                            'test_'+d[2]+'_'+r[1]), ]

                output = ((subprocess.check_output(args)).decode('ascii')).split()
                right = int(output.count(os.path.join(tmp_dir, color,
                                                      d[2]+'_'+r[1])))

                count = int(len(output))
                result.write('%d.  %s result: %d out of %d\n' % (repeat+1, d[2], right, count))

    result.close()


dirs = [
    (
        '/home/falcon/Workspace/Klyushin/analysis/FractalDimensionAnalysis/' +
        'patients/cancer/short/',
        '/home/falcon/Workspace/Klyushin/analysis/P-Statistics/patients/',
        'cancer',
    ), (
        '/home/falcon/Workspace/Klyushin/analysis/FractalDimensionAnalysis/' +
        'patients/control/short/',
        '/home/falcon/Workspace/Klyushin/analysis/P-Statistics/patients/',
        'control',
    ),
]

ratios = [
    (5/10, '5-10'),
    (6/10, '6-10'),
    (7/10, '7-10'),
    (8/10, '8-10'),
    (9/10, '9-10'),
]


cross_validate(dirs, ratios, 'red', 'patients', 6)
cross_validate(dirs, ratios, 'green', 'patients', 6)
cross_validate(dirs, ratios, 'blue', 'patients', 6)

cross_validate(dirs, ratios, 'red', 'patients', 7)
cross_validate(dirs, ratios, 'green', 'patients', 7)
cross_validate(dirs, ratios, 'blue', 'patients', 7)

cross_validate(dirs, ratios, 'red', 'patients', 8)
cross_validate(dirs, ratios, 'green', 'patients', 8)
cross_validate(dirs, ratios, 'blue', 'patients', 8)

cross_validate(dirs, ratios, 'red', 'patients', 9)
cross_validate(dirs, ratios, 'green', 'patients', 9)
cross_validate(dirs, ratios, 'blue', 'patients', 9)

cross_validate(dirs, ratios, 'red', 'patients', 10)
cross_validate(dirs, ratios, 'green', 'patients', 10)
cross_validate(dirs, ratios, 'blue', 'patients', 10)

cross_validate(dirs, ratios, 'red', 'patients', 12)
cross_validate(dirs, ratios, 'green', 'patients', 12)
cross_validate(dirs, ratios, 'blue', 'patients', 12)

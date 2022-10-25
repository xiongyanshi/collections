'''
Read qstat output from linux pipe, tail it with executing working dir and script name.

By: xiongyanshi
20200911 first implemant.

Contact: 15201921363, xiongyanshi@genomics.cn
Usage  : qstat | qx

History:
20210812: get qstat output from linux pipe.

'''

import sys
import re
import subprocess as sps
from subprocess import PIPE


def get_qstat():
        ''' from linux pipe or run qstat by itself'''
        if sys.stdin.isatty():
            qstat = sps.run('qstat', stdout=PIPE, shell=True, check=True, encoding='utf-8')
            return qstat.stdout
        qstat = sys.stdin.read()
        return qstat

def find_cwd_and_script(jobid):
    cwd     = 'NA'   # qsub -cwd
    script  = 'NA'   # bash.sh
    qstat_j = ''

    try:
        qstat_j = sps.run(['qstat','-j', jobid], stdout=PIPE, stderr=PIPE, check=True, encoding='utf-8').stdout
    except:
        pass

    for line in qstat_j.split('\n'):
        if line.startswith('cwd:'):
            cwd = line.strip().split()[-1]
            continue
        if line.startswith('script_file:'):
            script_file = line.strip().split()[-1]
            script = script_file.split('/')[-1]

    return cwd, script

def main():
    qstat_out = get_qstat()
    if '-g' in sys.argv:
        d = {}  # {cwd: [jobid, jobid, ...]}
        for line in qstat_out.strip('\n').split('\n'):
            match = re.match(r'^\d+', line.lstrip())
            if match:
                jobid = match.group()
                cwd, script = find_cwd_and_script(jobid)
                if cwd not in d:
                    d[cwd] = []
                else:
                    d[cwd].append(jobid)

        for cwd in sorted(d.keys()):
            l_jobid = d[cwd]
            if len(l_jobid) == 0:
                continue
            print('{}:\t{}'.format(cwd, ' '.join(l_jobid)))
    else:
        for line in qstat_out.strip('\n').split('\n'):
            ll = line.split()
            if line.startswith('job-ID'):
                print(' '.join(ll[:7] + ['cwd','script']))
                continue
            if line.startswith('-------------'):
                print(line)
                continue
            match = re.match(r'^\d+', line.lstrip())
            if match:
                jobid = match.group()
                cwd, script = find_cwd_and_script(jobid)
                if len(cwd) > 52:
                    cwd = '...' + cwd[-52:]
                else:
                    cwd = '{:>55s}'.format(cwd)
                print(' '.join(ll[:7] + [cwd, script]))
            else:
                print()

if __name__ == "__main__":
    main()


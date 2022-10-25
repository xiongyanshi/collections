#!/usr/bin/env python3

# python3 compatible only.

import os
import sys
import re
import shutil
import argparse
import subprocess as sps


def find_mark_index(n_line, b_line, target):
    '''
    find the index to place a ^ mark.
    input:
        number line, line 1 of samtools tview output
        reference bases line, line 2
        interested position
    '''

    match = re.search(r'\d+', n_line)
    n1 = int(match.group(0))        # number start the line.
    n1_index = match.start(0)

    n = n1
    index = n1_index
    for b in b_line[n1_index+1:]:
        index += 1
        if b in "acgtnACGTN":       # skip '*' in reference line.
            n += 1
        if n == target:
            break

    return index


def line_core(line, index):
    '''
    for each line, only keep reads which overlap with target index position.
    '''
    l = index
    r = index
    while l > 0:
        l -= 1
        if line[l].isspace():
            break
    while r < len(line):
        if line[r].isspace():
            break
        r += 1

    return ' '*l + line[l:r] + ' '*(len(line)-r)


def get_tview(chrm, pos, bam):

    cmd = "export COLUMNS=181; \
          {samtools} tview -d T -p {chrm}:{start} {bam} {ref}".format(
           samtools=SAMTOOLS, chrm=chrm, start=pos-90, bam=bam, ref=GENOME)
    print(cmd)
    run = sps.run(cmd, shell=True, stdout=sps.PIPE)
    run.check_returncode()
    tv_lines = run.stdout.decode('utf-8').strip('\n').split('\n')

    return tv_lines


def snv(chrm, pos, bam, view_all):

    tv_lines = get_tview(chrm, pos, bam)
    mark_index = find_mark_index(tv_lines[0], tv_lines[1], pos)

    print(tv_lines[0])
    print(tv_lines[1])
    print(tv_lines[2])
    print(' '*mark_index + '^')
    base_ref = tv_lines[1][mark_index]
    base_n = {'A':0, 'C':0, 'G':0, 'T':0,
              'a':0, 'c':0, 'g':0, 't':0,
              'N':0, '*':0, 'n':0}
    line_print = []
    for line in tv_lines[3:]:
        base_query = line[mark_index]

        if base_query == ' ':
            continue

        if base_query == '.':
            base_real = base_ref.upper()
        elif base_query == ',':
            base_real = base_ref.lower()
        else:
            base_real = base_query

        base_n[base_real] += 1

        if (not view_all) and (base_query in '.,' or base_query == base_ref):
            continue

        line = line_core(line, mark_index)
        line_print.append(line)

    print('\n'.join(sorted(line_print, reverse=True)))

    template = '      {}     {}     {}     {}     {}\n' + \
               '+ {:5d} {:5d} {:5d} {:5d} {:5d}\n' + \
               '- {:5d} {:5d} {:5d} {:5d} {:5d}'
    print(template.format(
           'A','C','G','T','N/*',
           base_n['A'], base_n['C'], base_n['G'], base_n['T'], base_n['N']+base_n['n'],
           base_n['a'], base_n['c'], base_n['g'], base_n['t'], base_n['*'],
           ))

    return 0


def indel(chrm, pos, bam):

    id_bam = os.path.basename(bam)+'.id.'+chrm+'_'+str(pos)+'.bam'
    cmd1 = "samtools view -h {bam} {chrm}:{start}-{end} | awk '($0 ~ /^@/) || ($6 ~ /[ID]/)' | samtools view -Sb - > {id_bam} ; \
            samtools index {id_bam}".format(bam=bam, chrm=chrm, start=pos-10, end=pos+10, id_bam=id_bam)
    sps.call(cmd1, shell=True)
    cmd2 = "export COLUMNS=181; samtools tview {bam} {hg19} -p {chrm}:{start} -d t ".format(
            bam=id_bam, hg19=GENOME, chrm=chrm, start=pos-90)

    run = sps.run(cmd2, shell=True, stdout=sps.PIPE)
    run.check_returncode()
    tv_lines = run.stdout.decode('utf-8')
    print(tv_lines)

    sps.call('rm {id_bam} {id_bam}.bai'.format(id_bam=id_bam), shell=True)

    return 0


def sv(chrm, pos, bam):

    sc_bam = os.path.basename(bam)+'.sc.'+chrm+'_'+str(pos)+'.bam'   # soft clip.
    cmd1 = "{samtools} view -h {bam} {chrm}:{start}-{end} | awk '($0 ~ /^@/) || ($6 ~ /S/)' | {samtools} view -Sb - > {outbam} ; {samtools} index {outbam}".format(
            samtools=SAMTOOLS, bam=bam, chrm=chrm, start=pos-10, end=pos+10, outbam=sc_bam)
    sps.call(cmd1, shell=True)
    cmd2 = "export COLUMNS=181; {samtools} tview {bam} {hg19} -p {chrm}:{start} -d t ".format(
            samtools=SAMTOOLS, bam=sc_bam, hg19=GENOME, chrm=chrm, start=pos-90)

    run = sps.run(cmd2, shell=True, stdout=sps.PIPE)
    run.check_returncode()
    tv_lines = run.stdout.decode('utf-8')
    print(tv_lines)

    return 0


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('func', choices=['snv','indel','sv'], nargs=1,
                        help='"snv" or "indel" or "sv"')
    parser.add_argument('bam', nargs=1, help='bam file')
    parser.add_argument('loc', nargs=1, help='eg: chr1:10000, 1-based')
    parser.add_argument('-a', '--all', action='store_true',
                        help='enable to view alt (default) and ref reads, snv only')
    parser.add_argument('-g','--genome', nargs='?', help='reference fasta, like: /path/to../hg19.fa')
    parser.add_argument('-s', '--samtools', nargs='?', help='samtools bin')
    parser.add_argument('--output', nargs='?', help='text file write into')

    args = vars(parser.parse_args())

    func = args['func'][0]
    bam = args['bam'][0]
    chrm, pos = args['loc'][0].split(':')
    pos = int(pos)
    view_all = args['all']    # default = False

    global SAMTOOLS
    if args['samtools']:
        SAMTOOLS = args['samtools'][0]
    else:
        SAMTOOLS = shutil.which('samtools')
        if SAMTOOLS is None:
            print("Error, set samtools first if not passed explicitly.")
            return 1

    global GENOME
    if args['genome']:
        GENOME = args['genome'][0]
    else:
        GENOME = os.getenv('hg19')
        if GENOME is None:
            print('Error, env hg19 should be set first if not passed explicitly.')
            return 1

    print(func)
    print(bam)
    print(chrm, pos)
    print(SAMTOOLS)
    print(GENOME)

    if func == 'snv':
        snv(chrm, pos, bam, view_all)
    elif func == 'indel':
        indel(chrm, pos, bam)
    elif func == 'sv':
        sv(chrm, pos, bam)

    return 0


if __name__ == '__main__':
    main()

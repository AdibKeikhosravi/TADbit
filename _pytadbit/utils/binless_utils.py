"""
19 april 2018

Wrapper for binless normalization https://github.com/3DGenomes/binless by Spill YG


"""

from subprocess import Popen, PIPE
from os import path

from numpy import genfromtxt

from pytadbit.utils.file_handling import which

def binless_normalization(interaction_files, tmp_dir='.', **kwargs):

    script_path = which('normalize_binless.R')
    proc_par = ["Rscript", "--vanilla", script_path]
    out_bias_csv = path.join(tmp_dir,'biases.csv')
    out_decay_csv = path.join(tmp_dir,'decay.csv')
    out_distance_csv = path.join(tmp_dir,'distance.csv')
    rdata = path.join(tmp_dir,'binless_optimized.RData')
    configfile = path.join(tmp_dir,'config_binless.r')
    with open(configfile, "w") as output:
        output.write('''action = 'normalize'\n''')
        output.write('''setwd('%s')\n''' % path.abspath(tmp_dir))
        output.write('''infiles = c(%s)\n'''
                     % ', '.join('"{0}"'.format(path.abspath(w)) for w in interaction_files))
        for key, val in kwargs.items():
            if key == "read_lens":
                output.write('''read_lens = c(%s)\n'''
                             % ', '.join('{0}'.format(sk) for sk in val))
                continue
            output.write('''%s = '%s'\n''' % (key,val))
        output.write('''output_bias = '%s'\n''' % path.abspath(out_bias_csv))
        output.write('''output_decay = '%s'\n''' % path.abspath(out_decay_csv))
        output.write('''output_distance = '%s'\n''' % path.abspath(out_distance_csv))
        output.write('''rdata = '%s'\n''' % path.abspath(rdata))
    proc_par.append(configfile)
    proc = Popen(proc_par, stderr=PIPE)
    err = proc.stderr.readlines()
    print '\n'.join(err)

    biases_binless = genfromtxt(out_bias_csv, delimiter=',', dtype=float)
    decay_binless = genfromtxt(out_decay_csv, delimiter=',', dtype=float)
    distance_binless = genfromtxt(out_distance_csv, delimiter=',', dtype=float)
    decay = dict((d,v) for d,v in zip(distance_binless,decay_binless))

    return biases_binless, decay, rdata

def binless_signal_detection(rdata, action, tmp_dir='.', **kwargs):
    
    script_path = which('normalize_binless.R')
    proc_par = ["Rscript", "--vanilla", script_path]
    out_signal_csv = path.join(tmp_dir,'signal.csv')
    configfile = path.join(tmp_dir,'config_binless.r')
    with open(configfile, "w") as output:
        output.write('''action = '%s'\n''' % action)
        output.write('''setwd('%s')\n''' % path.abspath(tmp_dir))
        output.write('''rdata = '%s'\n''' % path.abspath(rdata))
        for key, val in kwargs.items():
            output.write('''%s = '%s'\n''' % (key,val))
        output.write('''output_signal = '%s' ''' % path.abspath(out_signal_csv))
    proc_par.append(configfile)
    proc = Popen(proc_par, stderr=PIPE)
    err = proc.stderr.readlines()
    print '\n'.join(err)

    return path.abspath(out_signal_csv)
# Copyright 2015 CRI lab at PVAMU. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================

"""Pipline on Preprocessing RNAseq"""

import os
import sys
import getopt
import argparse
import multiprocessing


__doc__ = 'Pipline for Preprocessing RNAseq DATA'
__version__ = 'version 0.01'

SALMON = "/home/xishuang/Kim/salmon/bin/salmon"
SAMTOOLS = "/home/xishuang/Kim/samtools/samtools"
MAX_PROCESSES = 10

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def RNAseqPrerpocessing(input_dir, input_format, output_dir, tool, index_file):
    if input_format == 'fastq':
        quantify(input_dir, output_dir, tool, index_file)
    elif input_format == 'bam':
        if not os.path.exists('./fastq'):
            os.mkdir('./fastq')
            convert_bam_to_fastq(input_dir, './fastq')
            quantify('./fastq', output_dir, tool, index_file)

def merge_salmon_quant_result(output_dir):
    file_vector = []
    TPM = []
    for dirs, files in os.walk(output_dir):
        count = 0
        for file_name in files:
            if file_name == 'quant.sf'
                file_path = os.path.join(dirs, file_name)
                fc = open(file_path)
                d = fc.readlines()
                count =  count + 1
                if count == 1:
                    line_num = 0
                    for line in d:
                        line_num = line_num + 1
                        v = line.split('\t')
                        tmp_v = []
                        tmp_v.append(v[0])
                        tmp_v.append(v[1])
                        tmp_v.append(v[2])
                        if line_num == 1:
                            tmp_v.append('TPM' + str(count))
                        else:
                            tmp_v.append(v[3])
                        TPM.append(tmp_v)
                else:
                    line_num = 0
                    for line in d:
                        line_num = line_num + 1
                        v = line.split('\t')
                        if line_num == 1:
                            TPM[line_num - 1].append('TPM' + str(count))
                        else:
                            TPM[line_num - 1].append(v[3])



    return 0

def convert_bam_to_fastq(input_dir, tmp_dir):
    for dirs, files in os.walk(self.input_dir):
        for file_name in files:
            input_path = os.path.join(input_dir, file_name)
            output_path = os.path.join(tmp_dir, os.path.splitext(file_name)[0] + ".fastq")
            cmd = SAMTOOLS + ' fastq ' + input_path + ' > ' + output_path
            os.system(cmd)


def quantify(input_dir, output_dir, tool, index_file):
    if self.tool == 'salmon':
        self.salmon_quantify_mp(input_dir, output_dir, index_file)
    
def salmon_quantify(input_dir, output_dir, index_file,file_name):
    input_file_path = os.path.join(input_dir, file_name)
    quant_file = os.path.splitext(file_name)[0] + ".quant"
    cmd = SALMON + " quant -i " + index_file +" -l A -r " + input_file_path + " -o " + output_dir + "/" + quant_file
    os.system(cmd)

def salmon_quantify_mp(input_dir, output_dir, index_file):
    for root, dirs, files in os.walk(self.input_dir):
        pool = multiprocessing.Pool(processes = MAX_PROCESSES)
        for file_name in files:
            pool.apply_async(salmon_quantify, (input_dir, output_dir, index_file,file_name,))
        pool.close()
        pool.join()
        
        print 'Quantification finished!'
    return 0


def main(argv=None): 
    
    input_dir = ''
    output_dir = ''
    index_file = ''
    input_format = ''
    tool = ''

    if argv is None:
        argv = sys.argv
        print "Argument List: "
        print " -h --help"
        print " -v --version"
        print " -i --input_dir"
        print " -f --input_format Format: .bam/.fastq default: fastq"
        print " -o --output_dir"
        print " -d --index_file"
        print " -q --quantification Tool: Salmon/Kollisto/RSEM default: Salmon"
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "vhi:o:d:q:f:", ["help","input=","output="])
        except getopt.error, msg:
             raise Usage(msg)
            # process options
        for opt, arg in opts:
            if opt in ("-h", "--help"):
                print __doc__
                sys.exit(0)
            elif opt in ("-v", "--version"):
                print __version__
                sys.exit(0)
            elif opt in ("-q", "--quantification"):
                tool = arg
            elif opt in ("-i", "--input_dir"):
                input_dir = arg
            elif opt in ("-o", "--output_dir"):
                output_dir = arg
            elif opt in ("-d", "--index_file"):
                index_file = arg
            elif opt in ("-f", "--input_format"):
                input_format = arg

        if input_dir == '' or output_dir == '' or input_format == '' or tool == '' or index_file == '':
            print __doc__
            sys.exit(0)

        quantify(input_dir, input_format, output_dir, tool, index_file)


    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "for help use --help"
        

if __name__ == '__main__':
    main()

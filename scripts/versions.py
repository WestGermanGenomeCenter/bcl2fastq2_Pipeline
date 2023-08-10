#!/usr/bin/env python

import yaml
import glob
import os
import sys

# this file needs to be edited accordingly in case new features are added
# what mainly needs to be added are specific indirect tools and parts that should be read explicitly


# specify tools that may be in the yamls but are mostly supporting libraries or software
indirect_tools = ['perl',
                  'gffutils',
                  'pandas',
                  'libiconv',
                  'libgcc',
                  'zlib']

yamls = list()
software = list()

snakemakefolder = os.path.dirname(sys.argv[1])

def read_yaml(file):
    with open(file, 'r') as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

def get_yamls(directory):
    yamls = glob.glob(os.path.join(directory,"envs/*.yaml"))
    return yamls

def get_config(directory):
    config = glob.glob(os.path.join(directory,"config.yaml"))
    return config

def cutadapt_param(adapters, adapter_type):
    params = ""
    for i in range(len(adapters)):
        params+="-"+str(adapter_type[i])+" "+str(adapters[i])+" "
    return params

def pairwise(iterable):
    "s -> (s0, s1), (s2, s3), (s4, s5), ..."
    a = iter(iterable)
    return zip(a, a)

for file in get_yamls(snakemakefolder):
    list = read_yaml(file)["dependencies"][:]
    for i in list:
        out = i.split("=")[:2]
        if out[0].strip() in indirect_tools:
            continue
        print("Software: %s v%s" % (out[0].strip(), out[1]))
        software.append(out[0].strip())
        

print("")

for x in get_config(snakemakefolder):
    for y in read_yaml(x):
        if y in software or "bcl2fastq":
            if 'bcl2fastq' == y:
                print("bcl2fastq2 lief mit folgenden Parametern: barcode_mismatch %s %s" % (read_yaml(x)[y]['barcode_mismatch'], read_yaml(x)[y]['other_params']))
            elif 'cutadapt' == y:
                if read_yaml(x)[y]['cutadapt_active']:
                    print("cutadapt lief mit folgenden Parametern: %s %s" % (cutadapt_param(read_yaml(x)[y]['adapters'],read_yaml(x)[y]['adapter_type']),read_yaml(x)[y]['other_params']))
            else:
                 if 'other_params' in read_yaml(x)[y]:
                     active = [s for s in read_yaml(x)[y] if "active" in s]
                     if read_yaml(x)[y][active[0]]:
                         print("%s lief mit folgenden Parametern: %s" % (y, read_yaml(x)[y]['other_params']))
         
print("")
print("Alle sonstigen Tools die keine speziellen Parameterangaben haben wurden dementsprechend standardmäßig ohne besondere Parameter aufgerufen")


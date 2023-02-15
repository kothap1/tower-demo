#!/usr/bin/env python3

# Purpose          : Accepts a text file containing paths to JSON manifests for samples to run.
#                    Reads these manifests and extract sample information in a CSV file having
#                    3 columns: sampleID, read1 and read2. If run is executed for per-lane
#                    FASTQs, sampleID is <sample-name>_<lane_numbers>
# Usage            : parse_manifests.py [-h] --samples FILE [--output FILE] --processing STR
# Alternate usage  : parse_manifests.py [-h] --samplename FILE --manifestdir STR [--output FILE] --processing STR
# Date             : 2022

import argparse
import json
import os, sys
# import pandas as pd
from collections import defaultdict

# Get all manifests for a sample from a specified folder
def getSampleManifests(samplename, manifestFolder):
    sampleManifests = []
    for fn in os.listdir(manifestFolder):
        if not fn.endswith('.json'): continue
        manifest_fn = os.path.join(manifestFolder, fn)
        fastq_fn = json.load(open(manifest_fn))['file']['file_name']
        if not samplename in fastq_fn: continue
        sampleManifests.append(manifest_fn)
    return sampleManifests

# Read from input txt file
def readInput(filename):
    manifests = []
    with open(filename) as f:
        for line in f:
            manifests.append(line.strip())
    return(sorted(manifests))

# Read manifests for the entries in the txt file
def readManifests(manifests):
    collated = []
    sample_dir = os.path.dirname(manifests[1]) #TODO: tmp fix
    for json_file in manifests:
        content = json.load(open(json_file))
        file_path = content['file']['file_name']
        flow_cell = content['read_group']['flow_cell_barcode']

        base_dir = os.path.dirname(file_path)
        fastq_name = os.path.basename(file_path)
        sample_name = '_'.join(fastq_name.split('_')[:-3])
        lane_num = fastq_name.split('_')[-3]
        read_pair = fastq_name.split('_')[-2]
        extension = fastq_name.split('_')[-1]
        sample_info = {
        f'{flow_cell}/{sample_name}_{lane_num}': {
            'fastq_name': fastq_name,
            'sample_name':  sample_name,
            'lane_num': lane_num,
            'read_pair': read_pair,
            'extension': extension,
            'base_dir': base_dir,
            'flow_cell': flow_cell
            }
        }
        collated.append(sample_info)
    return collated, sample_dir

# Generate TSV from the data extracted from manifests
def generateTSV(collated, sample_dir, output_file, processing_way, fastq_rootdir):
    sample_grp = defaultdict(lambda: [])
    for x in collated:
        for k, v in x.items():
            if v not in x.values():
                sample_grp[k] = v
            else:
                sample_grp[k].append(v)

    with open(output_file,'w') as f:
        f.write('sampleID\tread1\tread2\n')
        for key, value in sample_grp.items():
            for info in value:
                if info['read_pair'] == 'R1':
                    read1 = os.path.join(fastq_rootdir, f"{info['base_dir']}/{info['fastq_name']}")
                elif info['read_pair'] == 'R2':
                    read2 = os.path.join(fastq_rootdir, f"{info['base_dir']}/{info['fastq_name']}")
            if processing_way == 'per-lane':
                sample_name = os.path.basename(key)
                f.write(f'{sample_name}\t{read1}\t{read2}\n')
            elif processing_way == 'per-sequencer':
                base_name = os.path.basename(key)
                flowcell = os.path.basename(os.path.dirname(key))
                sample_name = '_'.join(base_name.split('_')[:-1])
                f.write(f'{sample_name}_{flowcell}\t{read1}\t{read2}\n')
            else:
                sample_name = os.path.basename(key)
                flowcell = os.path.basename(os.path.dirname(key))
                f.write(f'{sample_name}_{flowcell}\t{read1}\t{read2}\n')

# get flowcells from redsheet
def get_flowcells(sample, redsheetfile):
    flowcells = redsheetfile.loc[sample,'flowcell_layout'].split(';')
    return flowcells

def getSamplesFromDirectory(samplename, sample_dir, flowcells):
	collated = []
	for flowcell in flowcells:
		for f in os.listdir(f'{sample_dir}/{flowcell}'):
			if f.endswith("fastq.gz"):
				if samplename in f:
					base_dir = f'{sample_dir}/{flowcell}'
					fastq_name = f
					sample_name = '_'.join(fastq_name.split('_')[:-3])
					lane_num = fastq_name.split('_')[-3]
					read_pair = fastq_name.split('_')[-2]
					extension = fastq_name.split('_')[-1]
					sample_info = {
					f'{flowcell}/{sample_name}_{lane_num}': {
						'fastq_name': fastq_name,
						'sample_name':  sample_name,
						'lane_num': lane_num,
						'read_pair': read_pair,
						'extension': extension,
						'base_dir': flowcell,
						'flow_cell': flowcell
						}
					}
					collated.append(sample_info)
	return collated

def get_sample_names(redsheetfile, key):
    samples = []
    headers = dict()
    for i, line in enumerate(open(redsheetfile).readlines()):
        if i<7: continue
        if i == 7:
            for j,k in enumerate(line.strip().split(',')):
                headers[k] = j
        line = line.split(',')
        if len(line)<8: continue
        if 'Control' in line[headers['eln_sample_id']]: continue
        if line[headers['analysis_type']] != 'Whole Genome Genomic': continue
        samples.append(line[headers[key]])
    return samples

print(sys.argv)
# Command line arguments -- {{{
parser = argparse.ArgumentParser(description='Parses reads information from manifest JSON files and writes parsed information in a CSV file.')
parser.add_argument('--samples', help='txt file containing path to JSON manifests to parse.', type=str, metavar='FILE')
parser.add_argument('--samplename', help='sample name pattern to identify manifests for.', type=str, metavar='STR')
parser.add_argument('--manifestdir', help='path to folder containing manifest files.', type=str, metavar='FILE')
parser.add_argument('--fastq_rootdir', help='path to folder containing fastqs folder.', type=str)
parser.add_argument('--sampledir', help='path to folder containing samples.', type=str, metavar='FILE')
parser.add_argument('--redsheet', help='redsheet for samples.', type=str, metavar='FILE')
parser.add_argument('--output', help='path to output TSV file', type=str, metavar='FILE', default='samples.csv')
parser.add_argument('--processing', help='per-lane or per-sequencer', type=str, metavar='STR', required=True)
args = parser.parse_args()
# }}} --


# Call to functions --{{{
if args.samplename and args.manifestdir:
	manifests = getSampleManifests(args.samplename, args.manifestdir)
	collated, sample_dir = readManifests(manifests)
elif args.redsheet and args.manifestdir:
	# redsheetfile = pd.read_csv(f'{args.redsheet}', skiprows=7, sep=',')
	# redsheetfile = redsheetfile[(redsheetfile.analysis_type=='Whole Genome Genomic')&~(redsheetfile.eln_sample_id.str.contains('Control'))]
	print(args.manifestdir, args.redsheet, args.fastq_rootdir)
	manifests = []
	for sample in get_sample_names(args.redsheet, 'name_alias'): # redsheetfile.name_alias:
		manifests += getSampleManifests(sample, args.manifestdir)
	collated, sample_dir = readManifests(manifests)
# elif args.redsheet and args.sampledir:
# 	sample_dir = args.sampledir
# 	redsheetfile = pd.read_csv(f'{args.redsheet}', skiprows=7, sep=',', index_col='name_alias')
# 	redsheetfile = redsheetfile[(redsheetfile.analysis_type=='Whole Genome Genomic')&~(redsheetfile.eln_sample_id.str.contains('Control'))]
# 	collated = []
# 	for sample in redsheetfile.index:
# 		flowcells = get_flowcells(sample, redsheetfile)
# 		collated += getSamplesFromDirectory(sample, args.sampledir, get_flowcells(sample, redsheetfile))
# else:
# 	manifests = readInput(args.samples)
# 	collated, sample_dir = readManifests(manifests)

generateTSV(collated, sample_dir, args.output, args.processing, args.fastq_rootdir)
# }}}--
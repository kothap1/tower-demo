#!/usr/bin/env python3


#TODO: cleanup and add comments

import os
import json
import pandas as pd
import argparse

def getSampleManifests(samplename, manifestFolder):
	sampleManifests = []
	for fn in os.listdir(manifestFolder):
		if not fn.endswith('.json'): continue
		manifest_fn = os.path.join(manifestFolder, fn)
		fastq_fn = json.load(open(manifest_fn))['file']['file_name']
		if not samplename in fastq_fn: continue	
		sampleManifests.append(manifest_fn)
	return sampleManifests[0]

def getFastqName(json_file):
	content = json.load(open(json_file))
	fastq_name = os.path.basename(content['file']['file_name'])
	sample_name = '_'.join(fastq_name.split('_')[:-3])
	return sample_name

def generateBamManifest(redsheetfile, redsheetfilename):
	for sample in redsheetfile.index:
		manifest_json = redsheetfile.loc[sample, 'manifest_json']
		with open(manifest_json) as f:
			json_decoded = json.load(f)
			json_decoded['file']['file_name'] = redsheetfile.loc[sample, 'bam']
			json_decoded['file']['description'] = f"RED Sheet: {redsheetfilename} | Flowcell: {redsheetfile.loc[sample, 'flowcell_layout']} | Sample Type: {redsheetfile.loc[sample, 'sample_type']}"
			json_decoded['user_id'] = redsheetfile.loc[sample, 'user_id']

			json_decoded['alignment_workflow'] = {
			"submitted_aligned_reads_files":f"{redsheetfile.loc[sample, 'bam']}",
			"workflow_link":"sambamba markdup",
			"submitter_id":"ECD Pipeline v1.0",
			"submitted_unaligned_reads_files": f"{redsheetfile.loc[sample, 'fastq']}",
			"type":"alignment_workflow"
   			}
			json_decoded['aligned_reads'] = {
			"data_format":"BAM",
			"data_type":"Aligned Reads",
			"data_category":"Sequencing Reads",
			"experimental_strategy":"WGS",
			}
			
			del json_decoded['unaligned_reads']
			del json_decoded['read_group']
			if 'file_id' in json_decoded['file']:
				del json_decoded['file']['file_id']
			if 'file_size' in json_decoded['file']:
				del json_decoded['file']['file_size']
			
			with open(f"{redsheetfile.loc[sample, 'bam']}.json", 'w') as f:
				f.write(json.dumps(json_decoded, indent=4, separators=(',', ': ')))
	return

def generateBaiManifest(redsheetfile):
	for bamjson in os.listdir("."):
		if bamjson.endswith('.bam.json'):
			for sample in redsheetfile.index:
				if sample in bamjson:
					with open(bamjson) as f:
						json_decoded = json.load(f)
						json_decoded['file']['file_name'] = redsheetfile.loc[sample, 'bai']

						json_decoded['aligned_reads_index'] = {
							"data_format":"BAI",
							"submitted_aligned_reads_files": redsheetfile.loc[sample, 'bam'],
							"data_category":"Sequencing Data",      
							"data_type":"Aligned Reads Index"
						}

						del json_decoded['aligned_reads']

						with open(f"{redsheetfile.loc[sample, 'bai']}.json", 'w') as f:
							f.write(json.dumps(json_decoded, indent=4, separators=(',', ': ')))
	return

def accessoryFilesTemplate(json_decoded, file_name, description, file_type, process):
	template = {
		"file":{
      		"file_name": file_name,
      		"description": description,
      		"file_type": file_type
		},
		"case":{
      		"case_submitter_id": json_decoded['case']['case_submitter_id'],
      		"case_id": json_decoded['case']['case_id'],
   		},
		"sample":{
			"sample_submitter_id": json_decoded['sample']['sample_submitter_id'],
			"sample_type": json_decoded['sample']['sample_type'],
			"date_of_collection": json_decoded['sample']['date_of_collection'],
			"sample_id": json_decoded['sample']['sample_id']
   		},
		"aliquot":{
			"analyte_type": json_decoded['aliquot']['analyte_type'],
			"aliquot_volume": json_decoded['aliquot']['aliquot_volume'],
			"source_center": json_decoded['aliquot']['source_center'],
			"aliquot_quantity": json_decoded['aliquot']['aliquot_quantity'],
			"aliquot_submitter_id": json_decoded['aliquot']['aliquot_submitter_id'],
			"aliquot_id": json_decoded['aliquot']['aliquot_id']
		},
		"clinical_supplement":{
      		"data_format": "TSV",
      		"data_type":"Clinical Supplement",
      		"data_category":"Clinical"
   		},
		"alignment_workflow":{
			"submitted_aligned_reads_files": json_decoded['alignment_workflow']['submitted_aligned_reads_files'],
			"workflow_link": process,
			"submitter_id": json_decoded['alignment_workflow']['submitter_id'],
			"submitted_unaligned_reads_files": json_decoded['alignment_workflow']['submitted_unaligned_reads_files'],
			"type": json_decoded['alignment_workflow']['type']
		},
		"type_of_operation": json_decoded['type_of_operation'],
		"file_category":"RAW-NONSEQ",
		"id":"manifest-raw-non-seq.json",
		"user_id": json_decoded['user_id']
		}
	return template

def generateTXTManifest(redsheetfile):
	for bamjson in os.listdir("."):
		if bamjson.endswith('.bam.json'):
			for sample in redsheetfile.index:
				if sample in bamjson:
					with open(bamjson) as f:
						json_decoded = json.load(f)

						# Fastp
						for lane in ['L001','L002','L003','L004']:
							file_name = f"{redsheetfile.loc[sample, 'fastq']}_{lane}.json"
							description = f"Output from fastp with the following options - Lane {lane}"
							template = accessoryFilesTemplate(json_decoded, file_name, description, ".json", "fastp")
							with open(f"{redsheetfile.loc[sample, 'fastq']}_{lane}.json.json", 'w') as f:
								f.write(json.dumps(template, indent=4, separators=(',', ': ')))

						# mosdepth
						for mosdepthreport in ['mosdepth.global.dist','mosdepth.summary']:
							file_name = f"{redsheetfile.loc[sample, 'fastq']}.{mosdepthreport}.txt"
							description = f"Output from mosdepth - {mosdepthreport}"
							template = accessoryFilesTemplate(json_decoded, file_name, description, ".txt", "mosdepth")
							with open(f"{redsheetfile.loc[sample, 'fastq']}_{mosdepthreport}.txt.json", 'w') as f:
								f.write(json.dumps(template, indent=4, separators=(',', ': ')))

						# Picard InsertSizeMetrics
						file_name = f"{redsheetfile.loc[sample, 'fastq']}.insert_size_metrics.txt"
						description = f"Output from Picard - CollectInsertSizeMetrics"
						template = accessoryFilesTemplate(json_decoded, file_name, description, ".txt", "picard CollectInsertSizeMetrics")
						with open(f"{redsheetfile.loc[sample, 'fastq']}.insert_size_metrics.txt.json", 'w') as f:
							f.write(json.dumps(template, indent=4, separators=(',', ': ')))

						# Picard MultipleMetrics
						for picardreport in ['alignment_summary_metrics','gc_bias.detail_metrics','gc_bias.summary_metrics','quality_distribution_metrics']:
							file_name = f"{redsheetfile.loc[sample, 'fastq']}.{picardreport}"
							description = f"Output from Picard - CollectMultipleMetrics"
							template = accessoryFilesTemplate(json_decoded, file_name, description, ".txt", "picard CollectMultipleMetrics")
							with open(f"{redsheetfile.loc[sample, 'fastq']}.{picardreport}.json", 'w') as f:
								f.write(json.dumps(template, indent=4, separators=(',', ': ')))
	return

def generateHTMLXLXSmanifest(redsheetfile, redsheetfilename):
	for bamjson in os.listdir("."):
		if bamjson.endswith('.bam.json'):
			for sample in redsheetfile.index:
				if sample in bamjson:
					with open(bamjson) as f:
						json_decoded = json.load(f)

						# Collated QC: HTML report
						for sample in redsheetfile.index:
							file_name = f"batchqc_{redsheetfile.loc[sample, 'fastq']}.html"
							description = f"Output of QC collation (Jupyter Notebook) with the following options - RED sheet: {redsheetfilename}"
							template = accessoryFilesTemplate(json_decoded, file_name, description, ".html", "Jupyter Notebook")
							with open(f"batchqc_{redsheetfile.loc[sample, 'fastq']}.html.json", 'w') as f:
								f.write(json.dumps(template, indent=4, separators=(',', ': ')))

						# Collated QC: XLSX report
						for sample in redsheetfile.index:
							file_name = f"batchqc_{redsheetfile.loc[sample, 'fastq']}.xlsx"
							description = f"Output of QC collation (Jupyter Notebook) with the following options - RED sheet: {redsheetfilename}"
							template = accessoryFilesTemplate(json_decoded, file_name, description, ".xlsx", "Jupyter Notebook")
							with open(f"batchqc_{redsheetfile.loc[sample, 'fastq']}.xlsx.json", 'w') as f:
								f.write(json.dumps(template, indent=4, separators=(',', ': ')))

						# Collated QC: MultiQC report
						for sample in redsheetfile.index:
							file_name = f"multiqc_{redsheetfile.loc[sample, 'fastq']}.html"
							description = f"Output of MultiQC with the following options - RED sheet: {redsheetfilename}; module: fastp"
							template = accessoryFilesTemplate(json_decoded, file_name, description, ".html", "MultiQC")
							with open(f"multiqc_{redsheetfile.loc[sample, 'fastq']}.html.json", 'w') as f:
								f.write(json.dumps(template, indent=4, separators=(',', ': ')))
	return
	

# Command line arguments --- {{{
parser = argparse.ArgumentParser(description='Creates JSON manifests for the outputs of the pipeline.')
parser.add_argument('--samplename', help='sample name to create manifests for.', type=str, metavar='STR')
parser.add_argument('--manifestdir', help='path to folder containing manifest files.', type=str, metavar='FILE')
parser.add_argument('--redsheet', help='redsheet for containing the sample', type=str, metavar='FILE')
parser.add_argument('--userid', help='User running the pipeline', type=str, metavar='STR')
args = parser.parse_args()
# }}} --

redheetpath = args.redsheet
redsheetname = os.path.basename(redheetpath)
redsheetname_noext = os.path.splitext(redsheetname)[0]
manifestdir = args.manifestdir
sample =  '_'.join(args.samplename.split('_')[:-1])

redsheetfile = pd.read_csv(redheetpath, skiprows=7, sep=',',index_col='name_alias')
redsheetfile = redsheetfile[redsheetfile.index.str.contains(sample)]
#redsheetfile = redsheetfile[(redsheetfile.analysis_type=='Whole Genome Genomic')&~(redsheetfile.eln_sample_id.str.contains('Control'))]


for sample in redsheetfile.index:
	manifest_json = getSampleManifests(sample, manifestdir)
	sample_name = getFastqName(manifest_json)
	redsheetfile.loc[sample, 'manifest_json'] = manifest_json
	redsheetfile.loc[sample, 'fastq'] = sample_name
	redsheetfile.loc[sample, 'bam'] = f'{sample_name}.markeddup.bam'
	redsheetfile.loc[sample, 'bai'] = f'{sample_name}.markeddup.bam.bai'
	redsheetfile.loc[sample, 'user_id'] = args.userid

generateBamManifest(redsheetfile, redsheetname)
generateBaiManifest(redsheetfile)
generateTXTManifest(redsheetfile)
generateHTMLXLXSmanifest(redsheetfile, redsheetname)

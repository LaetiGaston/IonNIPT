#!/usr/bin/env python

import argparse
import os
import subprocess 
from joblib import Parallel, delayed
import pickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def getArgs():
	parser = argparse.ArgumentParser(description="",version="1.0.0")
	parser.add_argument('-b',dest="bamDir",type=str,required=True,help='BAM folder')
	parser.add_argument('-p',dest="pluginDir",type=str,required=True,help='Path to IonNIPT folder')
	parser.add_argument('-j',dest="nbJob",type=int,required=False,help='Number of jobs for parallelisation', default=24)
	parser.add_argument('-o',dest="outDir",type=str,required=True,help='Output folder')
	arg = parser.parse_args()
	
	return arg

def main(args):

	### get BAM files (allows parallelisation)
	bam_lst = []
	for bam in os.listdir(args.bamDir):
		if bam.endswith(".bam"):
			if not os.path.exists(os.path.join(bam,".bai")):
				sys.stdout.write("Bam files must be indexed")
				sys.exit(1)
			bam_lst.append(bam)
			
	### Make Wisecondor analysis (in //)
	Parallel(n_jobs=args.nbJob)(delayed(wisecondor)(bam, args.bamDir, args.outDir, args.pluginDir) for bam in bam_lst)
	
	## Get smooth zscores
	z_wisecondor = {}
    	for bam in bam_lst:
        	name = os.path.basename(bam).split('.')[0]
        	z_wisecondor[name] = {}
        	with open(os.path.join("%s"%args.outDir,"%s.tested"%name)) as file:
            		data = pickle.loads(file.read())
            		for chrom in ["13","18","21"]:
                		zScores = data["zSmoothDict"][chrom]
                		score   = -1
                		try:
                    			score = sum(zScores) / len(zScores)
                		except:
                    			score = -1
                		z_wisecondor[name][chrom] = str(round(score,2))
	
	### Make Sanefalcon analysis
	readstarts = os.path.join(args.outDir, 'readstarts')
	os.makedirs(readstarts)
	
	## Get readstarts for each bam, in // per chromosome
	nbJob = args.nbJob>=4 and args.nbJob/4 or 1
	for bam in bam_lst:	
		Parallel(n_jobs=nbJob)(delayed(getReadStarts)(bam, id, args.bamDir, readstarts, args.pluginDir) for id in range(1,23))
	
	## Get nucleosome profiles
	profiles = = os.path.join(args.outDir, 'profiles')
	os.makedirs(profiles)
	
	Parallel(n_jobs=args.nbJob)(delayed(getProfiles)(readstart, args.pluginDir, readstarts, profiles) for readstart in os.listdir(readstarts))
			
	## Get foetal fraction
	
	predict = os.path.join(args.pluginDir, 'sanefalcon/predict.sh')
	trainModel = os.path.join(args.pluginDir, 'data/trainModel_BorBreCoc_defrag-rassf1.model')
	
	ff_sanefalcon = {}
	for bam in bam_lst:
		name = os.path.basename(bam).split('.')[0]
		bam = os.path.join(profiles, name)
		
		ff_cmd = "bash {predict} {trainModel} {bam}".format(predict=predict, trainModel=trainModel, bam=bam)
		jobLauncher(ff_cmd)
		
		ff = open(bam + '.ff','r')
		for line in ff:
			elem = line.split()
			if elem[0] == 'Fetal':
				ff_sanefalcon[name] = round(float(elem[2]), 2)
		
	#### defrag
	
	###
	
	defrag = os.path.join(args.pluginDir, 'wisecondor/defrag.py')
	male = os.path.join(args.pluginDir, 'data/male-defrag')
	female = os.path.join(args.pluginDir, 'data/female-defrag')
	sf = 0.688334125062
	percYonMales = 0.00146939199267
	
	###
	
	defrag_cmd = 'python {defrag} {male} {female} {test} --scalingFactor {sf} --percYonMales {percYmales} {graph}'.format(defrag=defrag, male=male, female=female , test=args.outDir, sf=sf, percYmales=percYonMales, graph="defrag")
	p = subprocess.Popen(defrag_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	stdout, stderr = p.communicate()
	if p.returncode != 0:
		raise Exception(stderr)

	## Get FF
	ff_defrag = {}
    	ff = open(outfile,'r')
    	for line in ff:
        	elem = line.split()
        	if elem[3] == "male" and elem[5] == "boys":
            		ff_defrag[elem[0]] = [round(float(elem[1]), 2), "male"]
        	elif elem[3] == "female" and elem[5] == "girls":
            		ff_defrag[elem[0]] = [round(float(elem[1]), 2), "female"]
        	else:
			# TODO: take into account coverage on specific genes on ChrY
            		ff_defrag[elem[0]] = [round(float(elem[1]), 2), "undetermined"]
		
	### Mimic InoNIPT plugin output on stdout
	print "sample\tchrom21\tchrom18\tchrom13\tff defrag\tff sanefalcon\tsex"
	for sample in ff_defrag.keys():
        print "{sample}\t{chr21}\t{chr18}\t{chr13}\t{ffdefrag}\t{ffsanefalcon}\t{sex}".format(sample=sample, chr21=z_wisecondor[sample]["21"], chr18=z_wisecondor[sample]["18"], chr13=z_wisecondor[sample]["13"], ffdefrag=ff_defrag[sample][0], ffsanefalcon=ff_sanefalcon[sample], sex=ff_defrag[sample][1])
	
def wisecondor(bam, input, output, plugin):
	
	### exe + data
	consampy = os.path.join(plugin, 'wisecondor/consam.py')
	testpy = os.path.join(plugin, 'wisecondor/test.py')
	plotpy = os.path.join(plugin, 'wisecondor/plot.py')
	gccpy = os.path.join(plugin, 'wisecondor/gcc.py')
	gccount = os.path.join(plugin, 'data/hg19.gccount')
	reftable =  os.path.join(plugin, 'data/reftable')
	
	### path + output names
	
	name = os.path.basename(bam).split('.')[0]
	bam = os.path.join(input, bam) # -> crush bam into new bam var, full path
	
	pickle = os.path.join(output, name + '.pickle')
	gcc = os.path.join(output, name + '.gcc')
	tested = os.path.join(output, name + '.tested')
	pdf = os.path.join(output, name)
		
	### wisecondor analysis
	
	pickle_cmd = "samtools view {bam} -q 1 | python {consampy} -outfile {pickle}".format(bam = bam, consampy = consampy, pickle = pickle)
	jobLauncher(pickle_cmd)
	
	gcc_cmd = "python {gccpy} {pickle} {gccount} {gcc}".format(gccpy = gccpy, pickle = pickle, gccount = gccount, gcc = gcc)
	jobLauncher(gcc_cmd)
	
	tested_cmd = "python {testpy} {gcc} {reftable} {tested}".format(testpy =  testpy, gcc = gcc, reftable = reftable, tested = tested)
	jobLauncher(tested_cmd)
	
	pdf_cmd = "python {plotpy} {tested} {pdf}".format(plotpy = plotpy, tested = tested, pdf = pdf)
	jobLauncher(pdf_cmd)

def getReadStarts(bam, id, input, output, plugin):
	
	name = os.path.basename(bam).split('.')[0]
	bam = os.path.join(input, bam)
	
	forward = os.path.join(output, name + '.' + str(id) + '.start.fwd')
	reverse = os.path.join(output, name + '.' + str(id) + '.start.rev')
	
	###
	
	retro = os.path.join(plugin, 'sanefalcon/retro.py')
	
	###
	
	cmd_forward = "samtools view {bam} chr{id} -F 20 -q 1 | python {retro} | awk ".format(bam = bam, retro = retro, id = id)
	cmd_forward += "'{print $4}' > "
	cmd_forward += "{forward}".format(forward = forward)
	
	cmd_reverse = "samtools view {bam} chr{id} -f 16 -F 4 -q 1 | python {retro} | awk ".format(bam = bam, retro = retro, id = id)
	cmd_reverse += "'{print ($4 + length($10) - 1)}' "
	cmd_reverse += "> {reverse}".format(reverse = reverse)
	
	jobLauncher(cmd_forward)
	jobLauncher(cmd_reverse)

def getProfiles(readstart, plugin, input, output):
	name, chr, type, strand = os.path.basename(readstart).split('.')
	
	###
	
	getProfile = os.path.join(plugin, 'sanefalcon/getProfile.py')
	nuclTrack = os.path.join(plugin, 'data')
	
	###
	
	if strand == 'fwd':
		cmd_0 = "python {getProfile} {nuclTrack}/nuclTrack.{chr} {input}/{name}.{chr}.start.fwd 0 {output}/{name}.{chr}.fwd".format(getProfile=getProfile, nuclTrack=nuclTrack, chr=chr, name=name, input=input, output=output)
		cmd_1 = "python {getProfile} {nuclTrack}/nuclTrack.{chr} {input}/{name}.{chr}.start.fwd 1 {output}/{name}.{chr}.ifwd".format(getProfile=getProfile, nuclTrack=nuclTrack, chr=chr, name=name, input=input, output=output)
		
		jobLauncher(cmd_0)
		jobLauncher(cmd_1)
		
	elif strand == 'rev':
		cmd_0 = "python {getProfile} {nuclTrak}/nuclTrack.{chr} {input}/{name}.{chr}.start.rev 1 {output}/{name}.{chr}.rev".format(getProfile=getProfile, nuclTrak=nuclTrak, chr=chr, name=name, input=input, output=output)
		cmd_1 = "python {getProfile} {nuclTrak}/nuclTrack.{chr} {input}/{name}.{chr}.start.rev 0 {output}/{name}.{chr}.irev".format(getProfile=getProfile, nuclTrak=nuclTrak, chr=chr, name=name, input=input, output=output)
		
		jobLauncher(cmd_0)
		jobLauncher(cmd_1)

def jobLauncher(cmd):
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	stdout, stderr = p.communicate()
	if p.returncode != 0:
		raise Exception(stderr)

if __name__ == '__main__':
	args = getArgs()
	main(args)

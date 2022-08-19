#!/usr/bin/python
#-*- coding: utf-8 -*-
#===============================================================================
#
#         FILE: ScanIndel.py
#
#        USAGE: ./ScanIndel.py -i sample.txt -c config.txt [opts]
#
#  DESCRIPTION: Indel detection for NGS data
#
#      OPTIONS: ---
# REQUIREMENTS: See README.md file
#         BUGS: ---
#        NOTES: 
#       AUTHOR: Rendong Yang (yang4414@umn.edu)
# ORGANIZATION:
#===============================================================================
import sys, os, time, subprocess, getopt
try: import pysam
except: sys.exit('pysam module not found.\nPlease install it before.')
try: import vcf
except: sys.exit('PyVCF module not found.\nPlease install it before.')
try: from Bio import SearchIO
except: sys.exit('Biopython module not found.\nPlease install it before.')
try: import numpy as np
except: sys.exit('numpy module not found.\nPlease install it before.')
try: from scipy.stats import binom
except: sys.exit('scipy module not found.\nPlease install it before.')

__version__ = '1.3'

def read_config_file(filename):
	path = {}
	ifp = open(filename)
	for line in ifp:
		if line[0] == '#' or line == '\n':
			continue
		item = line.rstrip().split('=')
		path[item[0]] = item[1]
	ifp.close()
	return path

def read_sample_file(filename):
	input = {}
	ifp = open(filename)
	for line in ifp:
		if line[0] == '#' or line == '\n':
			continue
		item = line.rstrip().split()
		input[item[0]] = ' '.join(item[1:])
	ifp.close()
	return input

def read_gfServer_log(infile):
	f = open(infile)
	for i in f:
		if 'Server ready for queries!' in i:
			f.close()
			return 1
		elif 'gfServer aborted' in i or 'error' in i:
			f.close()
			return -1
	f.close()
	return 0

def start_gfServer(reference,gfServer_port,output_dir):
	sys.stderr.write('Starting gfServer\n')
	cwd = os.getcwd()
	if os.path.isabs(output_dir):
		log_dir = output_dir
	else:
		log_dir = cwd+'/'+output_dir
	os.chdir(reference)
	try:
		subprocess.check_call('gfServer -canStop -log='+log_dir+'/gfserver.temp.log start localhost '+gfServer_port+' *.2bit &', shell=True)
	except subprocess.CalledProcessError as e:
		print("Execution failed for gfServer:", e, file=sys.stderr)
		sys.exit(1)
	os.chdir(cwd)
	time.sleep(5) # aspetta 5 sec
	gfready = 0
	while not gfready:
		gfready = read_gfServer_log(log_dir+'/gfserver.temp.log')
		if gfready == -1:
			print('gfServer start error! Check '+log_dir+'/gfserver.temp.log file, exit.')
			sys.exit(1)
		time.sleep(5)
	sys.stderr.write('gfServer ready\n')

def bwa_alignment(reference, input, output):
	try:
		subprocess.check_call("bwa mem -M -t8 "+reference['bwa']+" "+input+" >"+output+".bwa.sam", shell=True)
		subprocess.check_call("samtools view -bS "+output+".bwa.sam >"+output+".bwa.bam", shell=True)
		subprocess.check_call("samtools sort -o "+output+".bwa.sorted.bam "+output+".bwa.bam", shell=True)
		subprocess.check_call("mv "+output+".bwa.sorted.bam "+output, shell=True)
		subprocess.check_call("samtools index "+output, shell=True)
	except subprocess.CalledProcessError as e:
		print('Execution failed for BWA mem:', e, file=sys.stderr)
		sys.exit(1)

def psl2sam(hsp, query_seq_len):
	"""psl2sam try to implement the psl2sam.pl script and return the cigar and mapping position estimated from psl file"""
	cigar = ''
	query_start = hsp.query_start
	query_end = hsp.query_end
	strand = hsp.query_strand_all[0] # may need replace by qery_strand
	soft_len = 0
	if strand == -1:
		query_start = query_seq_len - hsp.query_end
		query_end = query_seq_len - hsp.query_start
	if query_start:
		# 5'-end clipping
		soft_len = query_start
		cigar += str(query_start)+'S'
	x = hsp.query_span_all
	if strand == -1:
		y = [query_seq_len - item[1] for item in hsp.query_range_all] # may need replace by query_start_all when the bug is fixed in Biopython
	else:
		y = [item[0] for item in hsp.query_range_all] # may need replace by query_start_all when the bug is fixed in Biopython
	z = hsp.hit_start_all
	y0, z0 = y[0], z[0]
	for i in range(1, len(hsp)):
		ly = y[i] - y[i-1] - x[i-1]
		lz = z[i] - z[i-1] - x[i-1]
		if ly < lz:
			# del: the reference gap is longer
			cigar += str(y[i] - y0) + 'M'
			cigar += str(lz - ly) + 'D'
			y0, z0 = y[i], z[i]
		elif lz < ly:
			# ins: the query gap is longer
			cigar += str(z[i] - z0) + 'M'
			cigar += str(ly - lz) + 'I'
			y0, z0 = y[i], z[i]
	cigar += str(query_end - y0) + 'M'
	if query_seq_len != query_end:
		# 3'-end clipping
		end3 = query_seq_len - query_end
		if end3 > soft_len:
			soft_len = end3
		cigar += str(end3) + 'S'
	return cigar, soft_len

def get_softclip_length(read):
	if not read.qual:
		qual = np.array([40] * read.rlen)
	else:
		qual = np.array(list(map(ord, list(read.qual))))
		qual = qual - 33
	if read.cigar[0][0] == 4:
		if read.cigar[-1][0] == 4:
			if read.cigar[0][1] > read.cigar[-1][1]:
				return read.cigar[0][1], qual[:read.cigar[0][1]], read.pos
			else:
				return read.cigar[-1][1], qual[read.rlen-read.cigar[-1][1]:], read.aend
		else:
			return read.cigar[0][1], qual[:read.cigar[0][1]], read.pos
	elif read.cigar[-1][0] == 4:
		return read.cigar[-1][1], qual[read.rlen - read.cigar[-1][1]:], read.aend
	else:
		return 0, qual, -1

def prob_of_indel_with_error(input, soft_chr, soft_pos, prob):
	alignment = pysam.Samfile(input,'rb')
	total = alignment.count(soft_chr,soft_pos,soft_pos+1)
	try:
		reads = [read for read in alignment.fetch(soft_chr, soft_pos - 1, soft_pos + 2)]
	except ValueError:
		reads = [read for read in alignment.fetch(soft_chr, soft_pos, soft_pos + 1)]
	num_soft = 0
	for each in reads:
		if each.is_secondary or each.is_unmapped:
			continue
		soft_len, soft_qual, soft_pos_read = get_softclip_length(each)
		if soft_len != 0 and abs(soft_pos_read - soft_pos) < 2: # +/- 1bp matching
			num_soft += 1
	return binom.sf(num_soft - 1, total, prob)

def blat_alignment(mapping, reference, scliplen_cutoff, lowqual_cutoff, min_percent_hq, mapq_cutoff, blat_ident_pct_cutoff, gfServer_port, hetero_factor, input, output):
	bwa_bam = pysam.Samfile(input, 'rb')
	blat_bam = pysam.Samfile(output + '.temp.bam', 'wb', template=bwa_bam)
	if hetero_factor != 'a':
		denovo = open(output+'.temp.fasta', 'w')
	putative_indel_cluster = set()
	for read in bwa_bam.fetch(until_eof=True):
		if read.is_secondary:
			# secondary alignment
			continue
		if read.is_unmapped:
			if hetero_factor != 'a':
				print('>' + read.qname, file=denovo)
				print(read.seq, file=denovo)
			continue
		soft_len, soft_qual, soft_pos = get_softclip_length(read)
		sclip_ratio = soft_len / float(read.rlen)
		if soft_pos != -1:
			sclip_hq_ratio = len(soft_qual[soft_qual >= lowqual_cutoff]) / float(len(soft_qual))
		else:
			sclip_hq_ratio = 0
		if sclip_ratio >= scliplen_cutoff and sclip_hq_ratio >= min_percent_hq and read.mapq >= mapq_cutoff:
			blat_aln = False
			soft_chr = bwa_bam.getrname(read.rname)
			if hetero_factor == 'a':
				blat_aln = True
			elif (soft_chr, soft_pos) in putative_indel_cluster:
				blat_aln = True
				print('>' + read.qname, file=denovo)
				print(read.seq, file=denovo)
				if not mapping:
					blat_bam.write(read)
					continue
			# estimate the probability of indels given the coverage and number of soft-clipping readss
			elif prob_of_indel_with_error(input, soft_chr, soft_pos, hetero_factor) < 0.05:
				putative_indel_cluster.add((soft_chr, soft_pos))
				blat_aln = True
				print('>' + read.qname, file=denovo)
				print(read.seq, file=denovo)
				if not mapping:
					blat_bam.write(read)
					continue
			if blat_aln:
				fa = open(output+'.temp.fa', 'w')
				print('>' + read.qname, file=fa)
				print(read.seq, file=fa)
				fa.close()
				try:
					subprocess.check_call('gfClient localhost ' + gfServer_port +' '+ reference['blat'] + ' ' + output + '.temp.fa ' + output + '.temp.psl >/dev/null 2>&1 ', shell=True)
				except subprocess.CalledProcessError as e:
					print('Execution failed for gfClient:', e, file=sys.stderr)
					sys.exit(1)
				try:
					blat = SearchIO.read(output+'.temp.psl', 'blat-psl')
					print('Blat aligned read:',read.qname, file=sys.stderr)
				except:
					print('No blat hit for read:',read.qname, file=sys.stderr)
					blat_bam.write(read)
					continue
				hsps = blat.hsps
				hsps.sort(key=lambda k: k.score, reverse=True)
				if hsps[0].ident_pct / 100 >= blat_ident_pct_cutoff and hsps[0].hit_id == bwa_bam.getrname(read.tid):
					# matching genomic coordinate
					if hsps[0].hit_start == read.pos or hsps[0].hit_end == read.aend:
						cigarstring, soft_len = psl2sam(hsps[0], blat.seq_len)
						if soft_len == 0:
							read.cigarstring, read.pos = cigarstring, hsps[0].hit_start
		blat_bam.write(read)
	bwa_bam.close()
	blat_bam.close()
	if hetero_factor != 'a':
		denovo.close()
	
	os.system('samtools sort -o ' + output + '.temp.sorted.bam ' + output + '.temp.bam')
	os.system('mv ' + output + '.temp.sorted.bam ' + output)
	os.system('samtools index ' + output)
	bwa_bam = pysam.Samfile(input, 'rb')
	# print input, bwa_bam
	readlen = 100
	empty = True
	for read in bwa_bam.fetch(until_eof=True):
		readlen = read.rlen
		empty = False
		break
	# print "readlen", readlen
	# readlen = bwa_bam.next().rlen
	bwa_bam.close()
	# os.system("cat "+ output+'.temp.fasta')
	return readlen, empty

def remove_assembly_fp(bam, input_vcf, output_vcf, len_cutoff, hetero_factor):
	"""remove false positives from assembly vcf file based on softclip reads enrichment"""
	input = vcf.Reader(open(input_vcf))
	output = vcf.Writer(open(output_vcf, 'w'), input)
	for record in input:
		if record.INFO['LEN'][0] < len_cutoff:
			continue
		else:
			chr = record.CHROM
			pos = record.POS
			if prob_of_indel_with_error(bam, chr, pos, hetero_factor) < 0.05:
				output.write_record(record)
	return

def usage():
	"""helping information"""
	print('Usage:')
	print(' python ScanIndel.py -p config.txt -i sample.txt [opts]')
	print('Opts:')
	print(' -o  :setting the output directory (default current working directory)')
	print(' -F  :setting min-alternate-fraction for FreeBayes (default 0.2)')
	print(' -C  :setting min-alternate-count for FreeBayes (default 2)')
	print(' -d  :setting min-coverage for Freebayes (default 0)')
	print(' -t  :setting --target for Freebayes to provide a BED file for analysis')
	print(' -s  :softclipping percentage triggering BLAT re-alignment (default 0.2)')

	print(' --min_percent_hq  :min percentage of high quality base in soft clipping reads, default 0.8')
	print(' --lowqual_cutoff  :low quality cutoff value, default 20')
	print(' --mapq_cutoff  :low mapping quality cutoff, default 1')
	print(' --blat_ident_pct_cutoff  :Blat sequence identity cutoff, default 0.8')
	print(' --gfServer_port  :gfServer service port number, default 50000')
	print(' --hetero_factor  :The factor about the indel\'s heterogenirity and heterozygosity, default 0.1')
	print(' --bam  :the input file is BAM format')
	print(' --rmdup  :exccute duplicate removal step before realignment')
	print(' --assembly_only  :only execute de novo assembly for indel calling without blat realignment (default False)')
	print(' --mapping_only  :only execute blat realignment withou de novo assembly for indel calling (default False)')
	print(' -h --help :produce this menu')
	print(' -v --version :show version of this tool')
	print('author: Rendong Yang <yang4414@umn.edu>, MSI, University of Minnesota, 2014')
	print('version: ' + __version__)

def external_tool_checking():
	"""checking dependencies are installed"""
	software = ['bedtools', 'bwa', 'samtools', 'gfServer', 'gfClient', 'freebayes', 'inchworm']
	cmd = "which"
	for each in software:
		try:
			path = subprocess.check_output([cmd, each], stderr=subprocess.STDOUT)
		except subprocess.CalledProcessError:
			print("Checking for '" + each + "': ERROR - could not find '" + each + "'", file=sys.stderr)
			print("Exiting.", file=sys.stderr)
			sys.exit(0)
		print("Checking for '" + each + "': found " + str(path))

def main():

	# parameters parsing
	cwd = os.getcwd()
	output_dir = cwd
	sample_file = ''
	config_file = ''
	freebayes_F = 0.2
	freebayes_C = 2
	softclip_ratio = 0.2
	depth = 0
	bam = False
	rmdup = False
	bedfile = ''
	assembly = True
	mapping = True
	min_percent_hq = 0.8
	lowqual_cutoff = 20
	hetero_factor = 0.1
	mapq_cutoff = 1
	blat_ident_pct_cutoff = 0.8
	gfServer_port = '50000'
	try:
		opts, args = getopt.getopt(sys.argv[1:], 'i:o:p:F:C:s:d:t:h:v', ['min_percent_hq=', 'lowqual_cutoff=', 'hetero_factor=', 'mapq_cutoff=', 'blat_ident_pct_cutoff=', 'gfServer_port=', 'bam', 'rmdup', 'assembly_only', 'mapping_only', 'help', 'version'])
		if not opts:
			print("Please use the -h or --help option to get usage information.")
			sys.exit(0)
	except getopt.GetoptError as err:
		print(str(err))
		usage()
		sys.exit(2)
	for o, a in opts:
		if o == '-i': sample_file = a
		elif o == '-o': 
			output_dir = a
			if not os.path.exists(output_dir):
				os.makedirs(output_dir)
		elif o == '-p': config_file = a
		elif o == '-F': freebayes_F = float(a)
		elif o == '-C': freebayes_C = int(a)
		elif o == '-s': softclip_ratio = float(a)
		elif o == '-d':	depth = int(a)
		elif o == '-t': bedfile = a
		elif o == '--min_percent_hq': min_percent_hq = float(a)
		elif o == '--lowqual_cutoff': lowqual_cutoff = int(a)
		elif o == '--hetero_factor': hetero_factor = float(a)
		elif o == '--mapq_cutoff': mapq_cutoff = int(a)
		elif o == '--blat_ident_pct_cutoff': blat_ident_pct_cutoff = float(a)
		elif o == '--gfServer_port': gfServer_port = a
		elif o == '--bam': bam = True
		elif o == '--rmdup': rmdup = True
		elif o == '--mapping_only': assembly = False
		elif o == '--assembly_only': mapping = False
		elif o in ('-v', '--version'):
			print('ScanIndel ' + __version__)
			sys.exit(0)
		elif o in ('-h', '--help'):
			usage()
			sys.exit(0)
		else:
			assert False, "unhandled option"
	if not assembly and not mapping:
		print('Cannot set both --mapping_only and --assembly_only.')
		usage()
		sys.exit(1)
	if not sample_file or not config_file:
		print('Please specify sample and config file.')
		usage()
		sys.exit(1)

	print('ScanIndel starts running: ' + time.strftime("%Y-%m-%d %H:%M:%S"))
	start = time.time()
	sample = read_sample_file(sample_file)
	reference = read_config_file(config_file)
	print('Loaded config and sample list')
	
	#check external tools used
	external_tool_checking()
	
	# start up BLAT server
	start_gfServer(reference['blat'], gfServer_port, output_dir)

	for each in sample:
		print("Analyzing sample:", each)
		if not bam:
			bwa_start = time.time()
			bwa_alignment(reference, sample[each], output_dir+'/'+each + '.temp.bwa.bam')
			blat_input = output_dir+'/'+each + '.temp.bwa.bam'
			bwa_end = time.time()
			print('BWA [ScanIndel] takes ' + str(bwa_end - bwa_start)+' seconds.')
		else:
			blat_input = sample[each]
			try:
				subprocess.check_call('samtools index ' + blat_input, shell=True)
			except subprocess.CalledProcessError as e:
				print("Execution failed for samtools index:", e, file=sys.stderr)
				sys.exit(1)
		if bedfile:
			try:
				subprocess.check_call("intersectBed -abam " + blat_input + " -b " + bedfile+  " > " + output_dir + "/" + each + ".temp.target.map.bam", shell=True)
			except subprocess.CalledProcessError as e:
				print("Execution failed for bedtools intersect:", e, file=sys.stderr)
				sys.exit(1)
			os.system("samtools view -b -f 4 " + blat_input + " >" + output_dir + "/" + each + ".temp.target.unmap.bam")
			os.system("samtools merge -f " + output_dir + "/" + each + ".temp.target.bam " + output_dir + "/" + each + ".temp.target.map.bam " + output_dir + "/" + each + ".temp.target.unmap.bam")
			os.system("samtools index " + output_dir + "/" + each + ".temp.target.bam")
			blat_input = output_dir + "/" + each + '.temp.target.bam'
		if rmdup:
			try:
				subprocess.check_call("samtools rmdup " + blat_input + " " + output_dir + "/" + each + ".temp.rmdup.bam", shell=True)
			except subprocess.CalledProcessError as e:
				print("Execution failed for samtools rmdup:", e, file=sys.stderr)
				sys.exit(1)
			os.system("samtools index " + output_dir + "/" + each + ".temp.rmdup.bam")
			blat_input = output_dir + "/" + each + '.temp.rmdup.bam'

		# extracting candidate soft-clipped reads for realignment and/or assembly
		blat_start = time.time()
		readlen, empty = blat_alignment(mapping, reference, softclip_ratio, lowqual_cutoff, min_percent_hq, mapq_cutoff, blat_ident_pct_cutoff, gfServer_port, hetero_factor, blat_input, output_dir + "/" + each + '.reads.bam')
		blat_end = time.time()
		print('BLAT [ScanIndel] takes ' + str(blat_end - blat_start) + ' seconds.')
		
		if assembly and os.path.getsize(output_dir + "/" + each + '.reads.bam.temp.fasta'):
			print("start assemby")
			assembly_start = time.time()
			# de novo assemble softclip reads with breakpoint evidence and unmapped reads"
			try:
				subprocess.check_call('inchworm --reads ' + output_dir + "/" + each + '.reads.bam.temp.fasta --run_inchworm --DS -L ' + str(readlen + 1) + ' >' + output_dir + "/" + each + '.temp.contig', shell=True)
				print(('inchworm --reads ' + output_dir + "/" + each + '.reads.bam.temp.fasta --run_inchworm --DS -L ' + str(readlen + 1) + ' >' + output_dir + "/" + each + '.temp.contig'))
			except subprocess.CalledProcessError as e:
				print("Execution failed for inchworm:", e, file=sys.stderr)
				sys.exit(1)
			bwa_alignment(reference, output_dir + '/' + each + '.temp.contig', output_dir + '/' + each + '.denovo.temp.bwa.bam')
			blat_alignment(mapping, reference, softclip_ratio, lowqual_cutoff, min_percent_hq, mapq_cutoff, blat_ident_pct_cutoff, gfServer_port, 'a', output_dir+'/'+each+'.denovo.temp.bwa.bam', output_dir+'/'+each + '.contigs.bam')
			assembly_end = time.time()
			print('Assembly [ScanIndel] takes ' + str(assembly_end - assembly_start) + ' seconds.')
		else:
			print ("We have not found softclip and unmapped reads.")
			sys.exit(1)

		freebayes_start = time.time()
		try:
			subprocess.check_call('freebayes -I -X -u -F ' + str(freebayes_F) + ' -C ' + str(freebayes_C) + ' --min-coverage ' + str(depth) + ' -f ' + reference['freebayes'] + ' ' + output_dir + '/' + each + '.reads.bam > ' + output_dir + '/' + each + '.mapping.indel.vcf', shell=True)
		except subprocess.CalledProcessError as e:
			print("Execution failed for freebayes:", e, file=sys.stderr)
			sys.exit(1)
		if assembly:
			try:
				subprocess.check_call('freebayes -I -X -u -F 0 -C 1 -f ' + reference['freebayes'] + ' ' + output_dir + '/' + each + '.contigs.bam > ' + output_dir + '/' + each + '.temp.indel.vcf', shell=True)
				# print 'freebayes -I -X -u -F 0 -C 1 -f ' + reference['freebayes'] + ' ' + output_dir + '/' + each + '.contigs.bam > ' + output_dir + '/' + each + '.temp.indel.vcf'
			except subprocess.CalledProcessError as e:
				print("Execution failed for freebayes in assembly:", e, file=sys.stderr)
				sys.exit(1)
			remove_assembly_fp(blat_input, output_dir + '/' + each + '.temp.indel.vcf', output_dir + '/' + each + '.assembly.indel.vcf', readlen * softclip_ratio, hetero_factor)
			path = os.path.dirname(os.path.realpath(__file__))
			os.system('python3 ' + path + '/vcf-combine.py ' + output_dir + '/'+ each + '.mapping.indel.vcf ' + output_dir + '/' + each + '.assembly.indel.vcf | bedtools sort -i stdin -header > ' + output_dir + '/' + each + '.merged.indel.vcf')
		freebayes_end = time.time()
		print('Freebayes [ScanIndel] takes ' + str(freebayes_end - freebayes_start) + ' seconds.')

	os.system('gfServer stop localhost '+gfServer_port)
	os.system('rm '+output_dir+'/*.temp.*')

	print("ScanIndel running done: " + time.strftime("%Y-%m-%d %H:%M:%S"))
	end = time.time()
	print('ScanIndel takes ' + str(end - start) + ' seconds.')

if __name__ == '__main__':
	main()

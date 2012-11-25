#!/software/bin/python

## Copyright alexandros.kanterakis@gmail.com

import operator


def Load_PEDMAP_user_Kantale(
	map_filename = None,
	ped_filename = None,
	key = "reference",  # key can be either "reference" or "position"
	verbose = True,
	):
	if verbose: print "Loading map: ", map_filename, " ..."
	snps= {}
	if key == "reference":
		key_index = 1
		value_index = 3
	elif key == "position":
		key_index = 3
		value_index = 1
	else:
		raise Exception("Unknown key parameter: " + str(key))
	
	snp_names = []
	map_file = open(map_filename, 'U')
	while True:
		line = map_file.readline()
		if not line:
			break
		
		line_s = line.replace("\n", "").split()
		snps[line_s[key_index]] = line_s[value_index]
		snp_names += [line_s[value_index]]
	map_file.close()
#	snp_names = [y[key_index] for y in [[operator.setitem(snps, x.split()[1], []), x.split()[1]] for x in open(map_filename, 'U').readlines()]]
	if verbose: print "..Done"
	
	if verbose: print "Loading PED: ", ped_filename, " ..."
	ped_file = open(ped_filename, 'U')
	index_lines = 0
	sample_names = []
	while True:
		line = ped_file.readline()
		if not line:
			break
		
		index_lines += 1
		
		line_s = line.replace("\n", "").split()
		sample_name = line_s[1]
		sample_names += [sample_name]
		genotypes = line_s[6:]
		
		if verbose:
			print "Reading ped sample:", index_lines, sample_name
		
		for z in xrange(0, len(genotypes), 2):
				if not snps.has_key(snp_names[z/2]):
					snps[snp_names[z/2]] = []

				snps[snp_names[z/2]] += [(genotypes[z], genotypes[z+1])]
		
	ped_file.close()
	if verbose: print "...Done"
	# The previous in one line:
	#_=[[ operator.iadd(snps[snp_names[z/2]], [(y[z], y[z+1])]) for z in xrange(0,len(y),2)] for y in [ x.split()[6:] for x in open(ped_filename, 'U').readlines()]]
	
	return snps, sample_names


def Convert_PEDMAP_to_VCF_user_Kantale(
	input_ped_filename = None,
	input_map_filename = None,
	output_filename = None,
	verbose = False,
	):

	output = open(output_filename, "w")
	if verbose:
		print "Loading PED / MAP files.."
	snps, sample_names = Load_PEDMAP_user_Kantale(map_filename=input_map_filename, ped_filename=input_ped_filename)
	if verbose:
		print "..DONE"

	header, header_str = Default_VCF_header_user_Kantale(header = None, sample_names = sample_names)

	output.write(header_str)
	filtered = 0

	for line in open(input_map_filename).readlines():
		line_splitted = line.replace("\n", "").split()
		
		to_print = [line_splitted[0], line_splitted[3], line_splitted[1]]

		al1, al2 = (None, None)
		observations_chr1 = []
		observations_chr2 = []
		al1 = None
		for gen in snps[line_splitted[1]]:
			if not al1 and gen[0] != "0": al1 = gen[0]
			if not al1 and gen[1] != "0": al1 = gen[1]

			if not al1: continue

			if al1 != gen[0] and gen[0] != "0": al2 = gen[0]
			if al1 != gen[1] and gen[1] != "0": al2 = gen[1]

			if gen[0] != "0": observations_chr1 += [gen[0]]
			if gen[1] != "0": observations_chr2 += [gen[1]]

		if al1 and al2: 
			(maf, minor_allele, obs) = Minor_allele_frequency_user_Kantale(
						allele1 = al1,
						allele2 = al2,
						observations_chr1 = observations_chr1,
						observations_chr2 = observations_chr2
									)
			alternative = minor_allele
			reference = al1 if minor_allele == al2 else al2
		else:
			if al1:
				reference = al1
				alternative = "."
			elif al2:
				reference = al2
				alternative = "."
			else:
				if verbose:
					print "warning: Cannot infer reference/alternatice information to SNP:", line_splitted[1]
					print snps[line_splitted[1]]
				filtered += 1
				continue
	 
		to_print += [ reference, alternative, "100", "PASS", "NS=" + str(len(snps[line_splitted[1]])), "GT"]
		for gen in snps[line_splitted[1]]:
			if gen[0] == reference: genotype = "0"
			elif gen[0] == alternative: genotype = "1"
			elif gen[0] == "0": genotype = "."
			else:
				print gen
				print "reference:", reference
				print "alternative:", alternative
				raise Exception("exmm..")
			if gen[1] == reference: genotype += "/0"
			elif gen[1] == alternative: genotype += "/1"
			elif gen[1] == "0": genotype += "/."
			else:
				print "genotype:", gen
				print "reference:", reference
				print "alternative:", alternative
				raise Exception("exmm..")
		
			to_print += [genotype]

		output.write(str.join("\t", to_print) + "\n")


	output.close()
	print "Filtered SNPs:", filtered, "Cannot infer reference/alternatice information"



if __name__ == '__main__':
    Convert_PEDMAP_to_VCF_user_Kantale(
	input_ped_filename = 'omni2.5-8_20120516_gwa_ugand_gtu.ped',
	input_map_filename = 'omni2.5-8_20120516_gwa_ugand_gtu.map',
	output_filename = 'omni2.5-8_20120516_gwa_ugand_gtu.vcf',
	verbose = True,
        )

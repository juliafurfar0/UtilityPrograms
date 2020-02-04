localrules: SelectSmallestFeature,concatenate,movingMaskedCrms,movingMaskedNegs
species=['Dana', 'Dere','Dgri', 'Dmoj', 'Dper', 'Dpse', 'Dsec', 'Dsim', 'Dvir', 'Dyak','Dmel']
inputFile="testDm6"
tsetName="mapping1.test"
# rule all:
# 	input:
# 		"out/testDm6_15_output.bed"
# 		"out/testDm6_15_output_sorted.bed"
		# "out/testDm6_15_output.bed"
		# "out/testDm6_15_output_sorted.bed"

rule SelectSmallestFeature:
	input:
		expand("in/{input}.bed",input=inputFile)
	output:
		expand("{input}_output.bed",input=inputFile)
	conda:
		"envs/mapping.yml"
	shell:
		"scripts/SelectSmallestFeature.py -i {input} -o {output}"
# rule sort:
# 	input:
# 		expand("out/{input}_output.bed",input=inputFile)
# 	output:
# 		expand("{input}_output_sorted.bed",input=inputFile)
# 	conda:
# 		"envs/mapping.yml"
# 	shell:
# 		"bedtools sort -i {input} > {output}"

rule liftOver:
	input:
		bed_file= expand("{input}_output.bed",input=inputFile),
		chain_files_dir= "/gpfs/scratch/mshalfon/Training_Set_Construction_Pipeline/new_chain_files/",
		genome_files_dir="/gpfs/scratch/mshalfon/Training_Set_Construction_Pipeline/new_fasta_files/",
		path_to_executable="/gpfs/scratch/mshalfon/Training_Set_Construction_Pipeline/"
	output:
	   	expand("{specie}.{input}_output.bed.fa", specie=species, input=inputFile)
	shell:
		#"cp out/testDm6_15_output_sorted.bed ."
		"scripts/liftover_script_modified.sh {input.bed_file} {input.path_to_executable} {input.chain_files_dir} {input.genome_files_dir}"

rule concatenate:
	input:
		out=expand("{specie}.{input}_output.bed.fa", specie=species, input=inputFile)
		#tsetName=expand("{tset}",tset=tsetName)
	output:	
 		fileName=expand("{tset}/crms.fa",tset=tsetName)
		#tsetName=expand("{tset}",tset=tsetName)
	shell:
		"""
		cat {input.out} > {output.fileName}
		mkdir fasta_log_files
		mv D* fasta_log_files/
		#$tsetName
		#echo $tsetName
		#mkdir {tsetName}
		#mv {output.fileName} {tsetName}

		"""

rule negative:
	input:
		infileName=expand("{tset}/crms.fa",tset=tsetName),
		genomedir="/projects/academic/mshalfon/RandomSameGC_files/finalchr",
		gene="/projects/academic/mshalfon/RandomSameGC_files/genes.goodchronly"

	output:
		outfileName=expand("{tset}/neg.fa",tset=tsetName)
	shell:
		"perl scripts/randomWithSameGC.pl --crm {input.infileName} --output {output.outfileName} --size 10 --genomedir {input.genomedir} --gene {input.gene} --suffix .fa "

rule maskingCrms:
	input:
		inCrms=expand("{tset}/crms.fa",tset=tsetName)

	output:
		#outCrms=expand("{tset}/crms.fasta",tset=tsetName)
		outCrms="crms.fa.2.7.7.80.10.50.500.mask"
	shell:
		"scripts/trf409.linux64 {input.inCrms} 2 7 7 80 10 50 500 -m -h || true"

rule movingMaskedCrms:
	input:
		maskedCrms="crms.fa.2.7.7.80.10.50.500.mask",
		unMaskedCrms=expand("{tset}/crms.fa",tset=tsetName)
	output:
		maskedCrmsName=expand("{tset}/crms.fasta",tset=tsetName)
	shell:
		"""
		mv {input.maskedCrms} {output.maskedCrmsName}
		rm {input.unMaskedCrms}
		rm crms.fa.2.7.7.80.10.50.500.dat 
		"""

rule masking_negs:
	input:
		inNegs=expand("{tset}/neg.fa",tset=tsetName)
	output:
		outNegs="neg.fa.2.7.7.80.10.50.500.mask"
	shell:
		"scripts/trf409.linux64 {input.inNegs} 2 7 7 80 10 50 500 -m -h || true"

rule movingMaskedNegs:
	input:
		maskedNegs="neg.fa.2.7.7.80.10.50.500.mask",
		unMaskedNegs=expand("{tset}/neg.fa",tset=tsetName)
	output:
		maskedNegsName=expand("{tset}/neg.fasta",tset=tsetName)

	shell:
		"""
		mv {input.maskedNegs} {output.maskedNegsName}
		rm {input.unMaskedNegs}
		rm neg.fa.2.7.7.80.10.50.500.dat 
		"""

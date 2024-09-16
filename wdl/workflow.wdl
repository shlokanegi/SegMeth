version 1.0

workflow SegMeth {
    meta {
		author: "Shloka Negi"
		email: "shnegi@ucsc.edu"
		description: "Segmentation of methylome using circular binary segmentation. By default, runs whole-genome."
		}

	parameter_meta {
		SAMPLE: "Sample name"
		TARGET_BED: "Presegmented Genome BED file or any other targets BED file (e.g. CpG Islands, CCREs, VISTA Enhancers, etc..)"
		MODKIT_BEDS: "Haplotype-specific modkit CpG BED files"
		REGION_TYPE: "Descriptor of the regions of interest in the bed. eg: promoters, CpG_Islands, CCREs, etc.."
	}

	input {
		Array[String] CHRS = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
		String SAMPLE
		File TARGET_BED
		String REGION_TYPE
		Array[File] MODKIT_BEDS
	}

	scatter (chr in CHRS) {
		call runsegmeth {
			input:
			chrom=chr,
			sample=SAMPLE,
			target_bed=TARGET_BED,
			region_type=REGION_TYPE,
			modkit_beds=MODKIT_BEDS
		}
	}

	call concatenateBed {
		input:
		bed1files=flatten(runsegmeth.BED1s),
		bed2files=flatten(runsegmeth.BED2s),
		sample=SAMPLE,
		region_type=REGION_TYPE
	}

	output {
		File outBed1 = concatenateBed.concatenatedBed1
		File outBed2 = concatenateBed.concatenatedBed2
	}
}

task runsegmeth {
	input {
		String chrom
		String sample
		Array[File] modkit_beds
		File target_bed
		String region_type
		Int memSizeGB = 64
	}

	Int diskSizeGB = round(5*(size(target_bed, "GB") + size(modkit_beds, "GB"))) + 50
	String suffix = "~{region_type}"+"_segments_offsets.bed"

	command <<<
		# Set the exit code of a pipeline to that of the rightmost command
		# to exit with a non-zero status, or zero if all commands of the pipeline exit
		set -o pipefail
		# cause a bash script to exit immediately when a command fails
		set -e
		# cause the bash shell to treat unset variables as an error and exit immediately
		set -u
		# echo each line of the script to stdout so we can see what is happening
		# to turn off echo do 'set +o xtrace'
		set -o xtrace

		# initilization  
		#rm -rf "~{sample}.hp1.~{suffix}" "~{sample}.hp2.~{suffix}"

		# Extract targets in "chr"
		echo "Running first command"
		cat ~{target_bed} | awk -v c=~{chrom} '$1==c {print $0}' | sort -k1,1 -k2,2n -k3,3n  > "~{sample}.target.chr.bed"

		# Extract haplotype-specific BED files and run CBS per chromosome
		FID=1
		for MODKIT_FULL_PATH in ~{sep=" " modkit_beds}
		do
			MODKIT_FILENAME=$(basename -- "$MODKIT_FULL_PATH")
			if [[ $MODKIT_FILENAME =~ \.gz$ ]]; then
				gunzip -c $MODKIT_FULL_PATH | sort -k1,1 -k2,2n -k3,3n >> ~{sample}.hp${FID}.bed
			else
				cat $MODKIT_FULL_PATH | sort -k1,1 -k2,2n -k3,3n >> ~{sample}.hp${FID}.bed
			fi
			## Extract CpG sites in "chr"
			cat ~{sample}.hp${FID}.bed | awk -v c=~{chrom} '$1==c {print $0}' > ~{sample}.hp${FID}.chr.bed
			## Run CBS
			python3 /opt/scripts/SegMeth-v1.0/cbs_hp.py -file ~{sample}.hp${FID}.chr.bed -o ~{sample}.~{chrom}.hp${FID}.~{suffix} -p 1E-3 -s ~{sample}.hp${FID} -t ~{sample}.target.chr.bed -plot no -merge yes

			FID=$((FID+1))
		
		done
	>>>

	output {
		Array[File] BED1s = glob("*.hp1.~{suffix}")
		Array[File] BED2s = glob("*.hp2.~{suffix}")
	}
	runtime {
		memory: memSizeGB + " GB"
		disks: "local-disk " + diskSizeGB + " SSD"
		docker: "quay.io/shnegi/segmeth@sha256:9a74121bf8059af720504cdcb74e7ee2248efd4e488af5286a1438971148051e"
		preemptible: 1
	}
}

task concatenateBed {
	input {
		Array[File] bed1files
		Array[File] bed2files
		String sample
		String region_type
		Int memSizeGB = 64
		Int diskSizeGB = 250
	}

	command {
		cat ${sep=" " bed1files} | sort -k1,1 -k2,2n -k3,3n > ~{sample}.~{region_type}.hp1.bed
		cat ${sep=" " bed2files} | sort -k1,1 -k2,2n -k3,3n > ~{sample}.~{region_type}.hp2.bed
	}

	output {
		File concatenatedBed1 = "~{sample}.~{region_type}.segments.hp1.bed"
		File concatenatedBed2 = "~{sample}.~{region_type}.segments.hp2.bed"
	}

	runtime {
		memory: memSizeGB + " GB"
		disks: "local-disk " + diskSizeGB + " SSD"
		docker: "quay.io/biocontainers/bcftools@sha256:f3a74a67de12dc22094e299fbb3bcd172eb81cc6d3e25f4b13762e8f9a9e80aa"
		preemptible : 1
	}
}
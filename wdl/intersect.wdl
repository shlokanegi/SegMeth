version 1.0

workflow SegMethIntersect {
    meta {
		author: "Shloka Negi"
		email: "shnegi@ucsc.edu"
		description: "Generates an intersection BED file that combines the SegMeth BED files for multiple samples"
		}

	parameter_meta {
		BEDFILES: "SegMeth output BED files"
		COHORT: "Name of cohort. (eg. HBCC, NABEC, etc.)"
		REGION_TYPE: "Descriptor of the regions of interest in the bed. eg: promoters, CpG_Islands, CCREs, etc.."
	}

	input {
			Array[File] BEDFILES   # List of BED input files
			String COHORT
			String REGION_TYPE
		}

		call intersect {
			input:
				bedfiles=BEDFILES,
				cohort=COHORT,
				region_type=REGION_TYPE
		}

		output {
			File intersect_bed = intersect.bed
		}
}


task intersect {

	input {
		Array[File] bedfiles
		String cohort
		String region_type
		Int memSizeGB = 64
	}

	Int diskSizeGB = round(3*(size(bedfiles, "GB"))) + 20
	String outfile = "~{cohort}.~{region_type}"+".segmeth_intersection.bed"

	command {

		python3 /opt/scripts/SegMeth-v1.0/intersection.py -i ${sep=" " bedfiles} -o ~{outfile}
	}

	output {
		File bed = "~{outfile}"
	}

	runtime {
		memory: memSizeGB + " GB"
		disks: "local-disk " + diskSizeGB + " SSD"
		docker: "quay.io/shnegi/segmeth:latest"
		preemptible: 1
	}
}
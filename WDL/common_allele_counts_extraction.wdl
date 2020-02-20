
workflow common_allele_counts_extraction {

    meta {
        description: "extracts allelic counts from bam file after first marking duplicates."
    }

    ## file containing the list of bam files with gs:// paths
    File input_bams_listing 
    

    ## common alleles / dbsnp
    File input_vcf
    File input_vcf_tbi

    ## ref genome info
    File ref_genome_fa
    File ref_genome_fai
    File ref_genome_dict

    Array[Array[String]] bam_files = read_tsv(input_bams_listing)

    call get_bam_bai_files as get_bams {
        input:
          input_bams_listing = input_bams_listing
    }
        
    scatter(bam_bai_files in get_bams.bam_bai_files) {

	    call mark_duplicates as md {
	        input:
	          bam_file = bam_bai_files[0],
	          bai_file = bam_bai_files[1]
	    }
        

	    call allelic_counts as ac {
	       input:
	          bam_file = md.bam_file_dups_marked,
              bai_file = md.bai_file_dups_marked,
	          vcf_file = input_vcf,
	          vcf_idx_file = input_vcf_tbi,
              ref_genome_fa = ref_genome_fa,
              ref_genome_fai = ref_genome_fai,
              ref_genome_dict = ref_genome_dict
	    }

	    call extract_covered_sites as ecs {
	        input:
	            allelic_counts_file = ac.allelic_counts_file
	    
	    }
    }
    

    call summarize_allelic_counts_outputs {
	    input:
            files = ecs.allelic_counts_covered_sites
    }


    ## pipeline output
    output {
	    File allelic_count_file_locations_tsv = summarize_allelic_counts_outputs.file_locations_tsv
    }
}




task mark_duplicates {

    File bam_file
    File bai_file

    # runtime values
    String docker = "quay.io/humancellatlas/secondary-analysis-picard:v0.2.2-2.10.10"
    Int machine_mem_mb = 8250
    # give the command 1 GiB of overhead
    Int command_mem_mb = machine_mem_mb - 1000
    Int cpu = 2
    # use provided disk number or dynamically size on our own, with 10GiB of additional disk
    Int disk = ceil(size(bam_file, "GiB") + 10)
    Int preemptible = 5
    
    String dups_marked_bam_file = basename(bam_file) + ".duplicates_marked.bam"
    String dups_marked_bai_file = "${dups_marked_bam_file}.bai"
    
    
    command <<<

       set -e
	
       java -Xmx${command_mem_mb}m -XX:ParallelGCThreads=${cpu}  -jar /usr/picard/picard.jar  MarkDuplicates \
          VALIDATION_STRINGENCY=SILENT  \
          INPUT=${bam_file} \
          OUTPUT="${dups_marked_bam_file}" \
          METRICS_FILE="${dups_marked_bam_file}.dups_marked.metrics" \
          ASSUME_SORTED=true \
          REMOVE_DUPLICATES=false


       java -Xmx${command_mem_mb}m -XX:ParallelGCThreads=${cpu}  -jar /usr/picard/picard.jar BuildBamIndex \
     	  I="${dups_marked_bam_file}" \
	      O="${dups_marked_bai_file}"
  
	
     >>>
  
     runtime {
        docker: docker
        memory: "${machine_mem_mb} MiB"
        disks: "local-disk ${disk} HDD"
        cpu: cpu
        preemptible: preemptible
     }
     
     output {
         File bam_file_dups_marked = "${dups_marked_bam_file}"
	     File bai_file_dups_marked = "${dups_marked_bai_file}"
     }
     
     
}

task allelic_counts {

    File bam_file
    File bai_file

    File vcf_file
    File vcf_idx_file

    File ref_genome_fa
    File ref_genome_fai
    File ref_genome_dict
    
    # runtime values
    String docker = "broadinstitute/gatk:4.1.4.1"
    Int machine_mem_mb = 8250
    # give the command 1 GiB of overhead
    Int command_mem_mb = machine_mem_mb - 1000
    Int cpu = 2
    # use provided disk number or dynamically size on our own, with 10GiB of additional disk
    Int disk = ceil(size(bam_file, "GiB") + size(vcf_file, "GiB") + 10)
    Int preemptible = 5
    
    String allelic_counts_filename =  basename(bam_file) + ".allelic_counts.tsv"
        
    command <<<

       set -e
	
        java -Xmx${command_mem_mb}m -XX:ParallelGCThreads=${cpu} -jar /gatk/gatk.jar \
             CollectAllelicCounts -I ${bam_file} \
                                  -R ${ref_genome_fa} \
                                  -L ${vcf_file} \
                                  -O ${allelic_counts_filename}
        
	
     >>>
  
     runtime {
        docker: docker
        memory: "${machine_mem_mb} MiB"
        disks: "local-disk ${disk} HDD"
        cpu: cpu
        preemptible: preemptible
     }
     
     output {
         File allelic_counts_file = "${allelic_counts_filename}"
     }
     
     
}


task extract_covered_sites {

    File allelic_counts_file

    String allelic_counts_covered_sites_filename = basename(allelic_counts_file) + ".covered_sites.tsv"

    Int machine_mem_mb = 8250
    # give the command 1 GiB of overhead
    Int command_mem_mb = machine_mem_mb - 1000
    Int cpu = 2
    # use provided disk number or dynamically size on our own, with 10GiB of additional disk
    Int disk = ceil(size(allelic_counts_file, "GiB") + 10)
    Int preemptible = 5

    
    command <<<

        set -eou pipefail
        
        cat ${allelic_counts_file} | egrep -v ^\@ | perl -lane ' if ($F[2] || $F[3]) { print; } ' > ${allelic_counts_covered_sites_filename}

    >>>
    
    runtime {
        memory: "${machine_mem_mb} MiB"
        disks: "local-disk ${disk} HDD"
        cpu: cpu
        preemptible: preemptible
    }    

    output {
        File allelic_counts_covered_sites = "${allelic_counts_covered_sites_filename}"
    }
    
}

task summarize_allelic_counts_outputs {

    Array[String] files

    String dollar = "$"

    String file_locations_tsv_filename = "allelic_counts_files_listing.tsv"

    command <<< 

        set -e

        for file in ${sep=" " files}; do
            echo ${dollar}{file} >> ${file_locations_tsv_filename}
        done

    >>>

    runtime {
        memory: "1000 MiB"
        disks: "local-disk 10 HDD"
        cpu: 1
        preemptible: 5
    }    
        
    output {
        File file_locations_tsv = "${file_locations_tsv_filename}"
    }
    
}


task get_bam_bai_files {

    File input_bams_listing

    command <<<

        set -eou pipefail

        cat ${input_bams_listing} | perl -lane ' print "$F[0]\t$F[0].bai"; '

    >>>

    output {
        Array[Array[String]] bam_bai_files = read_tsv(stdout())
    }

}

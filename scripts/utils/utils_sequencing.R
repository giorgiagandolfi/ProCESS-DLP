simulate_seq_resources <- function(c,tumour,ref_path,coverage,purity,singularity_version,output_local_dir){
  p_info <- ps::ps_handle()
  start_time <- Sys.time()
  initial_cpu <- ps::ps_cpu_times(p_info)
  initial_mem <- ps::ps_memory_info(p_info)["rss"] / 1024^3
  if (tumour){
    if (singularity_version=="old"){
      seq_res <- simulate_seq(phylo_forest, reference_genome = ref_path,
                              chromosomes = c,
                              coverage = coverage,
                              purity = purity, 
                              write_SAM = TRUE, 
                              read_size = 150,
                              sequencer = basic_seq,
                              insert_size_mean = 350,
                              insert_size_stddev = 10,
                              output_dir = output_local_dir,
                              include_non_sequenced_mutations = F,
                              update_SAM = TRUE,
                              with_normal_sample = FALSE)  
    } else {
      seq_res <- simulate_seq(
        phylo_forest,
        coverage = coverage,
        purity = purity, 
        write_SAM = TRUE,
        with_normal_sample = FALSE,
        chromosomes = c,
        read_size = 150,
        sequencer = no_error_seq,
        insert_size_mean = 350,
        insert_size_stddev = 10,
        missed_SID_statistics=F, germline_statistics=F,
        wide_format=FALSE,
        output_dir = output_local_dir,
        update_SAM = TRUE
      )
    }
  } else{
    cat("simulate_normal_seq")
    seq_res <- simulate_normal_seq(phylo_forest, reference_genome = ref_path,
                                   chromosomes = c,
                                   coverage = coverage,
                                   write_SAM = TRUE, 
                                   read_size = 150,
                                   sequencer = basic_seq,
                                   insert_size_mean = 350,
                                   include_non_sequenced_mutations = TRUE,
                                   insert_size_stddev = 10,
                                   filename_prefix = sam_filename_prefix,
                                   template_name_prefix = paste0(lot_name,'r'),
                                   output_dir = output_local_dir,
                                   with_preneoplastic = with_preneoplastic,
                                   update_SAM = TRUE)
  }
  
  end_time <- Sys.time()
  final_cpu <- ps::ps_cpu_times(p_info)
  final_mem <- ps::ps_memory_info(p_info)["rss"] / 1024^3
  elapsed_time <- end_time - start_time
  elapsed_time <- as.numeric(elapsed_time, units = "mins")
  cpu_used <- (final_cpu["user"] + final_cpu["system"]) - (initial_cpu["user"]+ initial_cpu["system"])
  mem_used <- final_mem - initial_mem
  resource_usage <- data.frame(
    elapsed_time_mins =  elapsed_time,
    cpu_time_secs = cpu_used,
    memory_used_MB = mem_used,
    chr = c,
    coverage = coverage,
    purity = purity,
    tumour = tumour
  )
  rownames(resource_usage) <- NULL
  out <- list(mutations=seq_res$mutations,parameters=seq_res$parameters,
              resource_usage = resource_usage)
  return(out)
}
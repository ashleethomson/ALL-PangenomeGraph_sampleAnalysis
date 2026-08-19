# Load required packages
library(GenomicRanges)
library(data.table)
library(tidyverse)


sample_ids <- c(
  "A2189-1026",
  "A2241-10841",
  "A2638-12600",
  "A4551-18835",                  
  "A5058-21220",
  "A5322-22937",
  "A5560-24589",
  "A5654-25234",
  "A6314-28788",
  "ADU-0705-201321",
  "ADU-0919-201785",
  "ADU-1262-202602",
  "ADU-1281-202610",
  "AYA-0846-201738",
  "AYA-0892-201728",
  "AYAII-0454-DIA1-BM",
  "CHI-0796-202173",
  "CHI-1383-202781")

# Define buffer size (bp on each side of homopolymer region)
buffer <- 5

# Read in the BED file of homopolymer regions (only needs to load once)
bed <- fread("Homopolymer/GRCh38_homopolymerRegions_6bp.bed", header = FALSE)
colnames(bed) <- c("chr", "start", "end", "base", "length")

# Helper to standardize chromosome names to "chr" prefix style
normalize_chr <- function(x) {
  ifelse(grepl("^chr", x), x, paste0("chr", x))
}

bed$chr <- normalize_chr(bed$chr)

# Create GRanges for homopolymer regions (only needs to run once)
gr_homopolymer <- GRanges(
  seqnames = bed$chr,
  ranges = IRanges(start = bed$start, end = bed$end)
)

# Create buffered GRanges for homopolymer regions
gr_homopolymer_buffered <- gr_homopolymer + buffer

# Loop through each sample
for (sample_id in sample_ids) {
  
  cat(sprintf("\nProcessing sample: %s\n", sample_id))
  
  # Build file paths
  multianno_path <- sprintf("Data/Linear/Input/%s_OG_INDEL.hg38_multianno.txt", sample_id)
  filtered_path  <- sprintf("Data/Linear/Input/tsv/%s_OG_INDEL.hg38_filteredIndels.tsv", sample_id)
  # Skip if files don't exist
  if (!file.exists(multianno_path)) {
    warning(sprintf("Multianno file not found for sample %s, skipping.", sample_id))
    next
  }
  if (!file.exists(filtered_path)) {
    warning(sprintf("Filtered file not found for sample %s, skipping.", sample_id))
    next
  }
  
  # Read in files
  multianno <- fread(multianno_path, header = FALSE) %>% select(1:5)
  colnames(multianno)[1:5] <- c("chr", "start", "end", "ref", "alt")
  multianno$chr <- normalize_chr(multianno$chr)
  
  filtered <- fread(filtered_path, header = TRUE) %>% select(1:5)
  colnames(filtered)[1:5] <- c("chr", "start", "end", "ref", "alt")
  filtered$chr  <- normalize_chr(filtered$chr)
  
  
  
  # Create GRanges for variants
  gr_variants <- GRanges(
    seqnames = multianno$chr,
    ranges = IRanges(start = multianno$start, end = multianno$end)
  )
  gr_filtered_variants <- GRanges(
    seqnames = filtered$chr,
    ranges = IRanges(start = filtered$start, end = filtered$end)
  )
  
  # Find overlaps with exact and buffered homopolymer regions
  hits_exact             <- findOverlaps(gr_variants, gr_homopolymer)
  hits_buffered          <- findOverlaps(gr_variants, gr_homopolymer_buffered)
  hits_filtered_exact    <- findOverlaps(gr_filtered_variants, gr_homopolymer)
  hits_filtered_buffered <- findOverlaps(gr_filtered_variants, gr_homopolymer_buffered)
  
  # Subset overlapping variants (using buffered to capture within + adjacent)
  variants_in_homopolymer          <- multianno[queryHits(hits_buffered), ]
  filtered_variants_in_homopolymer <- filtered[queryHits(hits_filtered_buffered), ]
  
  # Add homopolymer details
  variants_in_homopolymer$homopolymer_base            <- bed$base[subjectHits(hits_buffered)]
  variants_in_homopolymer$homopolymer_length          <- bed$length[subjectHits(hits_buffered)]
  filtered_variants_in_homopolymer$homopolymer_base   <- bed$base[subjectHits(hits_filtered_buffered)]
  filtered_variants_in_homopolymer$homopolymer_length <- bed$length[subjectHits(hits_filtered_buffered)]
  
  # Label as "within" or "adjacent"
  variants_in_homopolymer$overlap_type <- ifelse(
    queryHits(hits_buffered) %in% queryHits(hits_exact),
    "within",
    "adjacent"
  )
  filtered_variants_in_homopolymer$overlap_type <- ifelse(
    queryHits(hits_filtered_buffered) %in% queryHits(hits_filtered_exact),
    "within",
    "adjacent"
  )
  
  # Write output (one file per sample)
  fwrite(variants_in_homopolymer,
         sprintf("Homopolymer/Linear_%s_variants_in_homopolymer.txt", sample_id), sep = "\t")
  fwrite(filtered_variants_in_homopolymer,
         sprintf("Homopolymer/Linear_%s_filtered_variants_in_homopolymer.txt", sample_id), sep = "\t")
  
  # Print summary
  cat(sprintf("  Multianno variants within homopolymer regions:   %d\n",
              sum(variants_in_homopolymer$overlap_type == "within")))
  cat(sprintf("  Multianno variants adjacent to homopolymer regions: %d\n",
              sum(variants_in_homopolymer$overlap_type == "adjacent")))
  cat(sprintf("  Filtered variants within homopolymer regions:   %d\n",
              sum(filtered_variants_in_homopolymer$overlap_type == "within")))
  cat(sprintf("  Filtered variants adjacent to homopolymer regions: %d\n",
              sum(filtered_variants_in_homopolymer$overlap_type == "adjacent")))
}

cat("\nDone!\n")
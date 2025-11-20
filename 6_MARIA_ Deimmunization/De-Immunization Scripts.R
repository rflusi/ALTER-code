library(data.table)
library(dplyr)
library(parallel)
library(progress)

#General Notes:
#For all residue index index starts at 1, position indices are inclusive

#------------------------------------------------------------------------------
# Purpose:
#   Identifies all peptide windows of a specified length range within a protein
#   sequence and returns only those peptides that contain one or more
#   non-native residues or that span a defined native → non-native junction.
#   This function serves as the initial screening step for locating potentially
#   immunogenic peptides derived from engineered or non-human sequence content.
#
# Inputs:
#   sequence            — Full-length protein sequence (single-letter amino acids).
#
#   non_human_regions   — List of fully non-native intervals, where each element
#                         is a vector c(start, end) using 1-based indexing.
#
#   mutated_positions    — Integer vector of individual positions already known
#                          to encode non-native amino acids (engineered mutations).
#
#   junctions           — Integer vector indexing junction boundaries; each entry
#                         corresponds to the residue immediately *before* a
#                         native → non-native transition. A peptide spanning
#                         j and j+1 is considered junction-spanning.
#
#   min_len, max_len    — Lower and upper bounds of peptide lengths to generate.
#                         For this de-immunization study, only 15mers were used.
#
#   include_all         — Debugging flag. If TRUE, returns *all* peptides
#                         regardless of whether they contain non-native content.
#
# Behavior:
#   • Slides a peptide window across the sequence for all lengths in the range.
#   • Flags peptides if:
#         (i) they contain ≥1 non-native residue, OR
#        (ii) they span a junction between native and/or non-native sequences.
#   • Optionally bypasses filtering (include_all = TRUE).
#
# Returns:
#   A dataframe with one row per flagged peptide containing:
#       peptide          — Amino acid sequence of the peptide window
#       start, end       — 1-based coordinates of the window in the protein
#       length           — Peptide length
#       nonhuman_count   — Number of non-native residues in the peptide
#       spans_junction   — TRUE if peptide crosses a junction boundary
#
# Notes:
#   - Function is order-preserving with respect to sequence traversal.
#------------------------------------------------------------------------------

generate_nonNative_peptides <- function(sequence,
                              non_human_regions = list(),
                              mutated_positions = integer(0),
                              junctions = integer(0),
                              min_len = 15,
                              max_len = 15,
                              include_all = FALSE) {
  
  n <- nchar(sequence)
  
  # Flag non-human positions
  is_nonhuman <- rep(FALSE, n)
  for (region in non_human_regions) {
    is_nonhuman[region[1]:region[2]] <- TRUE
  }
  is_nonhuman[mutated_positions] <- TRUE
  
  peptides <- list()
  
  for (len in min_len:max_len) {
    for (i in 1:(n - len + 1)) {
      j <- i + len - 1
      pep <- substr(sequence, i, j)
      nh_count <- sum(is_nonhuman[i:j])
      
      # Check for junction-spanning
      spans_junction <- any(sapply(junctions, function(jx) {
        i <= jx && j >= (jx + 1)
      }))
      
      # Include all if requested, otherwise filter
      if (include_all || nh_count >= 1 || spans_junction) {
        peptides[[length(peptides) + 1]] <- data.frame(
          peptide = pep,
          start = i,
          end = j,
          length = len,
          nonhuman_count = nh_count,
          spans_junction = spans_junction,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(peptides) > 0) {
    result <- do.call(rbind, peptides)
  } else {
    result <- data.frame()
  }
  
  return(result)
}

### ---------------------------------------------------------------------------
### Helper functions (should not need to be called manually)
### ---------------------------------------------------------------------------

#Identifies  occurrences of peptide sequences within a protein and
#returns indexed span. Used to map peptides predicted immunogenic onto a given protein sequence
#
#sequence = Protein sequence (single-letter amino acids) to be searched.
#peptides = Character vector of peptide sequences to locate within the protein.
#
#Returns a list where each entry corresponds to one matched peptide and contains:
#  $peptide  — The matched peptide sequence
#  $start    — 1-based start position in the protein
#  $end      — 1-based end position in the protein
#
#Note: Multiple occurrences of the same peptide are each returned as separate entries.

find_peptide_spans <- function(sequence, peptides) {
  spans <- list()
  for (pep in peptides) {
    match_pos <- unlist(gregexpr(pep, sequence, fixed = TRUE))
    for (start in match_pos[match_pos != -1]) {
      spans[[length(spans) + 1]] <- list(
        peptide = pep,
        start   = start,
        end     = start + nchar(pep) - 1
      )
    }
  }
  return(spans)
}

#Merges overlapping peptide spans into contiguous immunogenic regions for mutagenisis
#
#spans = Output list from find_peptide_spans(), where each entry contains
#        start/end coordinates for a matched peptide.
#
#Returns a dataframe with merged regions, containing:
#  $start  — Start position of merged region
#  $end    — End position of merged region
#  $peptide (unused beyond merging; retained for compatibility)
#
#Regions are merged when their intervals overlap

merge_regions <- function(spans) {
  if (length(spans) == 0) {
    return(data.frame(start = integer(0), end = integer(0)))
  }
  
  # Bind spans into a single data.frame, explicitly avoiding factors
  df <- do.call(
    rbind,
    lapply(spans, function(x) as.data.frame(x, stringsAsFactors = FALSE))
  )
  
  # Coerce start/end to integers to avoid factor/character issues
  df$start <- as.integer(df$start)
  df$end   <- as.integer(df$end)
  
  # Order by start (then end) to make merging deterministic
  df <- df[order(df$start, df$end), , drop = FALSE]
  
  merged <- list()
  current <- df[1, , drop = FALSE]
  
  if (nrow(df) > 1) {
    for (i in 2:nrow(df)) {
      row <- df[i, , drop = FALSE]
      
      # If this span overlaps or directly touches the current region, merge
      if (row$start <= current$end) {
        current$end <- max(current$end, row$end)
      } else {
        merged[[length(merged) + 1]] <- current
        current <- row
      }
    }
  }
  
  # Add the last region
  merged[[length(merged) + 1]] <- current
  
  merged_df <- do.call(rbind, merged)
  rownames(merged_df) <- NULL
  
  return(merged_df)
}

#Classifies an immunogenic region based on the origin of its non-native residues. Primarily used for debugging
#
#region = A single merged region (start/end positions).
#nonhuman_regions = List of non-human sequence intervals (c(start,end)).
#mutated_positions = Integer vector of positions containing engineered mutations.
#junctions = Integer vector indicating native–non-native junction positions.
#
#Returns one of:
#  "junction" — Region spans a known junction boundary
#  "mutation" — Region contains at least one engineered mutation
#  "nonhuman" — Region overlaps a fully non-human sequence segment
#  NA         — Region contains only native residues

get_region_reason <- function(region, nonhuman_regions, mutated_positions, junctions) {
  positions <- region$start:region$end
  
  if (any(sapply(junctions, function(jx) region$start <= jx && region$end >= (jx + 1)))) {
    return("junction")
  } else if (any(positions %in% mutated_positions)) {
    return("mutation")
  } else if (any(sapply(nonhuman_regions, function(rg) any(positions %in% rg[1]:rg[2])))) {
    return("nonhuman")
  }
  return(NA)
}


#Determines whether a residue position in the protein is native (i.e., not
#non-human, mutated, or at a junction boundary).
#
#pos = Position in the protein sequence to evaluate.
#nonhuman_regions = List of non-human sequence intervals.
#mutated_positions = Positions of mutated residues.
#junctions = Positions immediately before a junction to non-native sequence.
#
#Returns TRUE if residue is native, FALSE if it is non-native.
is_native_position <- function(pos, nonhuman_regions, mutated_positions, junctions) {
  is_nonhuman <- any(sapply(nonhuman_regions, function(rg) pos %in% rg[1]:rg[2]))
  is_mutated  <- pos %in% mutated_positions
  is_junction <- pos %in% junctions | (pos - 1) %in% junctions
  return(!(is_nonhuman || is_mutated || is_junction))
}

#Builds upstream or downstream sequence flanks containing exactly 14 consecutive
#native residues surrounding a mutated/non-native region.
#Used to evaluate the effect of introducing potentially presenting-suppressing 
#mutations on overall protein sequence
#
#direction = "upstream" or "downstream".
#sequence = Wild-type protein sequence as a character vector.
#flank_start_pos = Position immediately outside the non-native span:
#                  upstream  → first_non_native - 1
#                  downstream → last_non_native + 1
#head_start = Number of consecutive native residues already present inside the
#             region boundary (counts toward the 14 total required).
#nonhuman_regions = List of non-human regions.
#mutated_positions = Positions of mutated residues.
#junctions = Positions corresponding to junctions.
#return_debug = If TRUE, records descriptive logs of flank assembly.
#
#Returns a list with:
#  $seq   — The flank sequence (string)
#  $debug — Debugging log of the flank-building process (if enabled)
#
#Flank building terminates once 14 consecutive native residues are accumulated
#or sequence boundaries are reached.
build_flank <- function(direction,
                        sequence,
                        flank_start_pos,
                        head_start,
                        nonhuman_regions,
                        mutated_positions,
                        junctions,
                        return_debug = TRUE) {
  
  stride       <- if (direction == "upstream") -1 else 1
  flank        <- character(0)
  native_streak <- head_start
  current_pos  <- flank_start_pos
  debug        <- character(0)
  n            <- length(sequence)
  
  while (current_pos >= 1 && current_pos <= n) {
    aa <- sequence[current_pos]
    flank <- c(flank, aa)
    
    if (is_native_position(current_pos, nonhuman_regions, mutated_positions, junctions)) {
      native_streak <- native_streak + 1
      if (return_debug) {
        debug <- c(
          debug,
          sprintf("  ✓ Native %s @ %d (%s) — streak: %d",
                  direction, current_pos, aa, native_streak)
        )
      }
      if (native_streak >= 14) break
    } else {
      native_streak <- 0
      if (return_debug) {
        debug <- c(
          debug,
          sprintf("  ✗ Non-native %s @ %d (%s) — streak reset",
                  direction, current_pos, aa)
        )
      }
    }
    
    current_pos <- current_pos + stride
  }
  
  if (return_debug && native_streak < 14) {
    debug <- c(
      debug,
      sprintf("⚠️ Did not find 14 consecutive native residues in %s flank. Final streak: %d",
              direction, native_streak)
    )
  }
  
  if (direction == "upstream") flank <- rev(flank)
  return(list(seq = paste0(flank, collapse = ""), debug = debug))
}

### ---------------------------------------------------------------------------
### Main: mutation generation and context building
### ---------------------------------------------------------------------------

#Generates all allowable single- and double-amino-acid substitutions within a
#defined immunogenic region and constructs context windows with 14-residue native
#flanks on each side.
#
#region = Dataframe row containing $start and $end defining the region.
#sequence = Full-length protein sequence (string).
#reason = Classification of region ("nonhuman", "mutation", or "junction").
#allowed_pos =  Logical vector indicating which positions in the region may be mutated.
#               (Non-native residues are generally considered essential and thus cannot be mutated)
#max_mutations = Maximum number of simultaneous mutations to generate (1 or 2).
#nonhuman_regions = List of regions containing non-human residues.
#mutated_positions = Positions already known to be mutated.
#junctions = Positions marking junction boundaries.
#return_debug = Include detailed mutation and flank-building logs.
#
#Returns a dataframe with:
#  region_start, region_end   — Boundaries of the region
#  mutated_region             — Resulting mutated segment
#  context_sequence           — Full window: upstream flank + mutated region + downstream flank 
#                               this will be used to generate new peptide lists for Maria analysis
#  mutations                  — Text representation of substitutions applied
#  reason                     — Region classification
#  debug                      — Optional debugging string

generate_mutations <- function(region, sequence, reason, allowed_pos, max_mutations,
                               nonhuman_regions, mutated_positions, junctions,
                               return_debug = TRUE) {
  aa_options <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N",
                  "P", "Q", "R", "S", "T", "V", "W", "Y")
  original         <- unlist(strsplit(sequence, ""))
  region_positions <- region$start:region$end
  allowed_positions <- region_positions[allowed_pos]
  
  variants <- list()
  
  ## Helper: count native residues in a vector of positions
  count_native <- function(pos_vec) {
    pos_vec <- as.integer(pos_vec)
    pos_vec <- pos_vec[!is.na(pos_vec)]
    if (length(pos_vec) == 0) return(0L)
    
    native_flags <- vapply(
      pos_vec,
      function(p) {
        p >= 1 &&
          p <= length(original) &&
          is_native_position(p, nonhuman_regions, mutated_positions, junctions)
      },
      logical(1L)
    )
    
    sum(native_flags)
  }
  
  for (k in 1:min(max_mutations, length(allowed_positions))) {
    combos <- combn(allowed_positions, k, simplify = FALSE)
    
    for (combo in combos) {
      mut_sets <- expand.grid(
        lapply(combo, function(pos) setdiff(aa_options, original[pos])),
        stringsAsFactors = FALSE
      )
      
      for (i in 1:nrow(mut_sets)) {
        new_seq   <- original
        mutations <- list()
        
        for (j in seq_along(combo)) {
          idx    <- combo[j]
          new_aa <- mut_sets[i, j]
          if (original[idx] != new_aa) {
            mutations[[length(mutations) + 1]] <- paste0(idx, ":", original[idx], "→", new_aa)
            new_seq[idx] <- new_aa
          }
        }
        
        region_seq <- paste0(new_seq[region$start:region$end], collapse = "")
        if (region_seq == paste0(original[region$start:region$end], collapse = "")) next
        
        region_range <- region$start:region$end
        non_native_pos <- region_range[
          sapply(region_range, function(pos) {
            new_seq[pos] != original[pos] ||
              any(sapply(nonhuman_regions, function(rg) pos %in% rg[1]:rg[2])) ||
              pos %in% mutated_positions ||
              pos %in% junctions || (pos - 1) %in% junctions
          })
        ]
        if (length(non_native_pos) == 0) next
        
        first_non_native <- min(non_native_pos)
        last_non_native  <- max(non_native_pos)
        
        up_region   <- if (first_non_native > region$start) region$start:(first_non_native - 1) else integer(0)
        down_region <- if (last_non_native  < region$end)  (last_non_native + 1):region$end  else integer(0)
        
        native_before <- count_native(up_region)
        native_after  <- count_native(down_region)
        
        # ✅ Build flanks from OUTSIDE the region, but include internal natives via head_start
        upstream <- build_flank(
          direction        = "upstream",
          sequence         = original,
          flank_start_pos  = region$start - 1,   # outside region
          head_start       = native_before,      # internal native streak
          nonhuman_regions = nonhuman_regions,
          mutated_positions = mutated_positions,
          junctions        = junctions,
          return_debug     = return_debug
        )
        
        downstream <- build_flank(
          direction        = "downstream",
          sequence         = original,
          flank_start_pos  = region$end + 1,     # outside region
          head_start       = native_after,       # internal native streak
          nonhuman_regions = nonhuman_regions,
          mutated_positions = mutated_positions,
          junctions        = junctions,
          return_debug     = return_debug
        )
        
        context_seq <- paste0(upstream$seq, region_seq, downstream$seq)
        
        variants[[length(variants) + 1]] <- data.frame(
          region_start     = region$start,
          region_end       = region$end,
          mutated_region   = region_seq,
          context_sequence = context_seq,
          mutations        = paste(mutations, collapse = "; "),
          reason           = reason,
          debug = if (return_debug) paste0(
            sprintf("🧬 Region %d–%d mutated: [%s], non-native span: %d–%d",
                    region$start, region$end, region_seq, first_non_native, last_non_native),
            "\n🔧 Building upstream flank from position ", region$start - 1,
            "\n", paste0(upstream$debug, collapse = "\n"),
            "\n🔧 Building downstream flank from position ", region$end + 1,
            "\n", paste0(downstream$debug, collapse = "\n"),
            "\n🔬 Final upstream flank: [", upstream$seq, "]",
            "\n🔬 Final downstream flank: [", downstream$seq, "]"
          ) else NA_character_,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(variants) > 0) {
    do.call(rbind, variants)
  } else {
    data.frame(
      region_start     = integer(0),
      region_end       = integer(0),
      mutated_region   = character(0),
      context_sequence = character(0),
      mutations        = character(0),
      reason           = character(0),
      debug            = if (return_debug) character(0) else NULL,
      stringsAsFactors = FALSE
    )
  }
}

### ---------------------------------------------------------------------------
### Main: region scanning + calling generate_mutations
### ---------------------------------------------------------------------------

#Main pipeline function that identifies immunogenic regions, applies targeted
#mutagenesis, and assembles context windows for downstream MHC prediction.
#
#sequence = Full protein sequence (string).
#peptides = Character vector of peptides previously flagged as immunogenic.
#nonhuman_regions = List of regions containing non-human residues.
#mutated_positions = Positions of engineered mutations in the input sequence.
#junctions = Positions of native/non-native sequence junctions.
#max_mutations = Maximum number of substitutions allowed per region (default 2).
#output_cap = Optional limit on total number of generated variants.
#return_debug = Include detailed debugging logs in results.
#
#Pipeline workflow:
#  (1) Locate all immunogenic peptides (find_peptide_spans()).
#  (2) Merge overlapping peptides into continuous regions (merge_regions()).
#  (3) Classify each region by non-native origin (get_region_reason()).
#  (4) Generate all allowed mutation combinations (generate_mutations()).
#  (5) Construct 14-residue native context flanks for each design.
#  (6) Optionally subsample variants if exceeding output_cap.
#
#Returns a combined dataframe of all mutated variants across all regions.
silence_immunogenic_peptides <- function(sequence, peptides, nonhuman_regions, mutated_positions, junctions,
                                         max_mutations = 2, output_cap = NULL, return_debug = TRUE) {
  
  spans <- find_peptide_spans(sequence, peptides)
  if (length(spans) == 0) return(data.frame())
  
  merged_regions <- merge_regions(spans)
  
  count_so_far <- 0
  results      <- list()
  
  stop_if_capped <- function() {
    !is.null(output_cap) && count_so_far >= output_cap
  }
  
  for (i in 1:nrow(merged_regions)) {
    if (stop_if_capped()) break
    
    region <- merged_regions[i, ]
    reason <- get_region_reason(region, nonhuman_regions, mutated_positions, junctions)
    if (is.na(reason)) next
    
    region_positions <- region$start:region$end
    allowed <- switch(
      reason,
      "nonhuman" = sapply(region_positions, function(p) {
        any(sapply(nonhuman_regions, function(rg) p %in% rg[1]:rg[2]))
      }),
      "mutation" = !(region_positions %in% mutated_positions),
      "junction" = rep(TRUE, length(region_positions))
    )
    
    region_df <- generate_mutations(
      region           = region,
      sequence         = sequence,
      reason           = reason,
      allowed_pos      = allowed,
      max_mutations    = max_mutations,
      nonhuman_regions = nonhuman_regions,
      mutated_positions = mutated_positions,
      junctions        = junctions,
      return_debug     = return_debug
    )
    
    if (nrow(region_df) > 0) {
      if (!is.null(output_cap) && count_so_far + nrow(region_df) > output_cap) {
        remaining <- output_cap - count_so_far
        region_df <- region_df[sample(nrow(region_df), remaining), ]
      }
      results[[length(results) + 1]] <- region_df
      count_so_far <- count_so_far + nrow(region_df)
    }
  }
  
  if (length(results) == 0) {
    return(data.frame())
  } else {
    return(do.call(rbind, results))
  }
}

# High-performance peptide generator for large protein datasets.
#
# This function takes a list of protein sequences and generates all possible
# peptides of specified lengths (min_len to max_len). It is optimized for
# very large inputs by: (1) processing sequences in batches, (2) parallelizing
# computation across all available CPU cores, and (3) streaming output
# directly to disk rather than storing all peptides in memory.
#
# ARGUMENTS:
#   sequences    – Character vector of full-length protein sequences.
#                  Should be the output context sequences from silence_immunogenic_peptides
#   output_file  – Path to a tab-separated file where peptides are written.
#   min_len      – Minimum peptide length to generate (default: 15).
#   max_len      – Maximum peptide length to generate (default: 15).
#   batch_size   – Number of sequences to process per batch before writing
#                  to disk. Controls memory usage (default: 1000 sequences).
#
# OUTPUT FORMAT:
#   The output file is written in tab-separated form with column names:
#       peptide           – Generated peptide sequence
#       start             – 1-based start position of peptide in source protein
#       end               – 1-based end position
#       length            – Peptide length
#       context_sequence  – The full protein sequence from which the peptide
#                           was extracted (useful downstream for mapping).
#
# INTERNAL WORKFLOW:
#   1. Compute total peptides expected across all sequences (for ETA reporting).
#   2. Split input sequences into batches of batch_size.
#   3. Divide each batch evenly across all physical CPU cores.
#   4. For each core, generate all peptides for its subset of sequences.
#   5. Collect results, bind into a data.table, and append to the output file.
#   6. Stream progress updates including ETA and total peptides written.
#
# PARALLELIZATION:
#   - Uses parLapply() across detectCores(logical = FALSE) physical cores.
#   - Each worker receives only the sequences necessary for its subset.
#   - Output is merged batch-wise to minimize RAM usage.
#
# NOTES:
#   • Memory usage scales with batch_size, not total number of sequences.
#   • Output is written incrementally for reliability and to avoid RAM overflow.

generate_peptides_from_list_parallel_write_live <- function(sequences, output_file, min_len = 15, max_len = 15, batch_size = 1000) {
  
  if (file.exists(output_file)) file.remove(output_file)
  
  total_peptides <- sum(sapply(sequences, function(seq) {
    n <- nchar(seq)
    sum(sapply(min_len:max_len, function(len) max(0, n - len + 1)))
  }))
  
  cat("Total peptides to generate:", total_peptides, "\n")
  flush.console()
  
  start_time <- Sys.time()
  total_written <- 0
  
  ncores <- detectCores(logical = FALSE)
  
  sequence_batches <- split(sequences, ceiling(seq_along(sequences) / batch_size))
  batch_offsets <- (seq_along(sequence_batches) - 1) * batch_size
  
  split_batch_for_cores <- function(batch, start_idx, ncores) {
    indices <- start_idx + seq_along(batch) - 1
    df <- data.frame(sequence = batch, index = indices, stringsAsFactors = FALSE)
    split(df, rep(1:ncores, length.out = nrow(df)))
  }
  
  generate_batch <- function(batch_df, min_len, max_len) {
    result <- vector("list", sum(sapply(batch_df$sequence, function(seq) {
      n <- nchar(seq)
      sum(sapply(min_len:max_len, function(len) max(0, n - len + 1)))
    })))
    
    k <- 1
    for (i in seq_len(nrow(batch_df))) {
      sequence <- batch_df$sequence[i]
      n <- nchar(sequence)
      for (len in min_len:max_len) {
        for (j in 1:(n - len + 1)) {
          result[[k]] <- list(
            peptide = substr(sequence, j, j + len - 1),
            start = j,
            end = j + len - 1,
            length = len,
            context_sequence = sequence  # ✅ Save the full sequence
          )
          k <- k + 1
        }
      }
    }
    return(result)
  }
  
  cl <- makeCluster(ncores)
  on.exit(stopCluster(cl))
  
  clusterExport(cl, varlist = c("generate_batch", "min_len", "max_len"), envir = environment())
  clusterEvalQ(cl, { library(data.table) })
  
  tryCatch({
    header_written <- FALSE
    
    for (i in seq_along(sequence_batches)) {
      batch_start <- Sys.time()
      
      batch_split <- split_batch_for_cores(sequence_batches[[i]], start_idx = batch_offsets[i], ncores = ncores)
      
      batch_data_chunks <- parLapply(cl, seq_along(batch_split), function(j) {
        generate_batch(batch_split[[j]], min_len = min_len, max_len = max_len)
      })
      
      batch_data <- do.call(c, batch_data_chunks)
      
      dt <- rbindlist(lapply(batch_data, as.data.frame), use.names = TRUE, fill = TRUE)
      
      # Write header if not yet written
      fwrite(dt, output_file, sep = "\t", append = header_written, col.names = !header_written)
      header_written <- TRUE  # Only set after first batch
      
      batch_written <- nrow(dt)
      total_written <- total_written + batch_written
      
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      eta <- round(elapsed / total_written * (total_peptides - total_written), 1)
      
      cat(sprintf("Batch %d/%d complete — %d peptides written (total so far: %d) — ETA: %.1f sec\n",
                  i, length(sequence_batches), batch_written, total_written, eta))
      flush.console()
    }
    
    total_time <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)
    cat("All batches complete in", total_time, "seconds. Output written to:", output_file, "\n")
    flush.console()
    
  }, interrupt = function(e) {
    cat("\n⚠️ Peptide generation interrupted by user. Output file may be partial.\n")
    flush.console()
  })
}


combine_chunked_outputs <- function(basename, folder_path) {
  expected_cols <- c(
    "Allele1",
    "Allele2 (Same as Allele1 if analyzing a single allele)",
    "Gene Symbol",
    "Peptide Sequence",
    "TPM (Optional)",
    "TPM estimated",
    "MARIA raw scores",
    "MARIA percentile scores",
    "15mer core",
    "Positive presenters"
  )
  
  pattern <- paste0("^", basename, "_chunk_\\d+\\.txt\\.output\\.txt$")
  files <- list.files(path = folder_path, pattern = pattern, full.names = TRUE)
  
  if (length(files) == 0) {
    stop("No matching files found.")
  }
  
  cat("Found", length(files), "files. Combining...\n")
  
  combined_df <- rbindlist(lapply(files, function(file) {
    dt <- fread(file, sep = "\t", header = TRUE, fill = TRUE)
    
    if (!all(expected_cols %in% colnames(dt))) {
      cat("⚠️ Warning: File", basename(file), "does not match expected columns. Forcing column names and skipping first row.\n")
      # Now read again, header = FALSE, col.names = expected_cols, skip first line
      dt <- fread(file, sep = "\t", header = FALSE, fill = TRUE, skip = 1, col.names = expected_cols)
    }
    
    dt
  }), use.names = TRUE, fill = TRUE)
  
  cat("Combined", nrow(combined_df), "rows from", length(files), "files.\n")
  
  return(combined_df)
}

count_immunogenic_per_sequence <- function(df) {
  
  # Ensure 'Positive presenters' is numeric
  df$`Positive presenters` <- suppressWarnings(as.numeric(df$`Positive presenters`))
  
  # Check for NAs introduced by coercion
  if (any(is.na(df$`Positive presenters`) & !is.na(df$`Positive presenters`))) {
    warning("⚠️ Warning: NAs introduced in 'Positive presenters' during numeric conversion.")
  }
  
  # Create a grouping variable that preserves order
  df$..row_order <- match(df$`Gene Symbol`, unique(df$`Gene Symbol`))
  
  # Summarize using dplyr and preserve order
  result <- df %>%
    group_by(`Gene Symbol`, ..row_order) %>%
    summarize(`Immunogenic Peptide Count` = sum(`Positive presenters`, na.rm = TRUE), .groups = "drop") %>%
    arrange(..row_order) %>%
    select(-..row_order)
  
  return(as.data.frame(result))
}

append_immunogenic_counts <- function(original_df, summary_df) {
  # Initialize the new column to 10 by default
  original_df$`Immunogenic Peptide Count` <- 10
  
  # Match using context sequence as key
  match_indices <- match(original_df$context_sequence, summary_df$`Gene Symbol`)
  
  # Fill only for matches found
  matched <- !is.na(match_indices)
  original_df$`Immunogenic Peptide Count`[matched] <- summary_df$`Immunogenic Peptide Count`[match_indices[matched]]
  
  return(original_df)
}

filter_lowest_immunogenic_mutations <- function(df, collect = TRUE) {
  # ensure the column is numeric
  df$`Immunogenic Peptide Count` <- as.numeric(df$`Immunogenic Peptide Count`)
  
  filtered_df <- df %>% 
    dplyr::group_by(region_start, region_end) %>% 
    dplyr::group_modify(~{
      min_count <- min(.x$`Immunogenic Peptide Count`, na.rm = TRUE)
      
      if (collect) {
        # NEW: return *all* rows that share the minimum count (0 if present)
        .x %>% dplyr::filter(`Immunogenic Peptide Count` == min_count)
        
      } else {
        # OLD logic: one representative row
        if (min_count == 0) {
          .x %>% 
            dplyr::filter(`Immunogenic Peptide Count` == 0) %>% 
            dplyr::slice(1)
        } else {
          .x %>% dplyr::slice_min(`Immunogenic Peptide Count`, with_ties = FALSE)
        }
      }
    }) %>% 
    dplyr::ungroup()
  
  as.data.frame(filtered_df)
}

#------------------------------------------------------------------------------
# Purpose:
#   Combines multiple MARIA output chunk files (e.g., *_chunk_1.txt.output.txt)
#   back into a single unified results table. This function is used after 
#   splitting large input files into chunks and running MARIA on each chunk.
#
# Inputs:
#   basename     — Base filename prefix used when generating chunk outputs 
#                   (e.g., "ALTER4_MV1").
#   folder_path  — Directory containing the chunked *.output.txt files.
#
# Behavior:
#   • Scans the directory for files matching "<basename>_chunk_<n>.txt.output.txt"
#   • Loads each file with fread()
#   • Ensures the expected MARIA columns are present
#     - If columns do not match, warns and forcibly assigns column names
#   • Row-binds all cleaned chunk files into one combined data.frame
#
# Returns:
#   A fully combined data.table containing all rows from all chunked outputs.
#------------------------------------------------------------------------------

combine_chunked_outputs <- function(basename, folder_path) {
  expected_cols <- c(
    "Allele1",
    "Allele2 (Same as Allele1 if analyzing a single allele)",
    "Gene Symbol",
    "Peptide Sequence",
    "TPM (Optional)",
    "TPM estimated",
    "MARIA raw scores",
    "MARIA percentile scores",
    "15mer core",
    "Positive presenters"
  )
  
  pattern <- paste0("^", basename, "_chunk_\\d+\\.txt\\.output\\.txt$")
  files <- list.files(path = folder_path, pattern = pattern, full.names = TRUE)
  
  if (length(files) == 0) {
    stop("No matching files found.")
  }
  
  cat("Found", length(files), "files. Combining...\n")
  
  combined_df <- rbindlist(lapply(files, function(file) {
    dt <- fread(file, sep = "\t", header = TRUE, fill = TRUE)
    
    if (!all(expected_cols %in% colnames(dt))) {
      cat("⚠️ Warning: File", basename(file), "does not match expected columns. Forcing column names and skipping first row.\n")
      # Now read again, header = FALSE, col.names = expected_cols, skip first line
      dt <- fread(file, sep = "\t", header = FALSE, fill = TRUE, skip = 1, col.names = expected_cols)
    }
    
    dt
  }), use.names = TRUE, fill = TRUE)
  
  cat("Combined", nrow(combined_df), "rows from", length(files), "files.\n")
  
  return(combined_df)
}

#------------------------------------------------------------------------------

# Purpose:
#   Computes the total number of immunogenic peptides per input sequence.
#   This corresponds to summing the “Positive presenters” field for each
#   gene/protein sequence (MARIA output column: "Gene Symbol").
#
# Inputs:
#   df — A MARIA results dataframe containing at least:
#          • "Gene Symbol"
#          • "Positive presenters" (numeric or character coercible to numeric)
#
# Behavior:
#   • Groups rows by Gene Symbol
#   • Sums all positive presenters assigned to each sequence
#   • Preserves original ordering of Gene Symbols
#
# Returns:
#   A dataframe with:
#       Gene Symbol
#       Immunogenic Peptide Count   — sum of `Positive presenters`
#
#------------------------------------------------------------------------------

count_immunogenic_per_sequence <- function(df) {
  
  # Ensure 'Positive presenters' is numeric
  df$`Positive presenters` <- suppressWarnings(as.numeric(df$`Positive presenters`))
  
  # Check for NAs introduced by coercion
  if (any(is.na(df$`Positive presenters`) & !is.na(df$`Positive presenters`))) {
    warning("⚠️ Warning: NAs introduced in 'Positive presenters' during numeric conversion.")
  }
  
  # Create a grouping variable that preserves order
  df$..row_order <- match(df$`Gene Symbol`, unique(df$`Gene Symbol`))
  
  # Summarize using dplyr and preserve order
  result <- df %>%
    group_by(`Gene Symbol`, ..row_order) %>%
    summarize(`Immunogenic Peptide Count` = sum(`Positive presenters`, na.rm = TRUE), .groups = "drop") %>%
    arrange(..row_order) %>%
    select(-..row_order)
  
  return(as.data.frame(result))
}

#------------------------------------------------------------------------------
# Purpose:
#   Adds an "Immunogenic Peptide Count" column to a mutation results dataframe 
#   by matching each mutated design to its total immunogenicity score.
#
# Inputs:
#   original_df  — Output of the mutation-generation pipeline containing:
#                    • context_sequence (sequence identifier)
#   summary_df   — Output of count_immunogenic_per_sequence() with:
#                    • Gene Symbol
#                    • Immunogenic Peptide Count
#
# Behavior:
#   • Initializes all rows to a default value of 1000 (easy debugging) if a
#     given row cannot be matched
#   • Matches context_sequence → Gene Symbol via match()
#   • Populates only the rows where a matching summary entry exists
#
# Returns:
#   The original_df with a new "Immunogenic Peptide Count" column added.
#
# Notes:
#   - Matching assumes each context_sequence corresponds to a unique gene identifier
#------------------------------------------------------------------------------

append_immunogenic_counts <- function(original_df, summary_df) {
  # Initialize the new column to 10 by default
  original_df$`Immunogenic Peptide Count` <- 1000
  
  # Match using context sequence as key
  match_indices <- match(original_df$context_sequence, summary_df$`Gene Symbol`)
  
  # Fill only for matches found
  matched <- !is.na(match_indices)
  original_df$`Immunogenic Peptide Count`[matched] <- summary_df$`Immunogenic Peptide Count`[match_indices[matched]]
  
  return(original_df)
}

#------------------------------------------------------------------------------
# Purpose:
#   Identifies the mutation designs within each region (region_start/region_end)
#   that minimize predicted immunogenicity. This is the last stage of the 
#   de-immunization pipeline, selecting the best variants to carry forward.
#
# Inputs:
#   df       — Mutation results dataframe containing:
#                • region_start / region_end
#                • Immunogenic Peptide Count (numeric)
#   collect  — If TRUE: return *all* variants that achieve the minimum 
#                       immunogenicity per region.
#              If FALSE: return only a single representative variant per region.
#
# Behavior:
#   • Groups mutations by region (start/end)
#   • Within each region, identifies the minimum immunogenicity score
#   • Returns:
#       - all minimum-score variants (collect=TRUE), OR
#       - one representative row (collect=FALSE)
#
# Returns:
#   A dataframe of the selected minimum-immunogenicity variants.
#
#------------------------------------------------------------------------------

filter_lowest_immunogenic_mutations <- function(df, collect = TRUE) {
  # ensure the column is numeric
  df$`Immunogenic Peptide Count` <- as.numeric(df$`Immunogenic Peptide Count`)
  
  filtered_df <- df %>% 
    dplyr::group_by(region_start, region_end) %>% 
    dplyr::group_modify(~{
      min_count <- min(.x$`Immunogenic Peptide Count`, na.rm = TRUE)
      
      if (collect) {
        # NEW: return *all* rows that share the minimum count (0 if present)
        .x %>% dplyr::filter(`Immunogenic Peptide Count` == min_count)
        
      } else {
        # OLD logic: one representative row
        if (min_count == 0) {
          .x %>% 
            dplyr::filter(`Immunogenic Peptide Count` == 0) %>% 
            dplyr::slice(1)
        } else {
          .x %>% dplyr::slice_min(`Immunogenic Peptide Count`, with_ties = FALSE)
        }
      }
    }) %>% 
    dplyr::ungroup()
  
  as.data.frame(filtered_df)
}

#------------------------------------------------------------------------------
# apply_mutations_from_filtered_df()
#
# Purpose:
#   Applies the set of selected mutations returned by
#   filter_lowest_immunogenic_mutations(collect = FALSE) to an initial
#   wild-type protein sequence. Updates the list of mutated positions by
#   adding any newly altered residues.
#
# Inputs:
#   sequence    — Full-length input protein sequence (single-letter AA).
#
#   mutations_df — Dataframe (typically from filter_lowest_immunogenic_mutations)
#                  containing a "mutations" column formatted as:
#                       "31:A→G; 35:F→Y"
#                  where each mutation is:
#                       <pos>:<oldAA>→<newAA>   (1-based positions)
#
#   original_mutated_positions — Vector of existing mutated positions. These
#                  will be carried forward and expanded with newly applied
#                  mutations.
#
# Behavior:
#   • Parses mutation strings.
#   • Verifies the reference amino acid when possible.
#   • Applies each mutation to the sequence.
#   • Tracks all successfully applied mutation positions.
#   • Merges the new mutated positions with the original list, de-duplicated.
#
# Returns a list with:
#   $mutated_sequence         — Final amino-acid sequence after mutation.
#   $mutation_log             — Dataframe of all parsed and applied mutations.
#   $updated_mutated_positions — Union of original_mutated_positions and all
#                                 newly mutated residue indices.
#------------------------------------------------------------------------------
apply_mutations_from_filtered_df <- function(sequence,
                                             mutations_df,
                                             original_mutated_positions = integer(0)) {
  
  # Convert input sequence to vector
  seq_vec <- unlist(strsplit(sequence, split = ""))
  seq_len <- length(seq_vec)
  
  log_list <- list()
  newly_mutated_positions <- integer(0)  # track new mutations
  
  # Parse mutation like "31:A→G"
  parse_one_mutation <- function(mut_str) {
    m <- regexec("^(\\d+):([A-Z])\u2192([A-Z])$", mut_str)
    parts <- regmatches(mut_str, m)[[1]]
    if (length(parts) != 4) {
      warning(sprintf("Could not parse mutation string: '%s'", mut_str))
      return(NULL)
    }
    list(
      position = as.integer(parts[2]),
      from_aa  = parts[3],
      to_aa    = parts[4]
    )
  }
  
  # Iterate through all selected region-level mutations
  for (i in seq_len(nrow(mutations_df))) {
    mut_field <- mutations_df$mutations[i]
    
    if (is.na(mut_field) || !nzchar(mut_field)) next
    
    mut_strings <- strsplit(mut_field, ";\\s*")[[1]]
    
    for (ms in mut_strings) {
      if (!nzchar(ms)) next
      
      parsed <- parse_one_mutation(ms)
      if (is.null(parsed)) next
      
      pos     <- parsed$position
      from_aa <- parsed$from_aa
      to_aa   <- parsed$to_aa
      
      if (pos < 1 || pos > seq_len) {
        warning(sprintf("Position %d out of bounds (len=%d). Skipping.", pos, seq_len))
        next
      }
      
      # check original AA
      if (seq_vec[pos] != from_aa) {
        warning(sprintf(
          "Mismatch at pos %d: sequence has '%s', mutation expects '%s'. Applying mutation anyway.",
          pos, seq_vec[pos], from_aa
        ))
      }
      
      # Apply mutation
      seq_vec[pos] <- to_aa
      
      # Track new mutation if not previously mutated
      newly_mutated_positions <- c(newly_mutated_positions, pos)
      
      # Record log entry
      log_list[[length(log_list) + 1]] <- data.frame(
        region_start    = mutations_df$region_start[i],
        region_end      = mutations_df$region_end[i],
        position        = pos,
        from_aa         = from_aa,
        to_aa           = to_aa,
        mutation_string = ms,
        stringsAsFactors = FALSE
      )
    }
  }
  
  mutated_sequence <- paste0(seq_vec, collapse = "")
  
  # Merge old + new mutated positions
  updated_mutated_positions <- sort(unique(c(original_mutated_positions,
                                             newly_mutated_positions)))
  
  mutation_log <- if (length(log_list) > 0) do.call(rbind, log_list) else data.frame()
  
  list(
    mutated_sequence          = mutated_sequence,
    mutation_log              = mutation_log,
    updated_mutated_positions = updated_mutated_positions
  )
}

#------------------------------------------------------------------------------
# Purpose:
#   Collects and merges a selected subset of MARIA output dataframes, assigns
#   each table a new Gene.Symbol value, and returns a unified dataframe. Used
#   when multiple input sequences are processed and must be relabeled for
#   downstream immunogenicity analysis.
#
# Inputs:
#   dfs            — Named list of dataframes containing MARIA outputs.
#   subset_names   — Character vector of list names specifying which dfs to extract.
#   gene_symbols   — Character vector of gene/sequence identifiers to assign
#                    to the corresponding dataframes. Must be same length as
#                    subset_names.
#
# Behavior:
#   • Normalizes column names using make.names()
#   • Extracts only the requested dataframes, preserving order
#   • Replaces the Gene.Symbol column of each dataframe with the specified
#     annotation from gene_symbols
#   • Row-binds all updated tables into a single combined dataframe
#
# Returns:
#   A merged dataframe containing all selected MARIA outputs with unified and
#   reassigned Gene.Symbol identifiers.
#------------------------------------------------------------------------------

dataCollect <- function(dfs, subset_names, gene_symbols) {
  # Validate inputs
  if (length(subset_names) != length(gene_symbols)) {
    stop("subset_names and gene_symbols must have the same length.")
  }
  
  dfs <- lapply(dfs, function(x){names(x) <- make.names(names(x)); x})
  
  # Extract the selected dataframes in the desired order
  selected_dfs <- dfs[subset_names]
  
  updated_dfs <- Map(function(df, symbol) {
    df$Gene.Symbol <- as.character(df$Gene.Symbol)  # convert column
    df$Gene.Symbol <- as.character(symbol)          # assign new values
    df
  }, selected_dfs, gene_symbols)
  
  # Combine all dataframes row-wise
  combined_df <- do.call(rbind, updated_dfs)
  
  return(combined_df)
}

#------------------------------------------------------------------------------
#
# Purpose:
#   Produces a jittered scatter plot of MARIA percentile scores for each
#   sequence variant, visualizing the distribution of predicted MHC-II
#   presentation values. This is used to compare immunogenic burden across
#   candidate de-immunized sequences.
#
# Inputs:
#   MARIA_sheet — Dataframe containing (at minimum):
#                   • Gene.Symbol
#                   • MARIA.percentile.scores
#
# Behavior:
#   • Extracts sequence identifiers and MARIA presentation percentiles
#   • Generates a ggplot scatter plot with jittered markers
#   • Draws a dotted vertical threshold line at x = 63
#   • Shades the high-presentation region (63–100%) for emphasis
#   • Applies a prism-style theme for publication-quality appearance
#
# Returns:
#   A ggplot object representing the immunogenicity distribution plot. Can be
#   displayed directly or modified/saved externally.
#------------------------------------------------------------------------------

plot_immunogenic <- function(MARIA_sheet){
  data <- MARIA_sheet %>% select(`Gene.Symbol`, `MARIA.percentile.scores`)
  
  colnames(data) <- c("y","x")
  
  p <- ggplot(data, aes(x = x, y = y)) +
    geom_jitter(shape = 23, width = 1, height = 0.1, fill = "black") +  # Jittered diamond points
    geom_vline(xintercept = 63, linetype = "dotted", color = "black") +  # Dotted line at x=63
    annotate("rect", xmin = 63, xmax = 100, ymin = -Inf, ymax = Inf,
             fill = "pink", alpha = 0.2) +  # Shaded area
    scale_x_continuous(limits = c(0, 100), name = "MARIA MHC Presentation Percentile") +
    ylab(NULL) +  # Hide y-axis label
    theme_prism()
  
  return(p)
}

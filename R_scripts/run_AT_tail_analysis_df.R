run_AT_tail_analysis_df <- function(df) {
  # Make sure required libraries are loaded
  library(stringr)
  library(stringi)
  
  extract_AT_info <- function(seq, flag, cigar, chrR2, seqnames, strand, start, end, posR2) {
    total_length <- nchar(seq)
    ##############################################################################
    # 1. calculate how many AT is in the sequence in %
    at_content <- (sum(str_count(seq, c("A", "T"))) / total_length) * 100
    # set treashold for high AT content
    high_AT_mix <- at_content > 93
    ##############################################################################
    #2. prepare R2 sequence and save it as seq2
    # Reverse sequence of R2 if flag is 16 (reverse mapping should map to genes in + strand)
    if (flag == 16) {
      seq <- stri_reverse(seq)
    } else {
      # else if flag is 4 (sense mapping should map to genes in - strand)
      # swap the bases so it is like complement sequence but do not reverse the read
      seq <- str_replace_all(seq, "A", "X")
      seq <- str_replace_all(seq, "T", "A")
      seq <- str_replace_all(seq, "X", "T")
    }
    # this way I have sequence of read where tail is always from the start of the sequence.
    # Save original sequence to extract estimated tail at the end
    seq2 <- seq
    ##############################################################################
    #3. set empty values for variables that will be extracted
    t_len <- 0
    g_len <- 0
    c_len <- 0
    a_count <- 0
    tail_len <- 0
    tail <- NA
    ##############################################################################
    # 4. look for tail composition
    # 4.1 First cases when R2 mapped so CIGAR is present (is not NA or * !!! value here
    # may depend on the aligner)
    # I'm putting 2 conditions 1. CIGAR for R2 present; 2. R2 mapped to same chromosome as R1
    if (cigar != "*" && chrR2 == seqnames) {
      # next line will establish gene stop position from annotation which depends
      # if gene is located on the - or + strand. Locations are from left to right.
      stop <- if (strand == "+") end else if (strand == "-") start else NA
      # I give distance treashold, I do not bother with R2 strendedness at this point cause 
      # distance is very big.
      if (!is.na(stop) && abs(stop - posR2) < 50000) {
        # now I will extract softclipping infor from CIGAR
        # if flag is 16 softclipping info is at the end of CIGAR
        if (flag == 16) {
          softclip_3prime <- str_extract(cigar, "\\d+(?=S$)")
        } else if (flag == 0) {
          # flag 0 means softclipping info of interest is at the beggining
          softclip_3prime <- str_extract(cigar, "^\\d+(?=S)")
        }
        # save number of softclipped nucleotides (this is my tail)
        softclip_3prime <- as.integer(softclip_3prime)
        if (!is.na(softclip_3prime)) {
          # this number I use to extract tail sequence from left side of seq2
          tail <- str_sub(seq2, 1, softclip_3prime)
          seq <- tail
          # end record number of bases total in the tail
          tail_len <- softclip_3prime
          
          # only now I am scanning extracted tail for composition
          # I always use same code for this, first I am scanning end for potential
          # non-A bases.
          #4.1.1 starting with T
          if (str_sub(seq, 1, 1) == "T") {
            t_match <- regmatches(seq, regexpr("^T{1,}", seq))
            t_len <- str_length(t_match)
            # if I found any Ts I will strip them from the sequence and keep the sequence
            seq <- str_sub(seq, t_len + 1)
            
          } else if (str_sub(seq, 1, 1) == "G") {
            # 4.1.2 Same for G
            g_match <- regmatches(seq, regexpr("^G{1,}", seq))
            g_len <- str_length(g_match)
            seq <- str_sub(seq, g_len + 1)
            
          } else if (str_sub(seq, 1, 1) == "C") {
            # 4.1.3 Same for C
            c_match <- regmatches(seq, regexpr("^C{1,}", seq))
            c_len <- str_length(c_match)
            seq <- str_sub(seq, c_len + 1)
            
          }
          # 4.1.4 At this point even seq is stripped from eventual first few non_As 
          # or if it did not met any of those conditions it is in its original form
          # now I will look for A streatch, but allowing mismatches in the tail
          a_match <- regmatches(seq, regexpr("^A{1,}(?:[^A]A{5,})*", seq))
          a_count <- ifelse(length(a_match) > 0, str_count(a_match, "A"), 0)
          
        } else {
          # this is is soft clip is NA
          tail <- "full match"
        }
      }
    } else {
      # fallback: now analyse tail if R2 did not map, so I run my scanning code
      # through R2 from the left hand side seq now is seq2
      if (str_sub(seq, 1, 1) == "T") {
        t_match <- regmatches(seq, regexpr("^T{1,}", seq))
        t_len <- str_length(t_match)
        seq <- str_sub(seq, t_len + 1)
        
      } else if (str_sub(seq, 1, 1) == "G") { g_match <- regmatches(seq, regexpr("^G{1,}", seq))
      g_len <- str_length(g_match)
      seq <- str_sub(seq, g_len + 1)
      
      } else if (str_sub(seq, 1, 1) == "C") {
        c_match <- regmatches(seq, regexpr("^C{1,}", seq))
        c_len <- str_length(c_match)
        seq <- str_sub(seq, c_len + 1)
        
      }
      a_match <- regmatches(seq, regexpr("^A{1,}(?:[^A]A{5,})*", seq))
      a_count <- ifelse(length(a_match) > 0, str_count(a_match, "A"), 0)
      tail_len <- ifelse(length(a_match) > 0, str_length(a_match), 0) + t_len + g_len + c_len
      tail <- str_sub(seq2, 1, tail_len)
    }
    
    return(c(a_count, t_len, g_len, c_len, tail_len, high_AT_mix, tail))
  }
  
  # Preallocate
  n <- nrow(df)
  A_count <- numeric(n)
  T_count <- numeric(n)
  G_count <- numeric(n)
  C_len <- numeric(n)
  tail_len <- numeric(n)
  high_AT_mix <- logical(n)
  tail_seq <- character(n)
  
  # Loop
  for (i in seq_len(n)) {
    res <- extract_AT_info(
      df$seq_R2[i],
      df$flagR2[i],
      df$CIGAR_R2[i],
      df$chrR2[i],
      df$seqnames[i],
      df$strand[i],
      df$start[i],
      df$end[i],
      df$posR2[i]
    )
    A_count[i] <- as.numeric(res[1])
    T_count[i] <- as.numeric(res[2])
    G_count[i] <- as.numeric(res[3])
    C_len[i]  <- as.numeric(res[4])
    tail_len[i] <- as.numeric(res[5])
    high_AT_mix[i] <- as.logical(res[6])
    tail_seq[i] <- as.character(res[7])
  }
  
  # Add to data frame
  df$A_count <- A_count
  df$T_count <- T_count
  df$G_count <- G_count
  df$C_len <- C_len
  df$tail_len <- tail_len
  df$high_AT_mix <- high_AT_mix
  df$tail <- tail_seq
  
  return(df)
}
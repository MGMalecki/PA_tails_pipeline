assign_type_and_class <- function(df) {
  df <- df %>%
    mutate(
      type = case_when(
        tail_len > 0 & CIGAR_R2 != "*" ~ "cigar",
        tail_len > 0 & CIGAR_R2 == "*" ~ "grep",
        tail_len == 0 & CIGAR_R2 != "*" ~ "no_tail",
        TRUE ~ "dubious"
      ),
      class = case_when(
        T_count > 0 & A_count == 0 ~ "oU",
        T_count > 0 & A_count > 0 ~ "AU",
        G_count > 0 & A_count == 0 ~ "G",
        G_count > 0 & A_count > 0 ~ "AG",
        C_len > 0 & A_count == 0 ~ "C",   # Same as G_count condition
        C_len > 0 & A_count > 0 ~ "AC",
        A_count > 0 ~ "A",
        TRUE ~ "check"
      )
    )
  return(df)
}
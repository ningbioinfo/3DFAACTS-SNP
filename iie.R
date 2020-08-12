interaction_intersecting_elements <- function(interaction, gr_elements, minoverlap = 0L){
  # left
  df <- interaction[1:6] %>%
    magrittr::set_colnames(c("chr1", "start1", "end1", "chr2", "start2", "end2"))
  
  df_left <- df[1:3] %>%
    magrittr::set_colnames(c("chr","start","end")) %>%
    makeGRangesFromDataFrame()
  
  ov_1 <- findOverlaps(df_left, gr_elements)
  
  df_1 <- df_left[ov_1@from] %>%
    as_tibble() %>%
    dplyr::select(-c(width,strand)) %>%
    magrittr::set_colnames(c("chr1","start1","end1")) %>%
    left_join(df) %>%
    unique()
  
  # right
  
  df_right <- df[4:6] %>%
    magrittr::set_colnames(c("chr","start","end")) %>%
    makeGRangesFromDataFrame()
  
  ov_2 <- findOverlaps(df_right, gr_elements)
  
  df_2 <- df_right[ov_2@from] %>%
    as_tibble() %>%
    dplyr::select(-c(width,strand)) %>%
    magrittr::set_colnames(c("chr2","start2","end2")) %>%
    left_join(df) %>%
    unique() %>%
    dplyr::select(colnames(df_1))
  
  # merge
  final <- rbind(df_1,df_2) %>%
    arrange_(~chr1, ~start1, ~chr2, ~start2) %>%
    unique()
  
  y <<- c(gr_elements[ov_1@to], gr_elements[ov_2@to]) %>%
    unique() %>%
    sort()
  
  return(final)
  
}

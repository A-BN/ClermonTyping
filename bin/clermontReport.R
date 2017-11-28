#' ---
#' title: ClermonTyping Report
#' author: ""
#' date: "`r paste0('Date: ', date())`"
#' output:
#'    html_document:
#'      toc: false
#'      highlight: textmate
#' ---

#' 
#+ echo = FALSE  
# packages loading
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(tidyr))

#+ echo = FALSE, message = FALSE, warning = FALSE
# data loading, TARTAMPION will be replaced by the rightful name of the report as inputed in report_calling_func()
clermonT <- 
  read_delim("TARTAMPION", 
             "\t", escape_double = FALSE, col_names = c("file", "internal", "quadruplex", "supp", "phylogroup", "mash_file"),
             col_types = paste0(rep("c", 6), collapse = ""), 
             trim_ws = TRUE)
# clermonT <-
#   read_delim("/home/abn/Documents/clermonTyping/Jeannot/antoine-21112017_043026_phylogroups.txt",
#             "\t", escape_double = FALSE, col_names = c("file", "internal", "quadruplex", "supp", "phylogroup", "mash_file"),
#              col_types = paste0(rep("c", 6), collapse = ""),
#              trim_ws = TRUE)


clermonT_out_dir <- dirname("TARTAMPION")
# clermonT_out_dir <- dirname("/home/abn/Documents/clermonTyping/Jeannot/antoine-21112017_043026_phylogroups.txt")

clermonT <-
  clermonT %>%
    mutate(mash_file = paste0(clermonT_out_dir,"/",mash_file))

#+ echo = FALSE, message = FALSE, warning = FALSE
# loading the mash data_frames as a list
mash_results <-
  lapply(X = clermonT$mash_file,
                       FUN = function(x) read_delim(x,
                                                    "\t", escape_double = FALSE,
                                                    trim_ws = TRUE,
                                                    col_names = c("score", "kmer", "med_mul", "pval", "strain", "empty")))

# this dplyr transformation will try to give a warning if the mash result is unclear 
# (i.e multiple results w/ score > 0.90 non concordant on the group)
mash_group <-
  lapply(X = mash_results, FUN = function(mash_df){
    mash_df %>%
      filter(score > 0.95) %>%
      arrange(desc(score)) %>%
      mutate(mash_raw = str_replace(string = strain, pattern = ".*_(.*).fasta", replacement = "\\1")) %>%
      rowwise()%>%
      mutate(mash_group = ifelse(test = length(unique(.$mash_raw)) == 1, yes = mash_raw, no = paste0(mash_raw, "*"))) %>%
      slice(1) %>%
      select(mash_group)
  })
mash_group <- as.character(unlist(mash_group))
clermonT_2 <-
  clermonT %>%
    mutate(internal = NULL)%>%
    mutate(quadruplex = str_replace_all(string = quadruplex, pattern = "[\\[\\]']", replacement = ""),
           quadruplex = str_replace_all(string = quadruplex, pattern = ",", replacement = "   ")) %>%
    mutate(supp = str_replace_all(string = supp, pattern = "[\\[\\]']", replacement = "")) %>%
    select(-mash_file)%>%
    cbind(mash_group)

kable(clermonT_2)

messaga <- c()
for(i in 1:nrow(clermonT_2)){
  curr_phylo <- clermonT_2[i, "phylogroup"]
  curr_mash <- clermonT_2[i, "mash_group"]
  curr_mash_noS <- str_replace_all(string = curr_mash, pattern = "\\*", replacement = "")
  curr_file <- clermonT_2[i, "file"]
  if(!(curr_phylo == curr_mash_noS)){
    messaga <-
      c(messaga,
        paste0("Warning in : ", curr_file,".\tThe Clermont phylogroup doesn't match the mash closest neighbor's group !"))
  }
  if(grepl(pattern = "\\*", x = curr_mash)){
    messaga <- 
      c(messaga, paste0("Warning in : ", curr_file,
        "\tA star next to the mash_group results indicates\nthat the input sequence might be a mix between multiple *E. coli* genomes\n(i.e a metagenome or a recombined genome)."))
    }
}
cat(x = messaga, sep = "\n")
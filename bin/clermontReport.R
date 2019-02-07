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
             trim_ws = TRUE, 
             quoted_na = FALSE)

clermonT_out_dir <- dirname("TARTAMPION")

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
                                                    col_names = c("score", "kmer", "med_mul", "pval", "strain", "empty"),
                                                    col_types = cols(score = col_double(), 
                                                                  kmer = col_character(), 
                                                                  med_mul = col_integer(), 
                                                                  pval = col_double(), 
                                                                  strain = col_character(), 
                                                                  empty = col_character())))

# this dplyr transformation will output "unknown" if the mash result is unclear 
# (i.e multiple results w/ score > 0.90 non concordant on the group)
mash_threshold <- 0.95
mash_group <-
  lapply(X = mash_results, FUN = function(mash_df){
    if(nrow(mash_df) != 0){
      mash_df %>%
        filter(score > mash_threshold) %>%
        arrange(desc(score)) %>%
        mutate(mash_short = str_replace(string = strain, pattern = ".*_(.*)_.*.fasta", replacement = "\\1")) %>%
        mutate(mash_full = str_replace(string = strain, pattern = ".*_(.*_.*).fasta", replacement = "\\1")) %>%
        mutate(mash_full = str_replace(string = mash_full,pattern = "_", replacement = "")) %>%
        rowwise()%>%
        mutate(mash_short = ifelse(test = length(unique(.$mash_short)) == 1, yes = mash_short, no = paste0(mash_short, "*"))) %>%
        slice(1) %>%
        select(mash_short)
    } else
      mash_df <- "NA"
  })

mash_group <- as.character(paste0("",unlist(mash_group)))
mash_group[mash_group == ""] <- "Unknown"

clermonT_2 <-
  clermonT %>%
    mutate(internal = NULL)%>%
    mutate(quadruplex = str_replace_all(string = quadruplex, pattern = "[\\[\\]']", replacement = ""),
           quadruplex = str_replace_all(string = quadruplex, pattern = ",", replacement = "   ")) %>%
    mutate(supp = str_replace_all(string = supp, pattern = "[\\[\\]']", replacement = "")) %>%
    select(-mash_file)%>%
    cbind(mash_group) %>%
    mutate(mash_group = as.character(mash_group))

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
        paste0("Warning in : ", curr_file,".\tThe Clermont phylogroup doesn't match the mash closest neighbor's group!\nThis could indicate a mutation affecting the binding of a primer"))
  }
  if(grepl(pattern = "\\*", x = curr_mash)){
    messaga <- 
      c(messaga, paste0("Warning in : ", curr_file,
        "\tThe input sequence might be a mix between multiple *E. coli* genomes\n(i.e a metagenome or a recombined genome)."))
  }
  if(grepl(pattern = "Unknown", x = curr_phylo)){
    messaga <- 
      c(messaga, paste0("Warning in : ", curr_file,
                        "\tA You might have discovered a new quadruplex genotype (or simply a mutation affecting the binding of a primer)! To be on the safe side you should consider doublechecking your data."))
  }
  if(grepl(pattern = "Unknown", x = curr_mash)){
    messaga <- 
      c(messaga, paste0("Warning in : ", curr_file,
                        "\tA There was no close match in our mash database. This should not happen with a complete E. coli genome."))
  }
}
cat(x = messaga, sep = "\n")

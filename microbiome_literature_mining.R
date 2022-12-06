#' TO DO LIST:
#' 
#' EXTRACTING THE DATA:
#' 
#' 0) Re-think on how to classify the papers in either of the cetegories: Review, meta-analysis or research article
#'   -> One alternative could be to extract it from the xml
#' 
#' 1) CHECK THAT THE ORDER OF THE PMCID GIVEN TO ID_TO_TABLE FUNCTION DOESN'T AFFECT THE OUTPUT [DONE]
#' 
#' 2) Check which studies are from human samples and which are from other animals samples
#' 
#' 
#' 
#' FIGURES FOR THE PAPER
#' 
#' 1) Make a table out of the count of how many papers per year are included in the analysis
#' 
#' 2) Tendency over the years to use clr transformation
#' Include : 
#'        https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224/full
#'        https://www.nature.com/articles/s41564-018-0337-x 
#'        https://pubmed.ncbi.nlm.nih.gov/31194939/ 
#'
#' 3) Alluvial plots
#' 
#' 4) https://cran.r-project.org/web/packages/ggupset/readme/README.html
#' 
#' 
#' 
#' GENERAL LOOKING
#' 
#' 0) Color palettes
#' https://mikemol.github.io/technique/colorblind/2018/02/11/color-safe-palette.html 
#' https://github.com/nanxstats/ggsci
#' https://github.com/EmilHvitfeldt/paletteer



# Setting everything for the project --------------------------------------

library(tidypmc)
library(tidyverse)
library(rentrez)
library(ComplexUpset)



# Custom functions --------------------------------------------------------

#' metadata_as_tibble(id)

#' Description: Takes the Pubmed Central id of a paper and produce a tibble with its metadata
#' 
#' Input: 
#' id: character vector of length 1 in the shape of "PMC[0-9]{7}"
#' 
#' Output: 
#' a tibble with the metadata associated to the paper annotated under the id provided

metadata_as_tibble <- function(id){
  
  # Get the xml file of the paper
  pmc_xml(id = id) %>% 
    # Get the metadata of the paper
    pmc_metadata() %>% 
    # Transform it into a tidy format
    as_tibble()

}



#' id_to_tibble(id)

#' Description: Takes the Pubmed Central id of a paper and produce a tibble with its text and metadata
#' 
#' Input: 
#' id: character vector of length 1 in the shape of "PMC[0-9]{7}"
#' 
#' Output: 
#' a tibble with the text associated to the paper annotated under the id provided. It also includes columns
#' with the metadata of the paper by wrapping the "metadata_as_tibble" function.

id_to_table <- function(id){
  
  # Get the xml file of the paper
  pmc_xml(id = id) %>% 
    # Transform it into a tidy format
    pmc_text() %>% 
    # Add the ID of the paper as a column to link the text with its metadata
    mutate(PMCID = id) %>% 
    # Binding together the columns with the text of the paper along with its metadata 
    inner_join(., metadata_as_tibble(id), by = "PMCID")
  
}





# Data retrieving ---------------------------------------------------------


# Getting the papers of interest's IDs
papers_entrez <- entrez_search(db = "pmc", 
                               term = "gut microbiota[ALL] AND 2014[PDAT] AND english[LANG]",
                               retmax = 99999)

#"gut microbiota[ALL] AND 2012:2021[PDAT] AND english[LANG]",


# Adding PMC at the beggining suits the format of functions in 'tidypmc' package
papers_ids <- paste0("PMC", papers_entrez$ids) 


# Generate the dataset for 1 paper
microbiome_papers_text_and_metadata <- map_dfr(papers_ids, possibly(.f = id_to_table, otherwise = NULL))


# Exporting the data to csv files
# write_csv(x = microbiome_papers_text_and_metadata, file = "data/microbiome_papers_text_and_metadata.csv") # Server





# Data handling -----------------------------------------------------------


# Strings used to further filter papers
regex_16s <- "16[Ss]"
regex_wgs <- "[Ww]hole.{1}[Gg]enome.{1}[Ss]equencing|WGS|wgs|[Mm]etagenomics|[Mm]etagenomic"


beta_diversity_regex <- "[Bb]eta.{1}[Dd]iversity|β.diversity"
dim_reduction_regex <- "PCA|pca|PCOA|PCoA|pcoa|MDS|mds|NMDS|nmds"


bc_regex <- "[Bb]ray.{1}[Cc]urtis|BRAY.{1}CURTIS"
jac_regex <- "[Jj]accard|JACCARD"
w_unifrac_regex <- "(?<![Uu]n)[Ww]eighted.{1}[Uu]ni[Ff]rac|[Ww]eighted[Uu]ni[Ff]rac"
uw_unifrac_regex <- "[Uu]nweighted.{1}[Uu]ni[Ff]rac|[Uu]nweighted[Uu]ni[Ff]rac"
aitchison_regex <- "[Aa]itchison"

euclidian_regex <- "[Eu]clidean"
clr_regex <- "clr|CLR|[Cc]entered.{1}[Ll]ogarithmic.{1}[Tt]ransformation|[Cc]entered.{1}[Ll]og.{1}[Tt]ransformation"




# Format the dataset

# microbiome_papers_text_and_metadata <- read.csv(file = "data/microbiome_papers_text_and_metadata.csv") # Server
# microbiome_papers_text_and_metadata <- read.csv(file = "microbiome_literature_mining/data/microbiome_papers_text_and_metadata.csv") # PC

formatted_papers_dataset <- microbiome_papers_text_and_metadata %>% 
  
  mutate(across(.cols = where(is.character), .fns = replace_na, ""),
         major_section = case_when(str_detect(string = section, pattern = "[Aa]bstract") ~ "abstract",
                                   str_detect(string = section, pattern = "[Mm]ethods") ~ "methods",
                                   str_detect(string = section, pattern = "[Rr]esults") ~ "results",
                                   TRUE ~ "other")) %>% 
  
  filter(str_detect(major_section, "abstract|methods|results")) %>% 
  
  group_by(PMCID, major_section) %>% 
  
  # Mutate instead of summarize because I wanted to keep the rest of the columns
  mutate(full_text = str_c(text, collapse = " ")) %>% 
  
  ungroup() %>% 
  
  # Here I'm removing the paragraph and sentence columns and relocating the other columns
  select(major_section, full_text, PMCID:Date.received) %>% 
  
  # distinct removes the repeated rows generated with mutate str_c. 
  # This process is suboptimal. Think on more computationally effective ways to do the same process.
  distinct() %>% 
  
  pivot_wider(names_from = major_section, 
              values_from = full_text, 
              values_fill = "") %>%

  
  # Check if the articles are reviews or metaanalysis
  mutate(review = ifelse(str_detect(Title, "[Rr]eview") | str_detect(abstract, "[Rr]eview") | str_detect(methods, "[Rr]eview"), TRUE, FALSE),
         metaanalysis = ifelse(str_detect(Title, "[Mm]eta.{1}[Aa]nalysis") | str_detect(abstract, "[Mm]eta.{1}[Aa]nalysis") | str_detect(methods, "[Rr]eview"), TRUE, FALSE),
         original_article = ifelse(!review & !metaanalysis, TRUE, FALSE)) %>% 
  
  
  # Get decomposed information about the content of the paper
  mutate(
    # Check in the methods which technology was used
    methods_16s = str_detect(methods, regex_16s),
    methods_wgs = str_detect(methods, regex_wgs),

    # Checking beta div or dimmensionality reduction in methods
    methods_beta_div = str_detect(methods, beta_diversity_regex),
    methods_dim_reduction = str_detect(methods, dim_reduction_regex),
    
    # Check distance matrix in methods
    methods_bray_curtis = str_detect(methods, bc_regex),
    methods_jaccard = str_detect(methods, jac_regex), 
    methods_weighted_unifrac = str_detect(methods, w_unifrac_regex),
    methods_unweighted_unifrac = str_detect(methods, uw_unifrac_regex),
    methods_aitchison = str_detect(methods, aitchison_regex),
    methods_euclidian = str_detect(methods, euclidian_regex),
    methods_clr = str_detect(methods, clr_regex),
    
    
    
    # Checking beta div or dimmensionality reduction in results
    results_beta_div = str_detect(results, beta_diversity_regex),
    results_dim_reduction = str_detect(results, dim_reduction_regex),
    
    # Check distance matrix in results
    results_bray_curtis = str_detect(results, bc_regex),
    results_jaccard = str_detect(results, jac_regex),
    results_weighted_unifrac = str_detect(results, w_unifrac_regex),
    results_unweighted_unifrac = str_detect(results, uw_unifrac_regex),
    results_aitchison = str_detect(results, aitchison_regex),
    results_euclidian = str_detect(results, euclidian_regex),
    results_clr = str_detect(results, clr_regex),
    
    
    # Aggregating distances matrices mentioned in either methods or results
    braycurtis_performed = ifelse(methods_bray_curtis | results_bray_curtis, TRUE, FALSE),
    jaccard_performed = ifelse(methods_jaccard | results_jaccard, TRUE, FALSE),
    weighted_unifrac_performed = ifelse(methods_weighted_unifrac | results_weighted_unifrac, TRUE, FALSE),
    unweighted_unifrac_performed = ifelse(methods_unweighted_unifrac | results_unweighted_unifrac, TRUE, FALSE),
    #' Aitchison is considered as performed if:
    #' Either it was detected in the methods or results, or in case clr and euclidian are detected either in the methods or results.
    aitchison_performed = ifelse(methods_aitchison | results_aitchison | (methods_clr & methods_euclidian) | (results_clr & results_euclidian), TRUE, FALSE),
    #' Euclidian is considered as performed if:
    #' it was mentioned either in methods or results, but in either case, neither aitchison nor clr was mentioned neither in methods nor results.
    euclidian_performed = ifelse((methods_euclidian & !(methods_aitchison | methods_clr | results_aitchison | results_clr)) | (results_euclidian & !(results_aitchison | results_clr | methods_aitchison | methods_clr)), TRUE, FALSE))
  
  
# write_csv(x = formatted_papers_dataset, file = "data/formatted_papers_dataset.csv") # Server
# write_csv(x = formatted_papers_dataset, file = "microbiome_literature_mining/data/formatted_papers_dataset.csv") # PC




# Plots -------------------------------------------------------------------

# Reading the formatted papers dataset
# formatted_papers_dataset <- read_csv(file = "microbiome_literature_mining/data/formatted_papers_dataset.csv")



# How many papers are we dealing with
# Make this a table
formatted_papers_dataset %>% 
  drop_na(Year) %>% 
  count(Year)



# How many papers per year of each type are we dealing with

# DATA
n_of_papers_through_years_data <- formatted_papers_dataset %>% 
  as_tibble() %>% 
  pivot_longer(cols = c(review, metaanalysis, original_article),
               values_to = "value",
               names_to = "type_of_publication") %>% 
  drop_na() %>% 
  group_by(Year, type_of_publication) %>% 
  summarize(n_of_articles = sum(value), .groups = "drop")

# write_csv(x = n_of_papers_through_years_data, file = "microbiome_literature_mining/data/n_of_papers_through_years_data.csv") # PC 


# PLOT
n_of_papers_through_years_data %>% 
  
  ggplot(aes(x = Year, y = n_of_articles, color = type_of_publication)) + 
  geom_line(size = 1.15) + 
  geom_point(size = 2) +
  
  labs(x = "Year",
       y = "Number of articles published") +
  
  
  # CHANGE MIN(...$YEARS) FOR A NAMED VARIABLE. Something like: "years_with_microbiome_publications" or something alike.
  scale_x_continuous(breaks = seq(min(n_of_papers_through_years_data[n_of_papers_through_years_data$n_of_articles>0 , "Year"]),
                                  max(n_of_papers_through_years_data[n_of_papers_through_years_data$n_of_articles>0 , "Year"]),
                                  by = 1), 
                     limits = c(min(n_of_papers_through_years_data[n_of_papers_through_years_data$n_of_articles>0 , "Year"]),
                                max(n_of_papers_through_years_data[n_of_papers_through_years_data$n_of_articles>0 , "Year"])), 
                     expand = c(0.01,0.01)) + 
  
  scale_y_continuous(breaks = seq(min(n_of_papers_through_years_data$n_of_articles),
                                  max(n_of_papers_through_years_data$n_of_articles), 
                                  by = 100), 
                     limits = c(min(n_of_papers_through_years_data$n_of_articles),
                                max(n_of_papers_through_years_data$n_of_articles))) +
  
  guides(color = guide_legend("Type of publication")) +
  scale_color_discrete(breaks = c("original_article","review", "metaanalysis"),
                       labels = c("Research article", "Review", "Meta analysis")) +
  
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))

  
# ggsave(filename = "microbiome_literature_mining/figures/n_of_papers_through_years_data.jpg",
#        plot = last_plot(),
#        units = "in", height = 8, width = 10) # PC







# Distance matrices used

# DATA
stacked_barplot_data <- formatted_papers_dataset %>%
  
  group_by(Year) %>% 
  summarize(cum_bray = sum(braycurtis_performed, na.rm = T),
            cum_jaccard = sum(jaccard_performed, na.rm = T),
            cum_w_unifrac = sum(weighted_unifrac_performed, na.rm = T),
            cum_uw_unifrac = sum(unweighted_unifrac_performed, na.rm = T),
            cum_aitchison = sum(aitchison_performed, na.rm = T),
            cum_euclidian = sum(euclidian_performed, na.rm = T),
            .groups = "drop") %>% 

  pivot_longer(cols = starts_with("cum"),
               values_to = "num_studies",
               names_to = "distance_matrix") %>% 
  
  filter(num_studies > 0)

# write_csv(x = n_of_papers_through_years_data, file = "microbiome_literature_mining/data/stacked_barplot_distance_matrix_used_through_years.csv") # PC



# PLOT
stacked_barplot_data %>% 
  
  ggplot(aes(x = Year, y = num_studies, fill = distance_matrix)) + 
  geom_col(position = "fill") + 
  
  # Breaks cannot be hardcoded in automatic version
  scale_x_continuous(breaks = c(2012:2022), expand = c(0,0)) + 
  scale_y_continuous(labels = paste0(seq(0,100, by = 25),"%"), expand = c(0,0)) + 
  
  labs(y = "",
       x = "Year of publication") + 
  guides(fill = guide_legend("Distance\nmatrix\nused")) +
  
  scale_color_manual(values = "black") +
  scale_fill_brewer(palette = "Dark2",
                    labels = paste(c("Aitchison", 
                                     "Bray-Curtis", 
                                     "Euclidian",
                                     "Jaccard",
                                     "Unweighted\nUnifrac",
                                     "Weighted\nUnifrac"), "distance")) +
  
  theme_bw() + 
  
  theme(panel.grid = element_blank(),
        legend.title = element_text(face = "bold", size = 16),
        legend.position = "bottom",
        legend.text = element_text(size = 12),
        
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 14))


# ggsave(filename = "microbiome_literature_mining/figures/stacked_barplot_distance_matrix_used_through_years.jpg",
#        plot = last_plot(),
#        units = "in", height = 8, width = 10) # PC






# CLR through the years

# DATA
aitchison_through_years_data <- formatted_papers_dataset %>% 
  drop_na(Year, aitchison_performed) %>% 
  select(PMCID, Title, aitchison_performed, Year) %>% 
  group_by(Year) %>% 
  summarize(cum_aitchison_performed = sum(aitchison_performed), .groups = "drop")
  
# write_csv(x = n_of_papers_through_years_data, file = "microbiome_literature_mining/data/aitchison_through_years_data.csv") # PC



# PLOT
aitchison_through_years_data %>% 
  ggplot(aes(x = Year, y = cum_aitchison_performed)) + 
  
  geom_vline(xintercept = 2012, color = "grey60") + # HMP
  geom_vline(xintercept = 2017, color = "grey60") + # Compositional
  geom_vline(xintercept = 2019, color = "grey60") + # Neuroactive potential
  
  geom_line(size = 1.2) + 
  geom_point(size = 2) +
  
  scale_x_continuous(breaks = seq(2012, 2022, by = 1)) +
  
  labs(y = "",
       x = "Year of publication") +
  
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
  
# ggsave(filename = "microbiome_literature_mining/figures/line_plot_aitchison_through_years_data.jpg",
#        plot = last_plot(),
#        units = "in", height = 8, width = 10) # PC




# UPSETER PLOT
library(ggupset)




# Playground --------------------------------------------------------------

toy_papers_ids <- c("PMC3959530", "PMC3877837", "PMC4073011", "PMC3904282", "PMC3877837", "PMC4073018")

pmc_xml(toy_papers_ids[1])

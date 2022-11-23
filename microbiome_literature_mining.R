
#' TO DO LIST:
#' 
#' 3) Do the inner join in another function that takes as arguments the 2 custom functions made before
#' 
#' 4) Error handling of the metadata retrieving function
#' 
#' 5) Papers of interest must have 16S or WGS in their methods or results. Filter those that don't
#' 
#' 6) Interesting things the see:
#' 
#' 6.1 From the total of papers, how many use the different distance matrices (there are some euclidean without clr)
#' 
#' 6.2 Tendency over the years to use clr transformation
#'   6.2.1 Include : 
#'        https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224/full
#'        https://www.nature.com/articles/s41564-018-0337-x 
#'        https://pubmed.ncbi.nlm.nih.gov/31194939/ 
#'
#' 6.3 Alluvial plots
#' 
#' 6.4 https://cran.r-project.org/web/packages/ggupset/readme/README.html
#' 
#' 8) Check which studies are from human samples and which are from other animals samples




# Setting everything for the project --------------------------------------

library(tidypmc)
library(tidyverse)
library(rentrez)


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

  #print(paste("Retrieveing metadata of : ", id))
  
  # Get the xml file of the paper
  pmc_xml(id = id) %>% 
    # Get the metadata of the paper
    pmc_metadata() %>% 
    # Transform it into a tidy format
    as_tibble()
  
# print(paste("Metadata of : ", id, " succesfully retrieved"))
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
  
  print(paste("Text of : ", id, " succesfully retrieved"))
  
}



#' bool_str_detect (my_string, my_pattern)
#' 
#' Description: Takes a string and a pattern. If the pattern is in the string, 
#' the output is 1 (numeric), otherwise is 0 (numeric)
#' 
#' Input: 
#' my_string: a character vector of length 1
#' my_pattern: a character vector of length 1 or a regex expression
#' 
#' Output: a numeric vector of lenght 1, either 1 or 0
#' 
bool_str_detect <- function(my_string, my_pattern){
  
  ifelse(str_detect(string = my_string, my_pattern), yes = 1, no = 0)
  
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
#write_csv(x = microbiome_papers_text_and_metadata, file = "microbiome_papers_text_and_metadata.csv")





# Data handling -----------------------------------------------------------


# Strings used to further filter papers
strings_16s_regex <- "16[Ss]"
strings_wgs_regex <- "[Ww]hole.{1}[Gg]enome.{1}[Ss]equencing|WGS|wgs|[Mm]etagenomics|[Mm]etagenomic"


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
#microbiome_papers_text_and_metadata <- read.csv(file = "microbiome_literature_mining/data/microbiome_papers_text_and_metadata.csv")

formatted_papers_dataset <- microbiome_papers_text_and_metadata %>% 
  
  mutate(major_section = case_when(str_detect(string = section, pattern = "[Aa]bstract") ~ "abstract",
                                   str_detect(string = section, pattern = "[Mm]ethods") ~ "methods",
                                   str_detect(string = section, pattern = "[Rr]esults") ~ "results")) %>% 
  
  filter(str_detect(major_section, "abstract|methods|results")) %>% 
  
  group_by(PMCID, major_section) %>% 
  
  # Mutate instead of group_by because I wanted to keep the rest of the columns
  mutate(full_text = str_c(text, collapse = " ")) %>% 
  
  ungroup() %>% 
  
  # Here I'm removing the paragraph and sentence columns and relocating the other columns
  select(major_section, full_text, PMCID:Publisher) %>% 
  
  unique() %>% 
  
  pivot_wider(names_from = major_section, 
              values_from = full_text) %>% 

  
  # Check if the articles are reviews or metaanalysis
  mutate(review = ifelse(str_detect(Title, "[Rr]eview") | str_detect(abstract, "[Rr]eview"), 1, 0),
         metaanalysis = ifelse(str_detect(Title, "[Mm]eta.{1}[Aa]nalysis") | str_detect(abstract, "[Mm]eta.{1}[Aa]nalysis"), 1, 0),
         original_article = ifelse(!review & !metaanalysis, 1, 0)) %>% 
  
  
  
  # Get decomposed information about the content of the paper
  mutate(
    # Check in the methods which technology was used
    methods_16s = bool_str_detect(methods, strings_16s_regex),
    methods_wgs = bool_str_detect(methods, strings_wgs_regex),
    
    # Checking beta div or dimmensionality reduction in methods
    methods_beta_div = bool_str_detect(methods, beta_diversity_regex),
    methods_dim_reduction = bool_str_detect(methods, dim_reduction_regex),
    
    # Check distance matrix in methods
    methods_bray_curtis = bool_str_detect(methods, bc_regex),
    methods_jaccard = bool_str_detect(methods, jac_regex), 
    methods_weighted_unifrac = bool_str_detect(methods, w_unifrac_regex),
    methods_unweighted_unifrac = bool_str_detect(methods, uw_unifrac_regex),
    methods_aitchison = bool_str_detect(methods, aitchison_regex),
    methods_euclidian = bool_str_detect(methods, euclidian_regex),
    methods_clr = bool_str_detect(methods, clr_regex),
    
    
    
    # Checking beta div or dimmensionality reduction in results
    results_beta_div = bool_str_detect(results, beta_diversity_regex),
    results_dim_reduction = bool_str_detect(results, dim_reduction_regex),
    
    # Check distance matrix in results
    results_bray_curtis = bool_str_detect(results, bc_regex),
    results_jaccard = bool_str_detect(results, jac_regex),
    results_weighted_unifrac = bool_str_detect(results, w_unifrac_regex),
    results_unweighted_unifrac = bool_str_detect(results, uw_unifrac_regex),
    results_aitchison = bool_str_detect(results, aitchison_regex),
    results_euclidian = bool_str_detect(results, euclidian_regex),
    results_clr = bool_str_detect(results, clr_regex),
    
    
    # Aggregating distances matrices mentioned in either methods or results
    braycurtis_performed = ifelse(methods_bray_curtis | results_bray_curtis, 1, 0),
    jaccard_performed = ifelse(methods_jaccard | results_jaccard, 1, 0),
    weighted_unifrac_performed = ifelse(methods_weighted_unifrac | results_weighted_unifrac, 1, 0),
    unweighted_unifrac_performed = ifelse(methods_unweighted_unifrac | results_unweighted_unifrac, 1, 0),
    #' Aitchison is considered as performed if:
    #' Either it was detected in the methods or results, or in case clr and euclidian are detected either in the methods or results.
    aitchison_performed = ifelse(methods_aitchison | results_aitchison | (methods_clr & methods_euclidian) | (results_clr & results_euclidian), 1, 0),
    #' Euclidian is considered as performed if:
    #' it was mentioned either in methods or results, but in either case, neither aitchison nor clr was mentioned neither in methods nor results.
    euclidian_performed = ifelse((methods_euclidian & !(methods_aitchison | methods_clr | results_aitchison | results_clr)) | (results_euclidian & !(results_aitchison | results_clr | methods_aitchison | methods_clr)), 1, 0))
  
  
#write_csv(x = formatted_papers_dataset, file = "data/formatted_papers_dataset.csv")





# Plots -------------------------------------------------------------------

# How many papers are we dealing with
formatted_papers_dataset %>% 
  count(Year)


# How many papers per year of each type are we dealing with
n_of_papers_through_years_data <- formatted_papers_dataset_with_reviews %>% 
  as_tibble() %>% 
  pivot_longer(cols = c(review, metaanalysis_abstract, original_article),
               values_to = "value",
               names_to = "type_of_publication") %>% 
  group_by(Year, type_of_publication) %>% 
  summarize(n_of_articles = sum(value), .groups = "drop")
  
  
n_of_papers_through_years_data %>% 
  
  ggplot(aes(x = Year, y = n_of_articles, color = type_of_publication)) + 
  geom_line(size = 1.15) + 
  geom_point(size = 2) +
  
  labs(x = "Year",
       y = "Number of articles published") +
  
  scale_x_continuous(breaks = c(2012:2022), expand = c(0.01,0.01)) + 
  scale_y_continuous(breaks = seq(0,600, by = 100)) +
  
  guides(color = guide_legend("Type of publication")) +
  scale_color_discrete(breaks = c("original_article","review", "metaanalysis_abstract"),
                       labels = c("Research article", "Review", "Meta analysis")) +
  
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))
  



formatted_papers_dataset <- read.csv("microbiome_literature_mining/data/formatted_papers_dataset.csv")
formatted_papers_dataset_with_reviews <- read.csv("microbiome_literature_mining/data/formatted_papers_dataset_with_reviews.csv")


# Distance matrices used
stacked_barplot_data <- formatted_papers_dataset %>%
  group_by(Year) %>% 
  summarize(cum_bray = sum(braycurtis_performed),
            cum_jaccard = sum(jaccard_performed),
            cum_w_unifrac = sum(weighted_unifrac_performed),
            cum_uw_unifrac = sum(unweighted_unifrac_performed),
            cum_aitchison = sum(aitchison_performed),
            cum_euclidian = sum(euclidian_performed),
            .groups = "drop") %>% 

  pivot_longer(cols = starts_with("cum"),
               values_to = "num_studies",
               names_to = "distance_matrix") %>% 
  
  filter(num_studies > 0)

#write.csv(x = stacked_barplot_data, file = "data/stacked_barplot_data.csv")


# Plot the figure
stacked_barplot_data %>% 
  
  ggplot(aes(x = Year, y = num_studies, fill = distance_matrix)) + 
  geom_col(position = "fill") + 
  
  # Breaks cannot be hardcoded in automatic version
  scale_x_continuous(breaks = c(2012:2022), expand = c(0,0)) + 
  scale_y_continuous(labels = paste0(seq(0,100, by = 25),"%"), expand = c(0,0)) + 
  
  labs(y = "",
       x = "Year of publication") + 
  guides(fill = guide_legend("Distance matrix used")) +
  
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


ggsave(filename = "microbiome_literature_mining/figures/stacked_barplot_distance_matrix_used_through_years.jpg", 
       plot = last_plot(), 
       units = "in", height = 8, width = 10)



# CLR through the years

clr_through_years_data <- formatted_papers_dataset %>% 
  # Not sure about this : Maybe there is value in considering papers just with abstract and methods
  #filter(methods_16s == 1 | methods_wgs == 1) %>% 
  mutate(clr_performed = ifelse(methods_clr + results_clr >= 1, 1, 0)) %>% 
  select(PMCID, Title, clr_performed, Year) %>% 
  group_by(Year) %>% 
  summarize(cum_clr_performed = sum(clr_performed), .groups = "drop")
  

#write.csv(x = clr_through_years_data, file = "data/clr_through_years_data.csv")


clr_through_years_data %>% 
  ggplot(aes(x = Year, y = cum_clr_performed)) + 
  geom_line(size = 1.2) + 
  geom_point(size = 2) +
  
  geom_vline(xintercept = 2012) + # HMP
  geom_vline(xintercept = 2017) + # Compositional
  geom_vline(xintercept = 2019) + # Neuroactive potential
  
  scale_x_continuous(breaks = seq(2012, 2022, by = 1)) +
  
  labs(y = "",
       x = "Year of publication") +
  
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
  






# Playground --------------------------------------------------------------

toy_papers_ids <- paste0("PMC", c(3877837, 4073011, 3959530, 4428553, 4610029))
toy_papers_ids <- paste0("PMC", c(4428553, 4610029))
toy_microbiome_papers_text_and_metadata <- map_dfr(toy_papers_ids, possibly(.f = id_to_table, otherwise = NA))

toy_microbiome_papers_text_and_metadata %>% tibble() %>% filter(section == "Title") %>% pull(PMCID)
  pull(Title) %>% unique()

  
total_dataset <- read.csv(file = "microbiome_literature_mining/data/microbiome_papers_text_and_metadata.csv")

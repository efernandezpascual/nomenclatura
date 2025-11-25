
# packages_needed <- c(
#   "dplyr", "tibble", "stringr", "tidyr", "purrr",
#   "readxl", "openxlsx",
#   "TNRS", "flora", "rgbif",
#   "ggplot2", "scales",
#   "curl", "httr", "jsonlite"
# )
# 
# packages_to_install <- packages_needed[!(packages_needed %in% installed.packages()[,"Package"])]
# 
# if(length(packages_to_install) > 0){
#   install.packages(packages_to_install)
# }



library(dplyr)
library(stringr)



checklist <- read.csv("data/lista-campo.csv", sep = ";", stringsAsFactors = FALSE)


names(checklist) #Here we work with a column named "Species_name"


clean_name <- function(x) {
  x %>%
    str_squish() %>%
    str_replace_all("\\s+", " ") %>%
    str_replace_all("(?i)\\b(cf|aff|nr|sp|spp|group|gr)\\.?\\b", "") %>%
    str_replace_all("[?]", "") %>%
    str_squish()
} # Clean names of common things such as sp., etc.


####tnrs_powo###########
library(TNRS)

standardize_names_powo <- function(x,
                                   source = "wcvp",
                                   min_score_auto = 0.9,
                                   min_score_keep = 0.8) {
  
  df_in <- tibble(
    original_name = x,                # Store original names
    cleaned_name = clean_name(x)      # Preprocess names for TNRS
  )
  
  tnrs_raw <- TNRS(df_in$cleaned_name, sources = source)   # Query TNRS (POWO/WCVP)
  tnrs_raw <- as_tibble(tnrs_raw)                          # Convert to tibble
  
  message("Colonne TNRS: ", paste(names(tnrs_raw), collapse = ", "))  # Print TNRS columns
  
  best <- tnrs_raw %>%
    dplyr::group_by(Name_submitted) %>%                     # Group results by submitted name
    dplyr::slice_max(order_by = Overall_score, n = 1, with_ties = FALSE) %>% # Keep top-scoring match
    dplyr::ungroup()
  
  df_in %>%
    left_join(best, by = c("cleaned_name" = "Name_submitted")) %>%  # Join best match to input names
    transmute(
      original_name,                                   # Original input name
      cleaned_name,                                    # Cleaned version used for matching
      matched_name       = Name_matched,               # Matched taxon name
      accepted_name      = Accepted_name,              # Accepted name (if synonym)
      accepted_author    = Accepted_name_author,       # Author of accepted name
      family             = Name_matched_accepted_family, # Family according to POWO/WCVP
      powo_uri           = Name_matched_url,           # URL for POWO taxon page
      taxonomic_status   = Taxonomic_status,           # Accepted / Synonym / etc.
      overall_score      = Overall_score,              # TNRS matching score
      
      # Assign match category based on score and taxonomic info
      match_type = case_when(
        is.na(Overall_score) ~ "unresolved_no_match",        # No TNRS match
        Overall_score < min_score_keep ~ "unresolved_low_score", # Match too weak to trust
        Taxonomic_status == "Accepted" ~ "accepted_or_direct",   # Direct accepted match
        Taxonomic_status == "Synonym" &
          !is.na(Accepted_name) ~ "synonym_to_accepted",     # Synonym linked to an accepted name
        Overall_score >= min_score_auto ~ "fuzzy_auto",       # High-score fuzzy match
        TRUE ~ "fuzzy_needs_review"                           # Lower-score fuzzy match needing review
      )
    )
}





######rgbif################
library(rgbif)
library(dplyr)
library(tibble)

standardize_names_rgbif <- function(x,
                                    min_conf_auto = 90,
                                    min_conf_keep = 70) {
  
  df_in <- tibble(
    original_name = x,               # Store original names
    cleaned_name  = clean_name(x)    # Clean names for matching
  )
  
  # Query GBIF Backbone via rgbif
  bb_raw <- rgbif::name_backbone_checklist(df_in$cleaned_name)
  
  # Detect which column contains the submitted name
  verb_col <- if ("verbatim_name" %in% names(bb_raw)) {
    "verbatim_name"                   # Case 1: verbatim_name exists
  } else if ("verbatim_scientificName" %in% names(bb_raw)) {
    "verbatim_scientificName"         # Case 2: alternative column name
  } else {
    stop("not able to find taxon or scientific name")  # Unexpected structure
  }
  
  bb_raw$verbatim <- bb_raw[[verb_col]]   # Create unified submitted-name column
  
  bb_best <- bb_raw[!duplicated(bb_raw$verbatim), ]   # Keep first match per name
  
  df_in %>%
    left_join(bb_best, by = c("cleaned_name" = "verbatim")) %>%  # Join GBIF results
    transmute(
      original_name,                          # Original input name
      cleaned_name,                           # Cleaned version
      matched_name   = scientificName,        # GBIF matched scientific name
      canonical_name  = canonicalName,         # Canonical form (not always true accepted name)
      family         = family,                # Family from GBIF Backbone
      status         = status,                # Taxonomic status from GBIF
      matchType      = matchType,             # GBIF match category
      confidence     = confidence,            # GBIF confidence score
      
      # Assign simplified match classes
      match_type = case_when(
        is.na(matchType) | matchType == "NONE" ~ "no_match",    # No usable match
        confidence <  min_conf_keep           ~ "low_score",    # Below minimum threshold
        confidence >= min_conf_auto           ~ "auto",         # High-confidence auto match
        TRUE                                  ~ "review"        # Intermediate: needs review
      )
    )
}


#######checklist_creation####################
checklist <- checklist %>%
  mutate(cleaned_name = clean_name(Species_name))

unique_species_test <- unique(checklist$cleaned_name)

tnrs_table_test    <- standardize_names_powo(unique_species_test)
flora_table_test   <- standardize_names_flora(unique_species_test)
rgbif_table_test   <- standardize_names_rgbif(unique_species_test)

tnrs_table_test <- tnrs_table_test %>%
  rename_with(~ paste0("tnrs_", .x), -cleaned_name)

rgbif_table_test <- rgbif_table_test %>%
  rename_with(~ paste0("gbif_", .x), -cleaned_name)

flora_table_test <- flora_table_test %>%
  rename_with(~ paste0("flora_", .x), -cleaned_name)


checklist_norm <- checklist %>%
     left_join(tnrs_table_test, by = c("cleaned_name" = "cleaned_name"))

checklist_norm <- checklist_norm %>%
     left_join(rgbif_table_test, by = c("cleaned_name" = "cleaned_name"))

checklist_norm <- checklist_norm %>%
     left_join(flora_table_test, by = c("cleaned_name" = "cleaned_name"))


########score####

comparison <- checklist %>%
  select(cleaned_name) %>%                 # Keep only cleaned names
  distinct() %>%                           # Remove duplicates
  left_join(
    tnrs_table_test %>%                    # Join TNRS results
      select(cleaned_name,
             tnrs_name = tnrs_accepted_name,
             tnrs_overall_score = tnrs_overall_score,
             tnrs_match_type = tnrs_match_type),
    by = "cleaned_name"
  ) %>%
  left_join(
    rgbif_table_test %>%                   # Join GBIF results
      select(cleaned_name,
             gbif_name = gbif_matched_name,
             gbif_confidence = gbif_confidence,
             gbif_matchType = gbif_matchType),
    by = "cleaned_name"
  )



comparison <- comparison %>%
  mutate(
    gbif_score_norm = case_when(
      gbif_matchType == "NONE" ~ NA_real_,        # No match → no score
      gbif_matchType == "HIGHERRANK" ~ NA_real_,  # Matched only at higher rank → ignore
      TRUE ~ gbif_confidence / 100                # Normalize GBIF confidence to 0–1
    ),
    tnrs_score_norm = tnrs_overall_score          # TNRS is already on 0–1 scale
  )


scores_long <- comparison %>%
  select(cleaned_name,
         tnrs_score = tnrs_score_norm,
         gbif_score = gbif_score_norm) %>%
  tidyr::pivot_longer(
    cols = c(tnrs_score, gbif_score),             # Convert wide to long
    names_to = "service",                         # Column name indicating source
    values_to = "score"                           # Unified score column
  )


na_count <- comparison %>%
  summarise(
    tnrs_NA = sum(is.na(tnrs_score_norm)),        # Count TNRS missing scores
    gbif_NA = sum(is.na(gbif_score_norm)),        # Count GBIF missing scores
    total = n(),                                  # Total species
    tnrs_NA_pct = tnrs_NA / total,                # % TNRS missing
    gbif_NA_pct = gbif_NA / total                 # % GBIF missing
  )



#########visualization###########

library(dplyr)
library(ggplot2)
library(tidyr)

ggplot(scores_long, aes(service, score, fill = service)) +
  geom_boxplot() +
  labs(title = "Distribuzione degli score: TNRS vs GBIF",
       x = "Servizio",
       y = "Score") +
  theme_minimal()

ggplot(scores_long, aes(x = score, color = service)) +
  stat_ecdf(na.rm = TRUE) +
  labs(
    title = "Robustezza degli score (ECDF)",
    x = "Score (0–1)",
    y = "F(x)"
  ) +
  theme_minimal(base_size = 14)


###NA??

na_count

na_long <- comparison %>%
  summarise(
    TNRS = mean(is.na(tnrs_score_norm)),
    GBIF = mean(is.na(gbif_score_norm))
  ) %>%
  pivot_longer(cols = everything(),
               names_to = "service",
               values_to = "na_pct")


na_df <- na_long %>%
  mutate(label = paste0(round(na_pct*100, 1), "%"))

ggplot(na_df, aes(x = 2, y = na_pct, fill = service)) +
  geom_col(color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            color = "white",
            size = 6) +
  xlim(0.5, 2.5) +
  theme_void() +
  labs(title = "Percentage of NA")


##########dataset_export##########

#write.xlsx(checklist_norm, file = "Cleaned.xlsx", rowNames = FALSE)
write.csv(checklist_norm, file = "output/Campo_full.csv", row.names = FALSE)


###########other_list#####

checklist_GPT <-  read.csv("data/test_species_100.csv", sep = ";", stringsAsFactors = FALSE)
checklist_asturias <- read.csv("data/catalogo-asturias.csv", sep = ";", stringsAsFactors = FALSE)
checklist_iNaturalist <- read.csv("data/iNaturalist.csv", sep = ";", stringsAsFactors = FALSE)
checklist_Anna <- read.csv("data/Metadata_and_checklist.csv", sep = ";", stringsAsFactors = FALSE)

# Function to standardize species names in a dataframe using TNRS (POWO)
process_with_tnrs <- function(df, species_col, source = "wcvp") {
  
  # 1. Add a new column "cleaned_name" by cleaning the species names
  #    from the specified column using clean_name()
  df2 <- df %>%
    mutate(cleaned_name = clean_name(.data[[species_col]]))
  
  # 2. Extract unique cleaned names to avoid redundant queries
  unique_species <- unique(df2$cleaned_name)
  
  # 3. Standardize names via POWO (TNRS) to get accepted names, families, URIs, and scores
  tnrs_tab <- standardize_names_powo(unique_species, source = source)
  
  # 4. Prefix all TNRS output columns with "tnrs_" except for cleaned_name
  tnrs_prefixed <- tnrs_tab %>%
    rename_with(~ paste0("tnrs_", .x), -cleaned_name)
  
  # 5. Merge the TNRS results back to the original dataframe by cleaned_name
  df2 %>%
    left_join(tnrs_prefixed, by = "cleaned_name")
}

# List of checklists to process, each with its dataframe and species column name
checklists <- list(
  GPT        = list(df = checklist_GPT,        species_col = "Species_name"),
  Asturias   = list(df = checklist_asturias,   species_col = "Species_name"),
  iNaturalist= list(df = checklist_iNaturalist,species_col = "Species_name"),
  Anna       = list(df = checklist_Anna,       species_col = "Species_name")
)

# Create output directory if it doesn't exist
dir.create("output", showWarnings = FALSE)

# Process each checklist and save results as CSV
results_tnrs <- mapply(
  FUN = function(x, nm) {
    message("Processing: ", nm)  # Print which checklist is being processed
    
    # Run the TNRS processing function
    res <- process_with_tnrs(x$df, x$species_col)
    
    # Define output CSV file path
    out_file <- file.path("output", paste0(nm, "_TNRS.csv"))
    
    # Write the results to CSV with UTF-8 encoding
    write.csv(
      res,
      file = out_file,
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )
    
    # Return the processed dataframe
    res
  },
  x   = checklists,       # List of dataframes to process
  nm  = names(checklists),# Names of the checklists for output files
  SIMPLIFY = FALSE        # Keep results as a list of dataframes
)



######utaxonstand#######
#bonus
# github (requires `remotes` or `devtools`)
devtools::install_github("ecoinfor/U.Taxonstand")
library(U.Taxonstand)

load(file = "data/Taxon_Data/Plants_WCVP.rdata")

checklist %>%
  pull(cleaned_name) %>%
  unique -> sps


res <- nameMatch(spList = sps, spSource = database, author = TRUE, max.distance = 1) -> resolved

#####flora###############
library(flora)
library(dplyr)
library(purrr)
library(tibble)
library(stringr)

standardize_names_flora <- function(x) {
  
  df_in <- tibble(
    original_name = x,               # Store original names
    cleaned_name  = clean_name(x)    # Clean and normalize names
  )
  
  # Identify names that are non-empty and start with a letter
  valid_mask  <- !is.na(df_in$cleaned_name) &
    grepl("^[A-Za-z]", df_in$cleaned_name)
  
  # Unique valid names to query
  valid_names <- unique(df_in$cleaned_name[valid_mask])
  
  # If no valid names, return empty structured output
  if (length(valid_names) == 0) {
    return(
      df_in %>%
        mutate(
          matched_name      = NA_character_,   # No match
          accepted_name     = NA_character_,   # No accepted name
          family            = NA_character_,   # No family info
          taxonomic_status  = NA_character_    # No taxonomic status
        )
    )
  }
  
  # Safe wrapper around flora::get.taxa to avoid breaking on errors
  safe_get_taxa <- function(nm_vec) {
    tryCatch(
      flora::get.taxa(
        taxa             = nm_vec,
        replace.synonyms = TRUE,    # Replace synonyms with accepted names
        suggest.names    = TRUE     # Suggest close matches
      ),
      error = function(e) {
        warning("flora::get.taxa error: ", e$message)
        tibble(                             # Return fallback tibble
          taxon           = nm_vec,
          scientific.name = NA_character_,
          accepted.name   = NA_character_,
          family          = NA_character_,
          taxon.status    = NA_character_
        )
      }
    )
  }
  
  flora_raw <- safe_get_taxa(valid_names)    # Query Brazilian Flora
  
  # Detect which column contains the submitted name
  name_col <- dplyr::case_when(
    "taxon"           %in% names(flora_raw) ~ "taxon",
    "scientific.name" %in% names(flora_raw) ~ "scientific.name",
    TRUE ~ ""                               # No valid column found
  )
  
  if (name_col == "") {
    stop("not able to find taxon or scientific name")   # Fail if structure unexpected
  }
  
  flora_raw <- flora_raw %>%
    mutate(submitted = .data[[name_col]])   # Create a unified submitted-name column
  
  flora_best <- flora_raw[!duplicated(flora_raw$submitted), ]  # Keep first match per name
  
  # Join flora results back to the original input order
  df_out <- df_in %>%
    left_join(
      flora_best,
      by = c("cleaned_name" = "submitted")
    ) %>%
    transmute(
      original_name,                         # Original input name
      cleaned_name,                          # Normalized name used for matching
      matched_name     = scientific.name,    # Name matched by Flora
      accepted_name    = accepted.name,      # Accepted name according to Flora
      family           = family,             # Family information
      taxonomic_status = taxon.status        # Taxonomic status from Flora
    )
  
  df_out
}







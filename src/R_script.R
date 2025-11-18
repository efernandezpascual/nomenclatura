setwd("") #insert the path of the downloaded folder

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


library(readxl)
library(dplyr)
library(stringr)
library(openxlsx)


checklist <- read_excel("data/lista-campo.xlsx")


names(checklist) 


clean_name <- function(x) {
  x %>%
    str_squish() %>%
    str_replace_all("\\s+", " ") %>%
    # togli cf, aff, nr, sp, spp, group, gr (maiusc/minusc, con o senza punto)
    str_replace_all("(?i)\\b(cf|aff|nr|sp|spp|group|gr)\\.?\\b", "") %>%
    str_replace_all("[?]", "") %>%
    str_squish()
}


####tnrs_powo###########
library(TNRS)

standardize_names_powo <- function(x,
                                   source = "wcvp",
                                   min_score_auto = 0.9,
                                   min_score_keep = 0.8) {
  
  df_in <- tibble(
    original_name = x,
    cleaned_name = clean_name(x)
  )
  
  tnrs_raw <- TNRS(df_in$cleaned_name, sources = source)
  tnrs_raw <- as_tibble(tnrs_raw)   
  
  message("Colonne TNRS: ", paste(names(tnrs_raw), collapse = ", "))
  
  best <- tnrs_raw %>%
    dplyr::group_by(Name_submitted) %>%   
    dplyr::slice_max(order_by = Overall_score, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  df_in %>%
    left_join(best, by = c("cleaned_name" = "Name_submitted")) %>%
    transmute(
      original_name,
      cleaned_name,
      matched_name       = Name_matched,
      accepted_name      = Accepted_name,
      accepted_author    = Accepted_name_author,
      family             = Name_matched_accepted_family,
      powo_uri           = Name_matched_url,
      taxonomic_status   = Taxonomic_status,
      overall_score      = Overall_score,
      match_type = case_when(
        is.na(Overall_score) ~ "unresolved_no_match",
        Overall_score < min_score_keep ~ "unresolved_low_score",
        Taxonomic_status == "Accepted" ~ "accepted_or_direct",
        Taxonomic_status == "Synonym" &
          !is.na(Accepted_name) ~ "synonym_to_accepted",
        Overall_score >= min_score_auto ~ "fuzzy_auto",
        TRUE ~ "fuzzy_needs_review"
      )
    )
}



#####flora###############
library(flora)
library(dplyr)
library(purrr)
library(tibble)
library(stringr)

standardize_names_flora <- function(x) {
  
  df_in <- tibble(
    original_name = x,
    cleaned_name  = clean_name(x)
  )
  
  
  valid_mask  <- !is.na(df_in$cleaned_name) &
    grepl("^[A-Za-z]", df_in$cleaned_name)
  
  valid_names <- unique(df_in$cleaned_name[valid_mask])
  
  
  if (length(valid_names) == 0) {
    return(
      df_in %>%
        mutate(
          matched_name      = NA_character_,
          accepted_name     = NA_character_,
          family            = NA_character_,
          taxonomic_status  = NA_character_
        )
    )
  }
  
  
  safe_get_taxa <- function(nm_vec) {
    tryCatch(
      flora::get.taxa(
        taxa             = nm_vec,
        replace.synonyms = TRUE,
        suggest.names    = TRUE
      ),
      error = function(e) {
        warning("flora::get.taxa error: ", e$message)
        tibble(
          taxon          = nm_vec,
          scientific.name = NA_character_,
          accepted.name   = NA_character_,
          family          = NA_character_,
          taxon.status    = NA_character_
        )
      }
    )
  }
  
  flora_raw <- safe_get_taxa(valid_names)
  
  
  name_col <- dplyr::case_when(
    "taxon"           %in% names(flora_raw) ~ "taxon",
    "scientific.name" %in% names(flora_raw) ~ "scientific.name",
    TRUE ~ ""
  )
  
  if (name_col == "") {
    stop("not able to find taxon or scientific name")
  }
  
  flora_raw <- flora_raw %>%
    mutate(submitted = .data[[name_col]])
  
  
  flora_best <- flora_raw[!duplicated(flora_raw$submitted), ]
  
  df_out <- df_in %>%
    left_join(
      flora_best,
      by = c("cleaned_name" = "submitted")
    ) %>%
    transmute(
      original_name,
      cleaned_name,
      matched_name     = scientific.name,
      accepted_name    = accepted.name,
      family           = family,
      taxonomic_status = taxon.status
    )
  
  df_out
}

######rgbif################
library(rgbif)
library(dplyr)
library(tibble)

standardize_names_rgbif <- function(x,
                                    min_conf_auto = 90,
                                    min_conf_keep = 70) {
  
  df_in <- tibble(
    original_name = x,
    cleaned_name  = clean_name(x)
  )
  
  
  bb_raw <- rgbif::name_backbone_checklist(df_in$cleaned_name)
  

  verb_col <- if ("verbatim_name" %in% names(bb_raw)) {
    "verbatim_name"
  } else if ("verbatim_scientificName" %in% names(bb_raw)) {
    "verbatim_scientificName"
  } else {
    stop("not able to find taxon or scientific name")
  }
  

  bb_raw$verbatim <- bb_raw[[verb_col]]
  bb_best <- bb_raw[!duplicated(bb_raw$verbatim), ]
  
  df_in %>%
    left_join(bb_best, by = c("cleaned_name" = "verbatim")) %>%
    transmute(
      original_name,
      cleaned_name,
      matched_name   = scientificName,
      accepted_name  = canonicalName,
      family         = family,
      status         = status,
      matchType      = matchType,
      confidence     = confidence,
      match_type = case_when(
        is.na(matchType) | matchType == "NONE" ~ "no_match",
        confidence <  min_conf_keep           ~ "low_score",
        confidence >= min_conf_auto           ~ "auto",
        TRUE                                  ~ "review"
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
  select(cleaned_name) %>%
  distinct() %>%
  left_join(
    tnrs_table_test %>% 
      select(cleaned_name,
             tnrs_name = tnrs_accepted_name,
             tnrs_overall_score = tnrs_overall_score,
             tnrs_match_type = tnrs_match_type),
    by = "cleaned_name"
  ) %>%
  left_join(
    rgbif_table_test %>%
      select(cleaned_name,
             gbif_name = gbif_accepted_name,
             gbif_confidence = gbif_confidence,
             gbif_matchType = gbif_matchType),
    by = "cleaned_name"
  )



comparison <- comparison %>%
  mutate(
    gbif_score_norm = case_when(
      gbif_matchType == "NONE" ~ NA_real_,           # nessun match → score NA
      gbif_matchType == "HIGHERRANK" ~ NA_real_,     # match al genere/famiglia → NA
      TRUE ~ gbif_confidence / 100                   # match reale → normalizza
    ),
    tnrs_score_norm = tnrs_overall_score             # già 0–1
  )


scores_long <- comparison %>%
  select(cleaned_name,
         tnrs_score = tnrs_score_norm,
         gbif_score = gbif_score_norm) %>%
  pivot_longer(cols = c(tnrs_score, gbif_score),
               names_to = "service",
               values_to = "score")

na_count <- comparison %>%
  summarise(
    tnrs_NA = sum(is.na(tnrs_score_norm)),
    gbif_NA = sum(is.na(gbif_score_norm)),
    total = n(),
    tnrs_NA_pct = tnrs_NA / total,
    gbif_NA_pct = gbif_NA / total
  )



#########visualization###########

library(dplyr)
library(ggplot2)
library(tidyr)

# long format per ggplot
df_long <- comparison %>%
  select(cleaned_name, gold, tnrs, gbif) %>%
  pivot_longer(cols = c(tnrs, gbif),
               names_to = "service",
               values_to = "matched") %>%
  mutate(correct = matched == gold)

accuracy_table <- df_long %>%
  group_by(service) %>%
  summarise(accuracy = mean(correct, na.rm = TRUE))

ggplot(accuracy_table, aes(x = service, y = accuracy, fill = service)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Accuracy of TNRS vs GBIF",
       x = "Service",
       y = "Accuracy (%)") +
  theme_minimal()



unresolved_table <- df_long %>%
  group_by(service) %>%
  summarise(unresolved = mean(is.na(matched)))

ggplot(unresolved_table, aes(x = service, y = unresolved, fill = service)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "UNRESOLVED for service",
       x = "Service",
       y = "Unresolved (%)") +
  theme_minimal()


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

high_threshold <- 0.8

scores_high <- scores_long %>%
  group_by(service) %>%
  summarise(
    n_total = sum(!is.na(score)),
    n_high  = sum(score >= high_threshold, na.rm = TRUE),
    pct_high = n_high / n_total
  )

ggplot(scores_high, aes(x = service, y = pct_high, fill = service)) +
  geom_col(width = 0.6) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = paste0("Quota di score ≥ ", high_threshold),
    x = "Servizio",
    y = "Percentuale di score alti"
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

###########other_list#####

checklist_GPT <- read_excel("data/test_species_100_more_errors.xlsx")
checklist_asturias <- read_excel("data/catalogo-asturias.xlsx")
checklist_iNaturalist <- read_excel("data/iNaturalist.xlsx")
checklist_Anna <- read_excel("data/Metadata_and_checklist.xlsx")

process_with_tnrs <- function(df, species_col, source = "wcvp") {
  df2 <- df %>%
    mutate(cleaned_name = clean_name(.data[[species_col]]))
  
  unique_species <- unique(df2$cleaned_name)
  
  tnrs_tab <- standardize_names_powo(unique_species, source = source)
  
  tnrs_prefixed <- tnrs_tab %>%
    rename_with(~ paste0("tnrs_", .x), -cleaned_name)
  
  df2 %>%
    left_join(tnrs_prefixed, by = "cleaned_name")
}

checklists <- list(
  GPT        = list(df = checklist_GPT,        species_col = "Species_name"),
  Asturias   = list(df = checklist_asturias,   species_col = "Species_name"),
  iNaturalist= list(df = checklist_iNaturalist,species_col = "Species_name"),
  Anna       = list(df = checklist_Anna,       species_col = "Species_name")
)

#dir.create("output", showWarnings = FALSE)

results_tnrs <- mapply(
  FUN = function(x, nm) {
    message("Processing: ", nm)
    res <- process_with_tnrs(x$df, x$species_col)
    out_file <- file.path("output", paste0(nm, "_TNRS.xlsx"))
    openxlsx::write.xlsx(res, file = out_file, rowNames = FALSE)
    res
  },
  x   = checklists,
  nm  = names(checklists),
  SIMPLIFY = FALSE
)









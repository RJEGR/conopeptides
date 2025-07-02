

# Estimates the effect size of diet and time and superfamilies (optional) over DEGs
# as preliminar does not display proper effect size 

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)


pub_dir <- "/Users/cigom/Documents/GitHub/conopeptides/PUBLICATION_DIR/"

dir <- "/Users/cigom/Documents/GitHub/conopeptides/06.Quantification/MATRIX_RSEM_dir/"

# data 1 ----

subdir <- "KALLISTO_Merged_polyA_hisat_SuperDuper.fasta.transdecoder_DIR/"

f <- "KALLISTO_Merged_polyA_hisat_SuperDuper.fasta.transdecoder.matrix"

f <- list.files(file.path(dir, subdir), f, full.names = T)

dim(datExpr <- round(readRDS(f)))

recode_to <- structure(c("Control","Shrimp", "Mollusk", "Polychaete", "Mixed"))

recode_to <- structure(recode_to, names = c("Ctrl","Cam","Lit","Pol","Mix"))

recode_time <- structure(c("","2months", "4months"))
recode_time <- structure(recode_time, names = c("ctrl","2","4"))

# data 2 ----

.colData <- list.files(dir, pattern = "Manifest", full.names = T) 

Manifest <- read_tsv(.colData) %>% 
  mutate(LIBRARY_ID = ifelse(grepl("cam2v_", LIBRARY_ID), NA, LIBRARY_ID)) %>%
  mutate(LIBRARY_ID = ifelse(grepl("cam6", LIBRARY_ID), NA, LIBRARY_ID)) %>%
  # mutate(LIBRARY_ID = ifelse(grepl("ctrl", LIBRARY_ID), NA, LIBRARY_ID)) %>%
  select(LIBRARY_ID, Time, Diatery) %>% drop_na(LIBRARY_ID) %>%
  dplyr::mutate(Diatery = dplyr::recode_factor(Diatery, !!!recode_to)) %>%
  dplyr::mutate(Time = dplyr::recode_factor(Time, !!!recode_time)) %>%
  mutate_if(is.character, as.factor)
    
# data 3 ----

DEGS <- read_rds(file.path(pub_dir, "EDGER_EffectSize_input.rds")) %>% filter(abs(logFC) > 4 & FDR < 0.05) 


cols_to_check <- c("sampleA", "sampleB")

DEGS %>%
  filter(
    if_any(all_of(cols_to_check), ~ str_detect(.x, "Ctrl"))) %>%
  filter()
  count(sampleA, sampleB, sampleX, sort = T)


genes <- DEGS %>% distinct(protein_id, pep_seq) %>% pull(pep_seq, name = protein_id)

# genes <- structure(genes, names = genes)

sum(keep <- rownames(datExpr) %in% sort(unique(names(genes))))/length(genes)

dim(datExpr <- datExpr[keep,])

genes <- genes[match(rownames(datExpr), names(genes))]

identical(rownames(datExpr), names(genes))

rownames(datExpr) <- genes


DATA <- datExpr %>%  as_tibble(rownames = 'pep_seq') %>% 
  distinct() %>%
  pivot_longer(-pep_seq, values_to = "expression", names_to = "LIBRARY_ID") %>%
  right_join(Manifest) %>% filter(expression > 0) %>%
  left_join(DEGS %>% distinct(pep_seq, sf))
# mutate_if(is.character, as.factor)



library(tidymodels)

# Wrap ------

fit_model <- function(data, formula) {
  
  
  require(tidymodels)
  
  # formula <- formula("expression ~ EDAD + gene")
  
  # formula <- c("expression ~ CONTRASTE_A + EDAD + gene")
  
  formula <- formula(formula)
  
  data_rec <- recipe(formula, data = data) %>%
    step_normalize(all_numeric_predictors()) %>%
    step_dummy(all_nominal_predictors())
  
  
  
  model_spec <- linear_reg() %>% set_engine("lm")
  
  data_workflow <- workflow() %>%
    add_recipe(data_rec) %>%
    add_model(model_spec)
  
  # Split training and test model set
  
  set.seed(123)
  
  data_split <- initial_split(data, prop = 0.75)
  data_train <- training(data_split)
  data_test <- testing(data_split)
  
  data_fit <- data_workflow %>% fit(data = data_train)
  
  # data_predictions <- data_fit %>%
  #   predict(new_data = data_test) %>%
  #   bind_cols(data_test)
  
  results <- data_fit %>%
    extract_fit_parsnip() %>%
    tidy()
  
  
  return(results)
}

# vars <- c("Condition","Site","design")

vars <- c("Time","Diatery")

fit_model(DATA, formula = "expression ~ sf")

fit_data <- list()

for(i in vars) {
  
  form <- paste0("expression ~ ", i)
  
  cat("\nCalculating ", form,"\n")
  
  fit_data[[i]] <- fit_model(DATA, form) %>% mutate(Intercept = i)
  
}

fitdf <- do.call(rbind, fit_data)

fitdf <- fitdf %>%
  mutate(star = ifelse(p.value <.001, "***", 
    ifelse(p.value <.01, "**",
      ifelse(p.value <.05, "*", "ns")))) 

recode_to <- c( 
  `Time_X4months` = "4 months",
  `Time_Ctrl`= "Control",
  `Diatery_Shrimp`= "Shrimp",
  `Diatery_Mollusk`= "Mollusk",
  `Diatery_Polychaete`= "Polychaete",
  `Diatery_Mixed`= "Mixed")

fitdf <- fitdf %>%
  mutate(term = ifelse(Intercept == "Time" & grepl("Intercept", term), "2 months", term)) %>%
  mutate(term = ifelse(Intercept == "Diatery" & grepl("Intercept", term), "Control", term)) %>%
  dplyr::mutate(term = dplyr::recode_factor(term, !!!rev(recode_to))) %>%
  # separate(term, into = c("facet", "term"), sep = "_") %>%
  mutate(x_star = estimate + (0.1+std.error) * sign(estimate)) %>%
  mutate(xmin = estimate-std.error, xmax = estimate+std.error) %>%
  mutate(label = paste0(term, " (", star,")"))


fitdf %>% count(term, Intercept)

# Setting random valies

fit_data <- list()

counter <- 1

for (i in vars) {
  form <- paste0("expression ~ ", i)
  cat("\nCalculating ", form, "\n")
  
  for (n in 1:3) {
    # Get 5 random samples each time
    sample_data <- readRDS(f) %>%
      as_tibble(rownames = 'protein_id') %>%
      distinct() %>%
      dplyr::slice_sample(n = 5) %>%
      pivot_longer(-protein_id, values_to = "expression", names_to = "LIBRARY_ID") %>%
      right_join(Manifest) %>%
      filter(expression > 0)
    
    fit <- fit_model(sample_data, form) %>%
      mutate(Intercept = i, RandomSet = n)
    
    fit_data[[counter]] <- fit
    counter <- counter + 1
  }
}


fitRandomdf <- do.call(rbind, fit_data)

fitRandomdf <- fitRandomdf %>%
  mutate(star = ifelse(p.value <.001, "***", 
    ifelse(p.value <.01, "**",
      ifelse(p.value <.05, "*", "ns")))) %>%
  mutate(term = ifelse(Intercept == "Time" & grepl("Intercept", term), "2 months", term)) %>%
  mutate(term = ifelse(Intercept == "Diatery" & grepl("Intercept", term), "Control", term)) %>%
  dplyr::mutate(term = dplyr::recode_factor(term, !!!rev(recode_to))) %>%
  # separate(term, into = c("facet", "term"), sep = "_") %>%
  mutate(x_star = estimate + (0.1+std.error) * sign(estimate)) %>%
  mutate(xmin = estimate-std.error, xmax = estimate+std.error) %>%
  mutate(label = paste0(term, " (", star,")"))

fitRandomdf %>% count(Intercept, RandomSet)

fitdf <- fitdf %>% mutate(RandomSet = "True") %>%
  select(names(fitRandomdf)) %>%
  rbind(fitRandomdf)


fitdf %>%
  mutate(Intercept = ifelse(RandomSet != "True", "Random", Intercept)) %>%
  ggplot(aes(y = term, x = estimate, group = Intercept)) + # color = Intercept
  facet_grid(Intercept~ ., scales = "free", space = "free") +
  geom_point(size = 2, position = position_dodge(0.5)) +
  geom_text(aes(x = x_star, label = star),
    vjust = 0.5, hjust = 1, size= 4,
    color="black",
    position=position_dodge(0.5),
    family =  "GillSans") +
  geom_errorbar(aes(xmin = xmin, xmax = xmax),
    width = 0.1, alpha = 0.3, 
    # color = "black",
    position=position_dodge(width = 0.5)
  ) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5, color = "black") +
  labs(
    y = "",
    x = "Effect Size")


p <- fitdf %>%
  filter(RandomSet == "True") %>%
  ggplot(aes(y = Intercept, x = estimate, color = term, fill = term)) + # color = Intercept
  # facet_grid(Intercept~ ., scales = "free", space = "free") +
  geom_text(aes(x = x_star, label = label),
    vjust = 0.5, hjust = 1, size= 4,
    color="black",
    position=position_dodge(0.5),
    family =  "GillSans") +
  geom_point(size = 2, position = position_dodge(0.5)) +
  geom_errorbar(aes(xmin = xmin, xmax = xmax),
    width = 0.1, alpha = 0.3, 
    # color = "black",
    position=position_dodge(width = 0.5)
  ) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5, color = "black") +
  labs(
    y = "",
    x = "Effect Size") +
  scale_x_continuous(limits = c(-1000,700)) +
  theme_bw(base_family = "GillSans", base_size = 12)  +
  theme(legend.position = "none",
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    strip.text = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  )

p

ggsave(p, filename = 'MultipleContrast_effectSize.png', 
  path = pub_dir, width = 2.5, height = 3, device = png, dpi = 700)



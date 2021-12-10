
get_labs <- function(data){
  sapply(data, attr, which = "label")
}

find_labs <- function(data, label){
  labs <- get_labs(data)
  which(grepl(label, labs, ignore.case = TRUE))
}

search_for_symptoms_vars <- function(pr_data){
  key_terms <- c("weakness", "jaundice", "breathing", "yellow", "bleed", "consc", "symptoms")
  labs <- get_labs(pr_data)
  find <- stringr::str_detect(labs, paste(key_terms, collapse = "|"))
  
  if(length(find) > 0){
    print(labs[find])
    return(names(labs[find]))
  } else {
    message("No matched found")
  }
  return(NULL)
}

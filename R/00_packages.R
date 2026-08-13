
.required_packages <- c(
  "data.table",
  "readr",
  "dplyr",
  "tidyr",
  "purrr",
  "stringr",
  "lubridate",
  "tibble",
  "mice",
  "bestNormalize",
  "depmixS4",
  "zoo",
  "ggplot2",
  "DBI",
  "duckdb"
)


.missing_required <- .required_packages[
  !vapply(.required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(.missing_required) > 0L) {
  stop(
    "Missing required R packages: ",
    paste(.missing_required, collapse = ", "),
    "\nInstall them before running the analysis.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  for (.pkg in .required_packages) library(.pkg, character.only = TRUE)
})

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
mutate <- dplyr::mutate
summarise <- dplyr::summarise
summarize <- dplyr::summarise
arrange <- dplyr::arrange
group_by <- dplyr::group_by
ungroup <- dplyr::ungroup
left_join <- dplyr::left_join
inner_join <- dplyr::inner_join
right_join <- dplyr::right_join
full_join <- dplyr::full_join
anti_join <- dplyr::anti_join
semi_join <- dplyr::semi_join
bind_rows <- dplyr::bind_rows
all_of <- tidyselect::all_of
any_of <- tidyselect::any_of

options(datatable.integer64 = "double")
rm(.pkg, .missing_required, .required_packages)

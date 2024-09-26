### Function to format p-values as stars ###
# Modified by Bolívar Aponte Rolón (@jibarozzo) from @rich-iannone's fmt_star function posted
# on GitHub: https://github.com/rstudio/gt/issues/187#issuecomment-465769854


fmt_stars <- function(data,
                      columns,
                      rows = NULL) {
  
  # Capture expression in `rows`
  rows <- rlang::enquo(rows)
  
  # Pass `data`, `columns`, `rows`, and the formatting
  # functions as a function list to `fmt()`
  fmt(
    data = data,
    columns = columns,
    rows = !!rows,
    fns = list(
      default = function(x) {
        
        x_str <- 
          dplyr::case_when(
            dplyr::between(x, 0, 0.0001) ~ "****",
            dplyr::between(x, 0, 0.001) ~ "***",
            dplyr::between(x, 0, 0.01) ~ "**",
            dplyr::between(x, 0, 0.05) ~ "*",
            TRUE ~ "ns"
          )
      }
    )
  )
}
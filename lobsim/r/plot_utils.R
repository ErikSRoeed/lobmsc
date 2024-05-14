#---- PLOT UTILS ---------------------------------------------------------------
#>                                                                            <#
#>    Erik Sandertun Roeed                                                    <#
#>                                                                            <#
#>    Utility functions for plotting-related data wrangling, etc.             <#
#>                                                                            <#
#>    Contents                                                                <#
#>    function: replace_underscores()                                         <#
#>                                                                            <#
#------------------------------------------------------------------------------#

replace_underscores <- function(string, with = " ")
{
  string_without_underscores <- string |>
    stringr::str_replace_all(pattern = "_", replacement = with)
  return(string_without_underscores)
}

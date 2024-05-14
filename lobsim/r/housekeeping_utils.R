#---- PROJECT_UTILS ------------------------------------------------------------
#>                                                                            <#
#>    Erik Sandertun Roeed                                                    <#
#>                                                                            <#
#>    Contents                                                                <#
#>    function: check_slimmr_installed()                                      <#
#>                                                                            <#
#------------------------------------------------------------------------------#

check_slimmr_installed <- function()
{
  if (! "slimmr" %in% utils::installed.packages())
  {
    AFFIRMATIVE_USER_INPUTS <- c("YES", "Yes", "Y", "y")
    user_input <- readline("Install ErikSRoeed/slimmr from GitHub? [Y/N] ")
    
    if (! user_input %in% AFFIRMATIVE_USER_INPUTS)
    {
      return(FALSE)
    }
    
    devtools::install_github("ErikSRoeed/slimmr")
  }
  
  return(TRUE)
}

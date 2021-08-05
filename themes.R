# Custom theme for flextable
f_theme_surial <- function(ftab) {
  ftab %>% 
    font(fontname = "Arial Narrow", part = "all") %>% 
    fontsize(size = 10, part = "all") %>% 
    align(j = 2, align = "right", part = "all") %>% 
    height(height = 0, part = "body") %>% 
    autofit() %>% 
    padding(padding.top = 2, padding.bottom = 2)
}

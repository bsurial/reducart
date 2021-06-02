# Custom theme for flextable
f_theme_surial <- function(ftab) {
  ftab %>% 
    font(fontname = "Arial Narrow", part = "all") %>% 
    fontsize(size = 9, part = "body") %>% 
    fontsize(size = 8, part = "footer") %>% 
    align(j = 2, align = "right", part = "all") %>% 
    height(height = 0, part = "body") %>% 
    line_spacing(space = 0.5) %>% 
    line_spacing(space = 0.75, part = "footer") %>% 
    autofit()
}

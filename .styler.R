# Configuration for styler
# This configuration follows Bioconductor style guidelines

# Use tidyverse style as base
style_transformers <- styler::tidyverse_style(
  # Indentation
  indent_by = 2,
  
  # Spacing
  strict = TRUE,
  
  # Line breaks
  math_token_spacing = styler::specify_math_token_spacing(
    zero = c("'^'", "'/'", "'*'", "'+'"),
    one = c("'-'", "'='", "'~'", "'<'", "'>'", "'<='", "'>='", "'&'", "'|'", "'=='", "'!='")
  ),
  
  # Always use explicit returns
  reindention = styler::specify_reindention(
    regex_pattern = NULL,
    indention = 0,
    comments_only = TRUE
  ),
  
  # Comments
  scope = "tokens"
)

# Ensure each file ends with a newline
style_transformers$line_break$add_line_break_after_source <- TRUE

# Never use tabs, always use spaces
style_transformers$indention$indent_character <- " "

# Return the transformer
style_transformers
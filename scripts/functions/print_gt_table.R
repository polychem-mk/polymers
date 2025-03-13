# Converts data frame into a gt_tbl table 

# required packages:
# library(gt)

# arguments:
#          df    data frame;
# table_title    character, title for a table; 
# col_padding    numeric, spaces between column names.


print_gt_table <- function(df,
                           table_title = "", 
                           col_padding = 25 )
  {
  gt( df) %>% 
    tab_header(md(
      paste('<span style="color:#708090; font-family:system-ui; ">',
            table_title,'</span>'))) %>%
    tab_options(table.font.size = 12,
                table.font.names = "Helvetica",
                heading.title.font.size = 14,
                heading.align = "left",
                column_labels.padding.horizontal = col_padding,
                data_row.padding = 4,
                table.border.top.color = "white",
                table_body.hlines.color = "white",
                table.border.bottom.width = 0.6,
                column_labels.font.weight = "bold",
                column_labels.border.top.width = 0.6,
                column_labels.border.bottom.width = 0.6 )
}

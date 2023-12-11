library(grid)
library(forestploter)

dt <- read.csv("反向IVW.csv",header = T)

dt$Exposure <- ifelse(is.na(dt$or), 
                      dt$Exposure,
                      paste0("   ", dt$Exposure))

dt$Case <- ifelse(is.na(dt$Case), "", dt$Case)
dt$Control <- ifelse(is.na(dt$Control), "", dt$Control)
dt$se <- (dt$uci-dt$or)/1.96
dt$'' <- paste(rep(" ", 20), collapse = " ")

# Create confidence interval column to display 
dt$'OR (95% CI)' <- ifelse(is.na(dt$se), "",
                           sprintf("%.3f (%.3f to %.3f)",
                                   dt$or, dt$lci, dt$uci))

# Define theme
tm <- forest_theme(base_size = 10, #基本字体大小
                   refline_col = "red",
                   arrow_type = "closed",
                   footnote_col = "blue",
                     ci_lwd = 1.5, # Set the line width of the confidence interval
                     ci_Theight = 0.2, # Set an T end at the end of CI 
)

p <- forest(dt[,c(1:3,8:9)],
            est = dt$or,
            lower = dt$lci, 
            upper = dt$uci,
            sizes = dt$se,
            ci_column = 4,
            ref_line = 1.0,
            xlim = c(0.9, 1.25),
            ticks_at = c(0.9,1.0,1.1,1.25),
            theme = tm)


p <- edit_plot(p, row = c(1,3,7,19,22,25), gp = gpar(fontfamily = "sans", fontface = "bold.italic"))

p <- edit_plot(p, row = c(2,8,9,10,12,13,14,16,17,23,24,26,27,4), which = "background",
               gp = gpar(fill = "darkolivegreen2"))
# Print plot
plot(p)

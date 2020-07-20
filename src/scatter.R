library(ggplot2)

df = read.csv("combined_corr_pca_MI_colocal.csv")

# set the max overlap to 10000 to fix the color
df[df$No.overlap > 10000, ]$No.overlap = 10000

ggplot(df, aes(x=MACMIC, y=MI)) + ylim(0,1) +
    geom_point(aes(color=No.overlap)) + scale_color_gradient(low="#1C8AF9", high="red") +
    theme_class() + theme(plot.title = element_blank(),
                          axis.title.x = element_blank(),
                          axis.title.y = element_blank())


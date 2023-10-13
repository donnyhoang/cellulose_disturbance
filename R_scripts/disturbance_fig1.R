library(ggplot2)
library(cowplot)
library(magick)
library(tiff)
library(jpeg)


dump <- readJPEG("refuse_dump.jpg")
schematic <- readTIFF("schematic.tif")

p1 <- ggdraw() +
  draw_image(dump, scale = 0.9)
p2 <- ggdraw() +
  draw_image(schematic, scale = 1)
fig <- plot_grid(p1, p2, labels = "AUTO", label_size = 15, rel_widths = c(5, 10))
fig

ggsave("Figure_1.tiff", device = "tiff", dpi = 700)

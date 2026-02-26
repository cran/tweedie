options(bitmapType = "cairo")
library(hexSticker)
library(ggplot2)
library(tweedie)
library(magick)

## ---- Generate a stylised Tweedie density ----
set.seed(1)
x  <- seq(0.001, 4, length.out = 400)
mu <- 1.2; phi <- 1; p <- 1.05

y  <- dtweedie(x, mu = mu, phi = phi, power = p)
y0 <- dtweedie(0, mu = mu, phi = phi, power = p)
df <- data.frame(x = x, y = y)

density_plot <- ggplot(df, aes(x, y)) +
  geom_line(linewidth = 1.2, colour = "#2C3E50") +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  annotate("point", x = 0, y = y0, color = "#2C3E50", size = 1) +
  theme_void()

## ---- Create hex sticker ----
options(bitmapType = "cairo")
library(hexSticker)
library(ggplot2)
library(tweedie)
library(magick)

## ---- Generate a stylised Tweedie density ----
set.seed(1)
x  <- seq(0.001, 4, length.out = 400)
mu <- 1.2; phi <- 1; p <- 1.05

y  <- dtweedie(x, mu = mu, phi = phi, power = p)
y0 <- dtweedie(0, mu = mu, phi = phi, power = p)
df <- data.frame(x = x, y = y)

density_plot <- ggplot(df, aes(x, y)) +
  geom_line(linewidth = 1.2, colour = "#2C3E50") +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  annotate("point", x = 0, y = y0, color = "#2C3E50", size = 1) +
  theme_void()


## ---- Create hex sticker ----
s <- sticker(
  subplot              = density_plot,  # the ggplot object to embed in the hex
  package              = "tweedie",     # text label printed on the sticker
  p_size               = 25,            # font size of the package label
  p_color              = "#2C3E50",     # colour of the package label
  p_y                  = 1.55,          # vertical position of label (1=centre, >1 = higher)
  s_x                  = 1,             # horizontal centre of the subplot within the hex
  s_y                  = 0.85,          # vertical centre of the subplot (pushed below label)
  s_width              = 1.1,           # width of the subplot as a fraction of the hex width
  s_height             = 0.85,          # height of the subplot as a fraction of the hex height
  h_fill               = "white",       # background fill colour of the hex
  h_color              = "#2C3E50",     # border colour of the hex
  h_size               = 2,             # thickness of the hex border
  dpi                  = 300,           # resolution of the output PNG
  white_around_sticker = TRUE,          # mask outside the hex to white (not transparent)
  filename             = "inst/hex/tweedie_hex.png"  # where to save the PNG
)


## ---- Fix left border clipping with wider canvas ----
ggsave(
  "inst/hex/tweedie_hex.png",
  plot   = s,
  width  = 55,
  height = 55,
  units  = "mm",
  dpi    = 300,
  bg     = "white"
)

## ---- Copy for pkgdown / README ----
file.copy("inst/hex/tweedie_hex.png", "man/figures/logo.png",
          overwrite = TRUE)

## ---- Fix left border clipping with wider canvas ----
ggsave(
  "inst/hex/tweedie_hex.png",
  plot   = s,
  width  = 55,
  height = 55,
  units  = "mm",
  dpi    = 300,
  bg     = "white"
)

## ---- Copy for pkgdown / README ----
file.copy("inst/hex/tweedie_hex.png", "man/figures/logo.png",
          overwrite = TRUE)

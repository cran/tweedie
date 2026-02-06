options(bitmapType = "cairo")

library(hexSticker)
library(ggplot2)
library(tweedie)

## ---- Generate a stylised Tweedie density ----

set.seed(1)

x <- seq(0.001, 4, length.out = 400)

## Tweedie parameters (chosen for visual shape, not pedagogy)
mu   <- 1.2
phi  <- 1
p    <- 1.05

y <- dtweedie(x, 
              mu = mu, 
              phi = phi, 
              power = p)

df <- data.frame(x = x, y = y)

density_plot <- ggplot(df, aes(x, y)) +
  geom_line(linewidth = 1.2, 
            colour = "#2C3E50") +
  coord_cartesian(
    xlim = c(0, max(x)),
    expand = FALSE
  ) +
  theme_void()


## ---- Create hex sticker ----
options(bitmapType = "cairo")

sticker(
  subplot    = density_plot,
  package    = "tweedie",
  p_size     = 20,
  p_color    = "#2C3E50",
  s_x        = 1,
  s_y        = 0.8,
  s_width    = 0.7,
  h_fill     = "white",
  h_color    = "#2C3E50",
  h_size     = 2,
  #white_around = 0.08,
  white_around_sticker = TRUE,
  filename   = "tweedie_hex.pdf"
)

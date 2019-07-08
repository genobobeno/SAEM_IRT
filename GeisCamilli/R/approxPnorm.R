approxPnorm <- function(u, h = 0.2) {
  x <- seq(from = -4, to = 4, by = h)
  approx(x, pnorm(x), yleft = 0, yright = 1, xout = u)$y
}
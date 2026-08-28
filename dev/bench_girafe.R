# Time ggiraph conversion for a lollipop-scale plot (demo: 20 patients x 2 clones).
suppressMessages({library(ggplot2); library(ggiraph)})
set.seed(1)
d <- expand.grid(patient = paste0("P", 1:20), clone = c("c1", "c2"))
d$v <- rnorm(nrow(d)); d$w <- runif(nrow(d))
p <- ggplot(d, aes(x = clone, y = v, size = w)) +
  geom_point() +
  facet_wrap(~ patient, scales = "free_y") +
  theme_bw()
t <- system.time(g <- girafe(ggobj = p))
cat("girafe conversion (20 panels, 40 pts):", t[["elapsed"]], "sec\n")
cat("svg widget size:", format(object.size(g), units = "KB"), "\n")

# Larger: 200 clones x 20 patients = 4000 pts
d2 <- expand.grid(patient = paste0("P", 1:20), clone = paste0("cl", 1:200))
d2$v <- rnorm(nrow(d2)); d2$w <- runif(nrow(d2))
p2 <- ggplot(d2, aes(x = clone, y = v, size = w)) + geom_point() +
  facet_wrap(~ patient, scales = "free_y") + theme_bw()
t2 <- system.time(g2 <- girafe(ggobj = p2))
cat("girafe conversion (20 panels, 4000 pts):", t2[["elapsed"]], "sec\n")
cat("svg widget size:", format(object.size(g2), units = "MB"), "\n")

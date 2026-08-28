md <- readLines("docs/PERCEPTION-shiny.md", warn = FALSE)
cat("lines:", length(md), "\n")
hits <- grep("remit", md)
cat("lines containing 'remit':", length(hits), "\n")
if (length(hits) > 0) cat(md[hits[1]], "\n", sep = "\n")
# extract URLs
pat <- "https://img\\.remit\\.ee/api/file/[A-Za-z0-9_-]+"
m <- gregexpr(pat, md, perl = TRUE)
u <- unique(unlist(regmatches(md, m)))
cat("URL count:", length(u), "\n")
cat(head(u, 2), sep = "\n")

# copy jdx from NIST to excel for teaching purposes

txt <- "copy ##PEAK TABLE=(XY..XY) text here or in console"
txt2 <- gsub("\n", " ", txt)
txt3 <- gsub(" ", "\n", txt2 )

sink(type = "message")
cat(txt3)
sink(type = "message")

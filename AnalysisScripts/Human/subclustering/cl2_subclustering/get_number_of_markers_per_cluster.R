args = commandArgs(trailingOnly=TRUE)
suppressWarnings(suppressMessages(library(dplyr)))

m.out <- readRDS(args[1])
m.out <- m.out$statistics

for (i in 1:length(m.out)) {
print(as.data.frame(m.out[[i]]) %>% dplyr::filter(cohen.mean>=0.5) %>% nrow)
}


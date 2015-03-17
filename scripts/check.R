library(dplyr)
library(readr)

in_fn = "out.table"

t = read_delim(in_fn, delim=" ")

x = t %>% filter(mapq > 20) %>% group_by(chrom, pos, orientation, insertion_end) %>% summarise(count=n())

write.table(x, file="ACH2-clone3-sites.txt", row.names=FALSE, quote=FALSE, col.names=TRUE,
            sep="\t")

dim(t)
length(unique(t$read_name))

# A file containing commonly used functions created or modified by Jared Flater.
# July 10 2019

## Functions
# Get a list of OTUs from a phyloseq object, requires subseting to desired samples first, this will remove OTUs observed less than 5 times in the group of samples
# Note the column Treatment_Response is unique to incubation microcosm data
GetOTUs <- function(physeq, samples) {
    prune_samples(sample_data(physeq)$Treatment_Response %in% c(samples), physeq) %>%
    filter_taxa(function(x) sum(x) >= 5, T) %>%
    tax_table() %>%
    row.names()
}

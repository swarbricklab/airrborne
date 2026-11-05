# Select an observed sequence using transferred support and stable content ties.
# IGH callers supply library-member support, without pooling identical sequences
# across libraries. Non-IGH callers supply the library clonotype observation.
observed_representative_rows <- function(group, support, junction, sequence,
                                          v_gene, j_gene, record_id) {
  ranked <- order(-as_read_support(support), as.character(junction),
                  as.character(sequence), as.character(v_gene),
                  as.character(j_gene), as.character(record_id))
  ranked[!duplicated(as.character(group)[ranked])]
}

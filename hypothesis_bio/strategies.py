from hypothesis import given, example, strategies as st


# TODO: we are only testing with fixed k for now
KSIZE = 13

# strategy for creating kmers. Alphabet is derived from nucleotides.
kmer = st.text("ACGT", min_size=KSIZE, max_size=KSIZE)

sequence = st.text("ACGT", min_size=KSIZE, max_size=1000)

# strategy for creating valid FASTA records
record = st.fixed_dictionaries({
                 'name': st.characters(min_codepoint=32, max_codepoint=126),
                 'sequence': sequence}
                                              )

invalid_record = st.fixed_dictionaries({
                 'name': st.characters(),
                 'sequence': st.characters()}

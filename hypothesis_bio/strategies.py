from hypothesis import strategies as st


# TODO: we are only testing with fixed k for now
KSIZE = 13


# strategy for creating kmers. Alphabet is derived from nucleotides.
kmer = st.text("ACGT", min_size=KSIZE, max_size=KSIZE)


@st.composite
def sequence(draw):
    """
    A :mod:`hypothesis` strategy for building nucleotide sequences

    Parameters
    ----------
    draw
        For internal hypothesis use.

    Returns
    -------
    string
        a string representing a nucleotide sequence
    """

    return draw(st.text("ACGT", min_size=KSIZE, max_size=1000))


# strategy for creating valid FASTA records
record = st.fixed_dictionaries(
    {"name": st.characters(min_codepoint=32, max_codepoint=126), "sequence": sequence}
)

invalid_record = st.fixed_dictionaries(
    {"name": st.characters(), "sequence": st.characters()}
)

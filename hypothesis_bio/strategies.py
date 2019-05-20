from hypothesis import assume
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


CANONICAL_GENE_CODE = dict(
    AAA="K",
    AAC="N",
    AAG="K",
    AAT="N",
    ACA="T",
    ACC="T",
    ACG="T",
    ACT="T",
    AGA="R",
    AGC="S",
    AGG="R",
    AGT="S",
    ATA="I",
    ATC="I",
    ATG="M",
    ATT="I",
    CAA="Q",
    CAC="H",
    CAG="Q",
    CAT="H",
    CCA="P",
    CCC="P",
    CCG="P",
    CCT="P",
    CGA="R",
    CGC="R",
    CGG="R",
    CGT="R",
    CTA="L",
    CTC="L",
    CTG="L",
    CTT="L",
    GAA="E",
    GAC="D",
    GAG="E",
    GAT="D",
    GCA="A",
    GCC="A",
    GCG="A",
    GCT="A",
    GGA="G",
    GGC="G",
    GGG="G",
    GGT="G",
    GTA="V",
    GTC="V",
    GTG="V",
    GTT="V",
    TAA="*",
    TAC="Y",
    TAG="*",
    TAT="Y",
    TCA="S",
    TCC="S",
    TCG="S",
    TCT="S",
    TGA="*",
    TGC="C",
    TGG="W",
    TGT="C",
    TTA="L",
    TTC="F",
    TTG="L",
    TTT="F",
)

ALL_CODONS = tuple(CANONICAL_GENE_CODE.keys())
CANONICAL_STOP_CODONS = ("TAA", "TGA", "TAG")


def stop_codon(stop_codons=CANONICAL_STOP_CODONS):
    """
    A :mod:`hypothesis` strategy for building stop codons

    Parameters
    ----------
    draw
        For internal hypothesis use.

    stop_codons: list[string]
            a list of stop codons (defaults to the canonical TAA, TGA, and TAG)

    Returns
    -------
    string
        a single stop codon
    ----------
    """
    return st.sampled_from(stop_codons)


def codon():
    """
    A :mod:`hypothesis` strategy for getting one codon 

    Parameters
    ----------
    draw
        For internal hypothesis use.

    Returns
    -------
    string
        a single codon
    ----------
    """
    return st.sampled_from(ALL_CODONS)


def non_stop_codon(stop_codons=CANONICAL_STOP_CODONS):
    """
    A :mod:`hypothesis` strategy for building non-stop codons

    Parameters
    ----------
    draw
        For internal hypothesis use.

    stop_codons: list[string]
        The list of accepted stop codons

    Returns
    -------
    string
        a single non-stop codon
    ----------
    """
    return codon().filter(lambda cdn: not cdn in stop_codons)


def coding_sequence(
    min_size=None,
    max_size=None,
    include_stop_codon=True,
    include_start_codon=True,
    allow_internal_stops=False,
):
    """
    A :mod:`hypothesis` strategy for building coding sequences

    Parameters
    ----------
    min_size: int
        The minimum number of codons that will be generated (not including the
        stop codon). This number must be greater than or equal to zero. A
        sequence of length 0 will be either an empty string or a single stop
        codon.

    max_size: int
        The maximum number of non-stop codons in the sequence. max_size must be
        greater than or equal to min_size.

    include_stop_codon: bool
        If True, append a stop codon

    include_start_codon: bool
        If True, append a start codon (ATG)

    allow_internal_stops: bool
        If True, allow stop codons inside the coding sequence

    Returns
    -------
    string 
        a single coding sequence
    ----------
    """
    if include_start_codon:
        if min_size:
            min_size -= 1
        if max_size:
            max_size -= 1
    return st.builds(
        "{}{}{}".format,
        st.just("ATG" if include_start_codon else ""),
        st.lists(non_stop_codon() if include_stop_codon else non_stop_codon(), min_size=min_size, max_size=max_size).map("".join),
        stop_codon() if include_stop_codon else st.just(""),
    )

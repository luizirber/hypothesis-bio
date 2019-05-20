from enum import Enum

from hypothesis import assume
from hypothesis import strategies as st

DEFAULT_KSIZE = 21


class Alphabet(Enum):
    DNA = 1
    DNA_N = 2
    DNA_IUPAC = 3
    RNA = 4
    RNA_N = 5
    RNA_IUPAC = 6
    PROTEIN = 7


SYMBOLS = {
    Alphabet.DNA: "ACGT",
    Alphabet.DNA_N: "ACGTN",
    Alphabet.DNA_IUPAC: "AGCTYRWSKMDVHBNZ",
    Alphabet.RNA: "ACGU",
    Alphabet.RNA_N: "ACGUN",
    Alphabet.RNA_IUPAC: "AGCUYRWSKMDVHBNZ",
    Alphabet.PROTEIN: "ARNDCEQGHILKMFPSTWYV",
}


def symbols(alphabet: Alphabet) -> str:
    return SYMBOLS[alphabet]


REVERSE_COMPLEMENTS = {
    Alphabet.DNA: dict(zip(symbols(Alphabet.DNA), "TGCA")),
    Alphabet.DNA_N: dict(zip(symbols(Alphabet.DNA_N), "TGCAN")),
    Alphabet.DNA_IUPAC: dict(zip(symbols(Alphabet.DNA_IUPAC), "TCGARYWSMKHBDVNZ")),
    Alphabet.RNA: dict(zip(symbols(Alphabet.RNA), "ACGU")),
    Alphabet.RNA_N: dict(zip(symbols(Alphabet.RNA_N), "ACGUN")),
    Alphabet.RNA_IUPAC: dict(zip(symbols(Alphabet.RNA_IUPAC), "UCGARYWSMKHBDVNZ")),
}


def reverse_complement(seq: str, alphabet: Alphabet = Alphabet.DNA) -> str:
    if alphabet == Alphabet.PROTEIN:
        raise ValueError("Proteins don't have reverse complements")

    RC = REVERSE_COMPLEMENTS[alphabet]
    return "".join(RC[c] for c in reversed(seq))


@st.composite
def kmer(draw, *, alphabet: Alphabet = Alphabet.DNA, k: int = DEFAULT_KSIZE):
    """
    A :mod:`hypothesis` strategy for creating k-mers, short sliding window
    substrings from a sequence.

    Parameters
    ----------
    draw
        For internal hypothesis use.
    alphabet: Alphabet
        Symbols used to generate each position. Defaults to DNA symbols (ACGT)
    ksize: int
        Substring size. Must be positive.

    Returns
    -------
    hypothesis.searchstrategy
        a strategy for generating a string with length k
    """

    return draw(st.text(symbols(alphabet), min_size=k, max_size=k))


@st.composite
def kmers(
    draw,
    seq: str,
    *,
    alphabet: Alphabet = Alphabet.DNA,
    k: int = DEFAULT_KSIZE,
    rc: bool = False
):
    """
    A :mod:`hypothesis` strategy for creating k-mers (short sliding window
    substrings) from a given sequence.

    Parameters
    ----------
    draw
        For internal hypothesis use.
    sequence: str
        Base sequence to sample k-mers from. Must be at least ksize long.
    alphabet: Alphabet
        Symbols used to generate each position. Defaults to DNA symbols (ACGT)
    ksize: int
        Substring size. Must be positive. Defaults to 21.
    rc: bool
        Allow reverse complements to be returned. Defaults to False.

    Returns
    -------
    hypothesis.searchstrategy
        a strategy for generating a string with length k
    """
    assert len(seq) >= k

    i = draw(st.integers(min_value=0, max_value=len(seq) - k))
    kmer = seq[i : i + k]

    return_rc = False
    if rc:
        return_rc = draw(st.booleans())

    if return_rc:
        return reverse_complement(kmer, alphabet)

    return kmer


@st.composite
def sequence(draw, *, alphabet: Alphabet = Alphabet.DNA, max_size: int = 1000):
    """
    A :mod:`hypothesis` strategy for building nucleotide sequences

    Parameters
    ----------
    draw
        For internal hypothesis use.
    alphabet: Alphabet
        Symbols used to generate each position. Defaults to DNA symbols (ACGT)
    max_size: int
        Length of the generated sequence. Must be positive.

    Returns
    -------
    hypothesis.searchstrategy
        a strategy for generating a string representing a nucleotide sequence
    """

    return draw(st.text(symbols(alphabet), max_size=max_size))


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

    stop_codons: list[string]
            a list of stop codons (defaults to the canonical TAA, TGA, and TAG)

    Returns
    -------
    hypothesis.searchstrategy
        a strategy for generating a single stop codon
    ----------
    """
    return st.sampled_from(stop_codons)


def codon():
    """
    A :mod:`hypothesis` strategy for getting one codon

    Returns
    -------
    hypothesis.searchstrategy
        a strategy for generating a single codon
    ----------
    """
    return st.sampled_from(ALL_CODONS)


def non_stop_codon(stop_codons=CANONICAL_STOP_CODONS):
    """
    A :mod:`hypothesis` strategy for building non-stop codons

    Parameters
    ----------
    stop_codons: list[string]
        The list of accepted stop codons

    Returns
    -------
    hypothesis.searchstrategy
        a strategy for generating a single non-stop codon
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
    hypothesis.searchstrategy
        a strategy for generating a single coding sequence
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
        st.lists(
            non_stop_codon() if include_stop_codon else non_stop_codon(),
            min_size=min_size,
            max_size=max_size,
        ).map("".join),
        stop_codon() if include_stop_codon else st.just(""),
    )

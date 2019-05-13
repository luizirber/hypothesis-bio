from enum import Enum

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
    string
        a string with length k
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
    string
        a string with length k
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
    string
        a string representing a nucleotide sequence
    """

    return draw(st.text(symbols(alphabet), max_size=max_size))


# strategy for creating valid FASTA records
record = st.fixed_dictionaries(
    {"name": st.characters(min_codepoint=32, max_codepoint=126), "sequence": sequence}
)

invalid_record = st.fixed_dictionaries(
    {"name": st.characters(), "sequence": st.characters()}
)

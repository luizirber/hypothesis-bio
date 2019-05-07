from enum import Enum, auto

from hypothesis import strategies as st


DEFAULT_KSIZE = 21


class Alphabet(Enum):
    DNA = auto()
    DNA_N = auto()
    DNA_IUPAC = auto()
    RNA = auto()
    RNA_N = auto()
    RNA_IUPAC = auto()
    PROTEIN = auto()


SYMBOLS = {
    Alphabet.DNA: "ACGT",
    Alphabet.DNA_N: "ACGTN",
    Alphabet.DNA_IUPAC: "ACGTN",
    Alphabet.RNA: "ACGU",
    Alphabet.RNA_N: "ACGUN",
    Alphabet.RNA_IUPAC: "ACGURYSWKMBDHVNZ",
    Alphabet.PROTEIN: "ARNDCEQGHILKMFPSTWYV"
}


def symbols(alphabet: Alphabet) -> str:
    return SYMBOLS[alphabet]


def complement(c: str, alphabet: Alphabet = Alphabet.DNA) -> str:
    raise NotImplementedError()



@st.composite
def kmer(draw, *, alphabet="ACGT", k=DEFAULT_KSIZE):
    """
    A :mod:`hypothesis` strategy for creating k-mers, short sliding window
    substrings from a sequence.

    Parameters
    ----------
    draw
        For internal hypothesis use.
    ksize: int
        Substring size. Must be positive.

    Returns
    -------
    string
        a string with length k
    """
    # strategy for creating kmers. Alphabet is derived from nucleotides.
    return draw(st.text(alphabet, min_size=k, max_size=k))


@st.composite
def sequence(draw, *, alphabet="ACGT", max_size=1000):
    """
    A :mod:`hypothesis` strategy for building nucleotide sequences

    Parameters
    ----------
    draw
        For internal hypothesis use.
    max_size: int
        Length of the generated sequence. Must be positive.

    Returns
    -------
    string
        a string representing a nucleotide sequence
    """

    return draw(st.text(alphabet, max_size=max_size))


# strategy for creating valid FASTA records
record = st.fixed_dictionaries(
    {"name": st.characters(min_codepoint=32, max_codepoint=126), "sequence": sequence}
)

invalid_record = st.fixed_dictionaries(
    {"name": st.characters(), "sequence": st.characters()}
)

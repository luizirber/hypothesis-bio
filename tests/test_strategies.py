from hypothesis import given, assume
from hypothesis import strategies as st

from hypothesis_bio.strategies import (
    sequence,
    Alphabet,
    symbols,
    kmers,
    reverse_complement,
)


@given(st.data(), st.sampled_from(Alphabet))
def test_sequence_random(data, alphabet):
    seq_string = data.draw(sequence(alphabet=alphabet))

    assert isinstance(seq_string, str)
    assert all(nt in symbols(alphabet) for nt in seq_string)


@given(st.data(), st.sampled_from(Alphabet), st.booleans())
def test_kmers(data, alphabet, rc):
    if rc:
        assume(alphabet != Alphabet.PROTEIN)

    seq = data.draw(sequence(alphabet=alphabet))
    ksize = data.draw(st.integers(0, len(seq)))
    kmer_s = data.draw(kmers(seq, alphabet=alphabet, k=ksize, rc=rc))

    if rc:
        assert kmer_s in seq or reverse_complement(kmer_s, alphabet) in seq
    else:
        assert kmer_s in seq

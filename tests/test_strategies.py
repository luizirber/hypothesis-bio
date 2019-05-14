import pytest
from hypothesis import assume, given
from hypothesis import strategies as st

from hypothesis_bio.strategies import (
    Alphabet,
    kmer,
    kmers,
    reverse_complement,
    sequence,
    symbols,
)


@given(st.data(), st.sampled_from(Alphabet))
def test_sequence(data, alphabet):
    seq_string = data.draw(sequence(alphabet=alphabet))

    assert isinstance(seq_string, str)
    assert all(nt in symbols(alphabet) for nt in seq_string)


@given(st.data(), st.sampled_from(Alphabet), st.integers(min_value=0, max_value=150))
def test_kmer(data, alphabet, ksize):
    kmer_s = data.draw(kmer(alphabet=alphabet, k=ksize))

    assert isinstance(kmer_s, str)
    assert all(nt in symbols(alphabet) for nt in kmer_s)
    assert len(kmer_s) == ksize


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


def test_invalid_rc():
    with pytest.raises(ValueError):
        reverse_complement("A", alphabet=Alphabet.PROTEIN)

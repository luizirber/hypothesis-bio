from hypothesis import given
from hypothesis import strategies as st

from hypothesis_bio.strategies import sequence


@given(data=st.data())
def test_sequence_random(data):
    seq_string = data.draw(sequence())

    assert isinstance(seq_string, str)
    assert all(nt in "ACGT" for nt in seq_string)

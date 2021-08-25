from duplex_tools.split_on_adapter import find_mid_adaptor
from hypothesis import strategies as st, given, settings


@given(seq=st.text(alphabet='ACGT', min_size=0, max_size=3),
       target=st.text(alphabet='ACGT', min_size=0, max_size=3))
@settings(max_examples=30)
def test_find_mid_adaptor_noerror_empty_or_short_sequences(seq, target):
    print(f'seq   :{seq}')
    print(f'target:{target}')
    res = find_mid_adaptor(seq, targets=[target])
    print(f'result:{res}')


def test_find_mid_adaptor_edlib():
    middle_seq = "AGTCGTGTCA"
    padding = "GTGTGGTGTG" * 20
    seq = f"{padding}{middle_seq}{padding}"
    res = find_mid_adaptor(seq, [middle_seq], print_alignment=True, print_threshold=12)
    assert res['editDistance'] == 0

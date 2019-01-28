from invariants.tube import Tube

def test_tube_exists():
    t = Tube([1,2,3], [1,0,0])
    assert t is not None
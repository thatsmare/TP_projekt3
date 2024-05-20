from __future__ import annotations

import scikit_build_example as m


def test_version():
    assert m.__version__ == "0.0.1"


def test_add():
    assert m.add(1, 2) == 3


def test_sub():
    assert m.subtract(1, 2) == -1

def test_signal_generate_sinusoidal():
    assert m.signal_generate_sinusoidal(5,1,0,15,1)                              

def test_signal_generate_square_wave():
    assert m.signal_generate_square_wave(5,2,19,50,1) 
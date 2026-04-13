import more_itertools
import pytest

from cursemol import canvas
from cursemol.state import ScreenDimensions
from cursemol.state import State


@pytest.fixture
def screen_dims():
    return ScreenDimensions(40, 20)


@pytest.fixture
def state(screen_dims):
    state, _ = State.fromSmiles('CCO', screen_dims)
    return state


def compare_chars(got, expected):
    """
    Assert if we got the drawing that we expected. `got` is a list of list of
    char as returned by fill_screen_buffer(), but `expected` is just a string
    for conciseness. Blank lines and trailing whitespace are ignored.
    """
    got = [''.join(line).rstrip() for line in got]
    expected = [line.rstrip() for line in expected.split('\n')]
    got = [line for line in got if line]
    expected = [line for line in expected if line]
    assert got == expected


def test_CCO(state, screen_dims):
    chars, colors = canvas.fill_screen_buffer(state, screen_dims)
    expected_chars = """
                   ·C·
               ····   ····
           ····           ···
         C·                  ·OH"""
    compare_chars(chars, expected_chars)

    # To keep it simple we'll just count how many red bold cells we have.
    num_bold_red = sum(c == 257 for c in more_itertools.flatten(colors))
    assert num_bold_red == 2


def test_CCO_scaled(state, screen_dims):
    state.scale = (4.0, 1.6)
    chars, colors = canvas.fill_screen_buffer(state, screen_dims)
    expected_chars = """
          ··C··
       C··     ··OH"""
    compare_chars(chars, expected_chars)


def test_shift_view(state):
    (xmin, ymin, zmin), (xmax, ymax, zmax) = state.box

    # These are the defaults; check just in case to make failures more obvious.
    assert state.scale == (8.0, 3.2)

    canvas.shift_view(state, +2, -1)

    assert state.box[0][0] == pytest.approx(xmin + 0.25)
    assert state.box[1][0] == pytest.approx(xmax + 0.25)
    assert state.box[0][1] == pytest.approx(ymin - 0.3125)
    assert state.box[1][1] == pytest.approx(ymax - 0.3125)

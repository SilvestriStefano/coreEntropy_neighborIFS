import pytest
from src import functions
from pytest import approx, raises
# import pytest

@pytest.mark.parametrize("test_num_exact,test_den_exact,expected",[
    (0,1,0.0),
    (1,2,2.0),
    (1,1,0.0)
    ])


def test_core_entropy_exact(test_num_exact,test_den_exact,expected):
    """
    check that the core_entropy returns the exact values for allowable input
    """

    assert functions.core_entropy(num=test_num_exact,den=test_den_exact) == expected

@pytest.mark.parametrize("test_num_approx,test_den_approx,expected",[
    (1,4,1.69562),
    (2,7,1.0),
    (1,5,1.7455312)
    ])
def test_core_entropy_approx(test_num_approx,test_den_approx,expected):
    """
    check that the core_entropy returns (approximately) correct values for allowable input
    """

    assert functions.core_entropy(num=test_num_approx,den=test_den_approx) == approx(expected)

@pytest.mark.parametrize("test_num_val,test_den_val",[
    ("numerator","denominator"),
    ("numerator","1"),
    ("numerator",1.0),
    ("numerator",1),
    ("1","denominator"),
    (1.0,"denominator"),
    (1,"denominator"),
    ("1","1.0"),
    (3+9j,1+0j),
    (1+0j,1+0j)
    ])
def test_core_entropy_ValueErrors(test_num_val,test_den_val):
    """
    check that the core_entropy raises ValueError for not allowed input
    """
    
    with raises(ValueError): 
        functions.core_entropy(num=test_num_val,den=test_den_val)
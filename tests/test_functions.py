from src.functions import core_entropy
from pytest import approx, raises, mark

@mark.parametrize("test_num_exact,test_den_exact,expected",[
    (0,1,1.0),
    (1,2,2.0),
    (1,1,1.0)
    ],ids=["0/1","1/2","1/1"])
def test_core_entropy_exact(test_num_exact,test_den_exact,expected):
    """
    check that the core_entropy returns the exact values for allowable input
    """

    assert core_entropy(num=test_num_exact,den=test_den_exact) == expected

@mark.parametrize("test_num_approx,test_den_approx,expected",[
    (1,4,1.69562),
    (2,7,1.0),
    (1,5,1.3953369),
    (3,14,1.6180339)
    ],ids=["1/4","2/7","1/5","3/14"])
def test_core_entropy_approx(test_num_approx,test_den_approx,expected):
    """
    check that the core_entropy returns (approximately) correct values for allowable input
    """

    assert core_entropy(num=test_num_approx,den=test_den_approx) == approx(expected)

@mark.parametrize("test_num_val,test_den_val",[
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
    ],ids=["str str","str char","str float","str int","char str","float str","int str","char char","3+9i 1+0i","1+0i 1+0i"])
def test_core_entropy_ValueErrors(test_num_val,test_den_val):
    """
    check that the core_entropy raises ValueError for not allowed input
    """
    
    with raises(ValueError): 
        core_entropy(num=test_num_val,den=test_den_val)
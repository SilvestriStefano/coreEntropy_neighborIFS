from src.functions import core_entropy, neighbor_graph
from src.angles import Angle
from pytest import approx, raises, mark

@mark.parametrize("test_num_exact,test_den_exact,expected",[
    (0,1,1.0),
    (1,2,2.0),
    (1,1,1.0)
    ],ids=["0/1","1/2","1/1"])
def test_core_entropy_numden_exact(test_num_exact,test_den_exact,expected):
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

@mark.parametrize("test_angle_exact,expected",[
    (Angle(1,2),2.0),
    (Angle(th="2/7"),1.0),
    (Angle(3,14),1.6180339)
    ],ids=["1/2","2/7","3/14"])
def test_core_entropy_angle(test_angle_exact,expected):
    """
    check that the core_entropy returns approx values for Angle input
    """

    assert core_entropy(angle=test_angle_exact) == approx(expected)

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
def test_core_entropy_numden_ValueErrors(test_num_val,test_den_val):
    """
    check that the core_entropy raises ValueError for not allowed input
    """
    
    with raises(ValueError): 
        core_entropy(num=test_num_val,den=test_den_val)

@mark.parametrize("test_angle",[
    "numerator",
    "1",
    1.0,
    1,
    3+9j
    ],ids=["str","char","float","int","complex"])
def test_core_entropy_angle_ValueErrors(test_angle):
    """
    check that the core_entropy raises ValueError for not allowed input
    """
    
    with raises(ValueError): 
        core_entropy(angle=test_angle)


lambdas = [
    (0.25+0.15*1j,"disconnected"),
    (0.5**(0.5)*(1+1j),"on the unit disk"),
    (0.5+0.*1j,"1/2"),
    (0.366875964264129394095471688839696517595563618765593247116005 + 0.520259438865200828930723528940227115174924275956965606630538*1j,"3/14"),
    (0.595743941976559373530677134113758646737345510937144045504891 + 0.254425889416368995243211219643495532619603711040766445706284*1j,"3/8"),
    (0.102784715200295155851014308891911713188554459798849720235276 + 0.665456951152813476706190611629077711878015849981853233566949*1j,"Kolakoski"),
    (0.+0.707106781186547524400844362104849039284835937688474036588340*1j,"A4 paper")
]
expected_nbhg = [
    {},
    {},
    {
        'h+': {'h+': 'pm'},
        'id': {'h+': 'mp'}
    },
    {
        'h+': {'h+-': '*'},
        'h+-': {'h+-+': 'mp'},
        'h+-+': {'h+-++': 'mp'},
        'h+-++': {'h+-+++': '*'},
        'h+-+++': {'h+-+++-': 'pm'},
        'h+-+++-': {'h+': 'pm'},
        'id': {'h+': 'mp'}
    },
    {
        'h+': {'h+-': '*'},
        'h+-': {'h+--': 'pm', 'h+-0': '*'},
        'h+--': {'h+--': 'mp'},
        'id': {'h+': 'mp'}
    },
    {
        'h+': {'h+-': 'pm', 'h+0': '*'},
        'h+0': {'h+0+': 'mp'},
        'h+0+': {'h+0++': '*'},
        'h+0++': {'h+0++-': 'pm'},
        'h+0++-': {'h+0': '*'},
        'id': {'h+': 'mp'}
    },
    {
        'h+': {'h++': 'mp', 'h+-': 'pm', 'h+0': '*'},
        'h++': {'h+-': 'mp'},
        'h+-': {'h+-+': 'mp'},
        'h+-+': {'h+-+-': 'pm'},
        'h+-+-': {'h++': 'pm'},
        'h+0': {'h+0+': '*'},
        'h+0+': {'h+-+': 'pm', 'h+-+-': 'mp', 'h+0+0': '*'},
        'h+0+0': {'h+': '*'},
        'id': {'h+': 'mp'}
    }

]
@mark.parametrize("test_param,test_depth,expected",
                  [(lam[0],8,ex_ang) for lam, ex_ang in zip(lambdas,expected_nbhg)],
                  ids=[lam[1] for lam in lambdas])
def test_neighbor_graph(test_param,test_depth,expected):
    """
    check it creates the correct neighbor graph
    """
    test_nbh = neighbor_graph(test_param,test_depth)
    assert test_nbh == expected

@mark.parametrize("test_param,test_depth",
                  [("string","string"),
                   ("32","1"),
                   (32,1.0),
                   (3+9j,0.5),
                   (1+0j,1+0j)],
                  ids=["str str","char char","int float","complex float","complex complex"])
def test_neighbor_graph_ValueError(test_param,test_depth):
    """
    check that the neighbor_graph raises ValueError for not allowed input
    """
    
    with raises(ValueError): 
        neighbor_graph(test_param,test_depth)

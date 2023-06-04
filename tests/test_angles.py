from src.angles import Angle
import pytest 

## see documentation at https://docs.pytest.org/en/stable/example/parametrize.html#parametrizing-test-methods-through-per-class-configuration
def pytest_generate_tests(metafunc):
    # called once per each test function
    funcarglist = metafunc.cls.params[metafunc.function.__name__]
    funcidslist = metafunc.cls.ids[metafunc.function.__name__]
    argnames = sorted(funcarglist[0])
    metafunc.parametrize(
        argnames, [[funcargs[name] for name in argnames] for funcargs in funcarglist],
        ids=funcidslist
    )

class TestAngle:
    params = {
        "test_period": [
            {'test_angle': Angle(0,1), 'expected': (1,0)},
            {'test_angle': Angle(1,2), 'expected': (1,1)},
            {'test_angle': Angle(1,6), 'expected': (2,1)}
        ],
        "test_ks_from_angle": [
            {'test_angle': Angle(0,1), 'expected': '*'},
            {'test_angle': Angle(1,2), 'expected': '10'},
            {'test_angle': Angle(1,6), 'expected':  '110'}
        ],
        "test_attr_itin_from_ks": [
            {'test_angle': Angle(0,1), 'expected': '+'},
            {'test_angle': Angle(1,2), 'expected': '+-'},
            {'test_angle': Angle(1,6), 'expected': '+-++--'}
        ],
        "test_period_length_itin": [
            {'test_angle': Angle(0,1), 'expected': 1},
            {'test_angle': Angle(1,2), 'expected': 1},
            {'test_angle': Angle(1,6), 'expected': 4}
        ],
        "test_itin_to_rat": [
            {'test_angle': Angle(0,1), 'expected': '1'},
            {'test_angle': Angle(1,2), 'expected': '(+x^0)*(1-x) +(-x^1)'},
            {'test_angle': Angle(1,6), 'expected': '(+x^0-x^1)*(1-x^4) +(+x^2+x^3-x^4-x^5)'}
        ]
    }
    ids = {
        "test_period":["0/1","1/2","1/6"],
        "test_ks_from_angle":["0/1","1/2","1/6"],
        "test_attr_itin_from_ks":["0/1","1/2","1/6"],
        "test_period_length_itin":["0/1","1/2","1/6"],
        "test_itin_to_rat":["0/1","1/2","1/6"]
    }

    def test_period(self, test_angle,expected):
        assert test_angle.period() == expected

    def test_ks_from_angle(self, test_angle,expected):
        assert test_angle.ks_from_angle() == expected
    
    def test_attr_itin_from_ks(self,test_angle,expected):
        assert test_angle.attr_itin_from_ks() == expected
    
    def test_period_length_itin(self,test_angle,expected):
        assert test_angle.period_length_itin() == expected
    
    def test_itin_to_rat(self,test_angle,expected):
        assert test_angle.itin_to_rat(pow_symb='^') == expected
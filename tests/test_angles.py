from src.angles import Angle
import pytest 


def pytest_generate_tests(metafunc):
    # called once per each test function
    funcarglist = metafunc.cls.params[metafunc.function.__name__]
    argnames = sorted(funcarglist[0])
    metafunc.parametrize(
        argnames, [[funcargs[name] for name in argnames] for funcargs in funcarglist]
    )

class TestAngle:
    params = {
        "test_period": [
            {'test_angle': Angle(0,1), 'expected': (1,0)},
            {'test_angle': Angle(1,2), 'expected': (1,1)} 
        ],
        "test_ks": [
            {'test_angle': Angle(0,1), 'expected': '*'},
            {'test_angle': Angle(1,2), 'expected': '10'}
        ]
    }

    def test_period(self, test_angle,expected):
        assert test_angle.period() == expected

    def test_ks(self, test_angle,expected):
        assert test_angle.ks_from_angle() == expected
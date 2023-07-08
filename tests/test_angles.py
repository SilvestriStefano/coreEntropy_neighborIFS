from src.angles import Angle
from numpy import sqrt
from pytest import approx

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

preperiodic_angles = {
    "1/2": {
        "angle": Angle(1,2),
        "binary": "1p0",
        "period": (1,1),
        "ks_from_angle": "10",
        "attr_itin_from_ks": "+-",
        "period_length_itin": 1,
        "itin_to_rat": "(+x^0)*(1-x) +(-x^1)",
        "assoc_lambda": 0.5+0.j
    },
    "1/6": {
        "angle": Angle(1,6), # Tame Twin Dragon
        "binary": "0p01",
        "period": (2,1),
        "ks_from_angle": "110",
        "attr_itin_from_ks": "+-++--",
        "period_length_itin": 4,
        "itin_to_rat": "(+x^0-x^1)*(1-x^4) +(+x^2+x^3-x^4-x^5)",
        "assoc_lambda": 0.25+0.661437827766148j
    },
    "3/14": {
        "angle": Angle(3,14),
        "binary": "0p011",
        "period": (3,1),
        "ks_from_angle": "1100",
        "attr_itin_from_ks": "+-+++---",
        "period_length_itin": 6,
        "itin_to_rat": "(+x^0-x^1)*(1-x^6) +(+x^2+x^3+x^4-x^5-x^6-x^7)",
        "assoc_lambda": 0.366875964264129 + 0.520259438865201j
    },
    "3/8": {
        "angle": Angle(3,8), #same as 21/56
        "binary": "011p0",
        "period": (1,3),
        "ks_from_angle": "1010",
        "attr_itin_from_ks": "+--+",
        "period_length_itin": 1,
        "itin_to_rat": "(+x^0-x^1-x^2)*(1-x) +(+x^3)",
        "assoc_lambda": 0.595743941976559 + 0.254425889416369j
    },
    "7/16": {
        "angle": Angle(7,16), #same as 49/112
        "binary": "0111p0",
        "period": (1,4),
        "ks_from_angle": "10010",
        "attr_itin_from_ks": "+---+",
        "period_length_itin": 1,
        "itin_to_rat": "(+x^0-x^1-x^2-x^3)*(1-x) +(+x^4)",
        "assoc_lambda": 0.636009824757034 + 0.106924311121288j
    },
    "55/256": {
        "angle": Angle(55,256), #same as 385/1792 [CKW buried parameter]
        "binary": "00110111p0",
        "period": (1,8),
        "ks_from_angle": "110010010",
        "attr_itin_from_ks": "+-+++---+",
        "period_length_itin": 1,
        "itin_to_rat": "(+x^0-x^1+x^2+x^3+x^4-x^5-x^6-x^7)*(1-x) +(+x^8)",
        "assoc_lambda": 0.371858680074136 + 0.519411153747943j
    }
} 

class TestPreperiodicAngle:
    params = {
        "test_period": [{"test_angle": v["angle"], "expected": v["period"]} for _,v in preperiodic_angles.items()],
        "test_binary": [{"test_angle": v["angle"], "expected": v["binary"]} for _,v in preperiodic_angles.items()],
        "test_ks_from_angle": [{"test_angle": v["angle"], "expected": v["ks_from_angle"]} for _,v in preperiodic_angles.items()],
        "test_attr_itin_from_ks": [{"test_angle": v["angle"], "expected": v["attr_itin_from_ks"]} for _,v in preperiodic_angles.items()],
        "test_period_length_itin": [{"test_angle":v["angle"], "expected": v["period_length_itin"]} for _,v in preperiodic_angles.items()],
        "test_itin_to_rat": [{"test_angle":v["angle"], "expected": v["itin_to_rat"]} for _,v in preperiodic_angles.items()],
        "test_assoc_lambda": [{"test_angle":v["angle"], "expected": v["assoc_lambda"]} for _,v in preperiodic_angles.items()]
    }
    ids = {func : [key for key in preperiodic_angles] for func in params}

    def test_period(self, test_angle,expected):
        assert test_angle.period() == expected

    def test_binary(self, test_angle,expected):
        assert test_angle.to_binary() == expected

    def test_ks_from_angle(self, test_angle,expected):
        assert test_angle.ks_from_angle() == expected
    
    def test_attr_itin_from_ks(self,test_angle,expected):
        assert test_angle.attr_itin_from_ks() == expected
    
    def test_period_length_itin(self,test_angle,expected):
        assert test_angle.period_length_itin() == expected
    
    def test_itin_to_rat(self,test_angle,expected):
        assert test_angle.itin_to_rat(pow_symb='^') == expected

    def test_assoc_lambda(self,test_angle,expected):
        lam = complex(test_angle.assoc_lambda()) # need to cast from <class 'sympy.core.add.Add'>
        lam = lam if lam.imag>0 else lam.conjugate() # always want positive imaginary part
        assert lam == approx(expected, rel=1e-15)

periodic_angles = {
    "0/1": {
        "angle": Angle(0,1),
        "binary": "p0",
        "period": (1,0),
        "ks_from_angle": "*",
        "attr_itin_from_ks": None,
        "period_length_itin": None,
        "itin_to_rat": None,
        "assoc_lambda": None
    },
    "1/7": {
        "angle": Angle(1,7),
        "binary": "p001",
        "period": (3,0),
        "ks_from_angle": "11*",
        "attr_itin_from_ks": None,
        "period_length_itin": None,
        "itin_to_rat": None,
        "assoc_lambda": None
    },
    "2/5": {
        "angle": Angle(2,5),
        "binary": "p0110",
        "period": (4,0),
        "ks_from_angle": "101*",
        "attr_itin_from_ks": None,
        "period_length_itin": None,
        "itin_to_rat": None,
        "assoc_lambda": None
    },
}

class TestPeriodicAngle:
    params = {
        "test_period": [{"test_angle": v["angle"], "expected": v["period"]} for _,v in periodic_angles.items()],
        "test_binary": [{"test_angle": v["angle"], "expected": v["binary"]} for _,v in periodic_angles.items()],
        "test_ks_from_angle": [{"test_angle": v["angle"], "expected": v["ks_from_angle"]} for _,v in periodic_angles.items()],
        "test_attr_itin_from_ks": [{"test_angle": v["angle"], "expected": v["attr_itin_from_ks"]} for _,v in periodic_angles.items()],
        "test_period_length_itin": [{"test_angle":v["angle"], "expected": v["period_length_itin"]} for _,v in periodic_angles.items()],
        "test_itin_to_rat": [{"test_angle":v["angle"], "expected": v["itin_to_rat"]} for _,v in periodic_angles.items()],
        "test_assoc_lambda": [{"test_angle":v["angle"], "expected": v["assoc_lambda"]} for _,v in periodic_angles.items()]
    }
    ids = {func : [key for key in periodic_angles] for func in params}
    
    def test_period(self, test_angle,expected):
        assert test_angle.period() == expected

    def test_binary(self, test_angle,expected):
        assert test_angle.to_binary() == expected

    def test_ks_from_angle(self, test_angle,expected):
        assert test_angle.ks_from_angle() == expected
    
    def test_attr_itin_from_ks(self,test_angle,expected):
        assert test_angle.attr_itin_from_ks() == expected
    
    def test_period_length_itin(self,test_angle,expected):
        assert test_angle.period_length_itin() == expected
    
    def test_itin_to_rat(self,test_angle,expected):
        assert test_angle.itin_to_rat() == expected

    def test_assoc_lambda(self,test_angle,expected):
        assert test_angle.assoc_lambda() == expected

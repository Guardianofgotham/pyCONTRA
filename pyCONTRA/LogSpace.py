from pyCONTRA.Utilities import *
import math

NEG_INF = -2e20


# def Fast_Exp(x):
#     if (x <= NEG_INF/2):
#         return 0
#     return math.exp(x)


def Fast_Exp(x):
    # Bounds for tolerance of 4.96e-05: (-9.91152, 0)
    # Approximating interval: (-9.91152, -5.86228) --> ((T(0.0000803850)*x+T(0.0021627428))*x+T(0.0194708555))*x+T(0.0588080014);
    # Approximating interval: (-5.86228, -3.83966) --> ((T(0.0013889414)*x+T(0.0244676474))*x+T(0.1471290604))*x+T(0.3042757740);
    # Approximating interval: (-3.83966, -2.4915) --> ((T(0.0072335607)*x+T(0.0906002677))*x+T(0.3983111356))*x+T(0.6245959221);
    # Approximating interval: (-2.4915, -1.48054) --> ((T(0.0232410351)*x+T(0.2085645908))*x+T(0.6906367911))*x+T(0.8682322329);
    # Approximating interval: (-1.48054, -0.672505) --> ((T(0.0573782771)*x+T(0.3580258429))*x+T(0.9121133217))*x+T(0.9793091728);
    # Approximating interval: (-0.672505, -3.9145e-11) --> ((T(0.1199175927)*x+T(0.4815668234))*x+T(0.9975991939))*x+T(0.9999505077);
    # 6 polynomials needed.

    if (x < float(-2.4915033807)):
        if (x < float(-5.8622823336)):
            if (x < float(-9.91152)):
                return float(0)
            return ((float(0.0000803850)*x+float(0.0021627428))*x+float(0.0194708555))*x+float(0.0588080014)
        if (x < float(-3.8396630909)):
            return ((float(0.0013889414)*x+float(0.0244676474))*x+float(0.1471290604))*x+float(0.3042757740)
        return ((float(0.0072335607)*x+float(0.0906002677))*x+float(0.3983111356))*x+float(0.6245959221)
    if (x < float(-0.6725053211)):
        if (x < float(-1.4805375919)):
            return ((float(0.0232410351)*x+float(0.2085645908))*x+float(0.6906367911))*x+float(0.8682322329)
        return ((float(0.0573782771)*x+float(0.3580258429))*x+float(0.9121133217))*x+float(0.9793091728)
    if (x < float(0)):
        return ((float(0.1199175927)*x+float(0.4815668234))*x+float(0.9975991939))*x+float(0.9999505077)
    return (float(1e20)) if x > float(46.052) else (math.exp(x))


# def Fast_LogExpPlusOne(x):
#     assert 0 <= x and x <= 30, "Argument out-of-range."
#     return math.log(math.exp(x) + float(1))


def Fast_LogExpPlusOne(x):

    # Bounds for tolerance of 7.05e-06: (0, 11.8625)
    # Approximating interval: (0, 0.661537) --> ((T(-0.0065591595)*x+T(0.1276442762))*x+T(0.4996554598))*x+T(0.6931542306);
    # Approximating interval: (0.661537, 1.63202) --> ((T(-0.0155157557)*x+T(0.1446775699))*x+T(0.4882939746))*x+T(0.6958092989);
    # Approximating interval: (1.63202, 2.49126) --> ((T(-0.0128909247)*x+T(0.1301028251))*x+T(0.5150398748))*x+T(0.6795585882);
    # Approximating interval: (2.49126, 3.37925) --> ((T(-0.0072142647)*x+T(0.0877540853))*x+T(0.6208708362))*x+T(0.5909675829);
    # Approximating interval: (3.37925, 4.42617) --> ((T(-0.0031455354)*x+T(0.0467229449))*x+T(0.7592532310))*x+T(0.4348794399);
    # Approximating interval: (4.42617, 5.78907) --> ((T(-0.0010110698)*x+T(0.0185943421))*x+T(0.8831730747))*x+T(0.2523695427);
    # Approximating interval: (5.78907, 7.81627) --> ((T(-0.0001962780)*x+T(0.0046084408))*x+T(0.9634431978))*x+T(0.0983148903);
    # Approximating interval: (7.81627, 11.8625) --> ((T(-0.0000113994)*x+T(0.0003734731))*x+T(0.9959107193))*x+T(0.0149855051);
    # 8 polynomials needed.
    
    assert float(0.0000000000) <= x and x <= float(
        11.8624794162), f"Argument out-of-range. x: {x}"
    if (x < float(3.3792499610)):

        if (x < float(1.6320158198)):

            if (x < float(0.6615367791)):
                return ((float(-0.0065591595)*x+float(0.1276442762))*x+float(0.4996554598))*x+float(0.6931542306)
            return ((float(-0.0155157557)*x+float(0.1446775699))*x+float(0.4882939746))*x+float(0.6958092989)

        if (x < float(2.4912588184)):
            return ((float(-0.0128909247)*x+float(0.1301028251))*x+float(0.5150398748))*x+float(0.6795585882)
        return ((float(-0.0072142647)*x+float(0.0877540853))*x+float(0.6208708362))*x+float(0.5909675829)

    if (x < float(5.7890710412)):

        if (x < float(4.4261691294)):
            return ((float(-0.0031455354)*x+float(0.0467229449))*x+float(0.7592532310))*x+float(0.4348794399)
        return ((float(-0.0010110698)*x+float(0.0185943421))*x+float(0.8831730747))*x+float(0.2523695427)

    if (x < float(7.8162726752)):
        return ((float(-0.0001962780)*x+float(0.0046084408))*x+float(0.9634431978))*x+float(0.0983148903)
    return ((float(-0.0000113994)*x+float(0.0003734731))*x+float(0.9959107193))*x+float(0.0149855051)


# def Fast_LogExpMinusOne(x):
#     assert float(0) <= x and x <= float(30), "Argument out-of-range."
#     return math.log(math.exp(x) - float(1))


# def Fast_LogExpMinusOne(x):
#     raise Exception("Not Implemented")


def Fast_LogPlusEquals(x, y):
    # if(x!=y):
    #     print(x,y)
    if x < y:
        temp = x
        x = y
        y = temp
    if y> float(NEG_INF/2) and x-y < float(11.8624794162):
        x= Fast_LogExpPlusOne(x-y)+y
    return x

def UPDATE_MAX(bs,bt,s,t):
	if(s>bs):
		bs,bt=s,t
	return bs,bt

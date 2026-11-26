import sympy as sp
x,y,u,v,t=sp.symbols('x y u v t')
for F,(X,Y) in [(y+x*sp.tan(t)-2*sp.sin(t),(x,y)),(v+u*sp.tan(t)-3*sp.sin(t),(u,v))]:
    s=sp.solve([F,sp.diff(F,t)],[X,Y],dict=True)
    for d in s:
        px,py=sp.simplify(d[X]),sp.simplify(d[Y])
        print(px,py,'at t=1/3:',sp.N(px.subs(t,sp.Rational(1,3)),15),sp.N(py.subs(t,sp.Rational(1,3)),15))

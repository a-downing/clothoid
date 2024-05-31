var("x")

R = RealField(100)

P5 = 0 == 1/8*(63*x^5 - 70*x^3 + 15*x)
P9 = R(0) == R(1)/R(128)*(R(12155)*x^R(9) - R(25740)*x^R(7) + R(18018)*x^R(5) - R(4620)*x^R(3) + R(315)*x)
P10 = 0 == 1/256*(46189*x^10 - 109395*x^8 + 90090*x^6 - 30030*x^4 + 3465*x^2 - 63)

P5_prime = derivative(P5.right(), x)
P9_prime = derivative(P9.right(), x)
P10_prime = derivative(P10.right(), x)

P = P9
P_prime = P9_prime

solutions = solve(P, x)

print(f"Sample points:")
for solution in solutions:
    print(f"{solution.right().n(prec=100)}")

print(f"Weights:")
for solution in solutions:
    xi = solution.right()
    wi = R(2) / ((R(1) - xi^R(2))*(P_prime(x=xi))^R(2))
    print(f"{wi.n(prec=100)}")
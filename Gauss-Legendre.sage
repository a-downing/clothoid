var("x")

P5 = 0 == 1/8*(63*x^5 - 70*x^3 + 15*x)
P5_prime = derivative(P5.right(), x)

print(P5_prime)

solutions = solve(P5, x)

print(f"Sample points:")
for solution in solutions:
    print(f"{float(solution.right())}")

print(f"Weights:")
for solution in solutions:
    xi = solution.right()
    wi = 2 / ((1 - xi^2)*(P5_prime(x=xi))^2)
    print(f"{float(wi)}")
from sage.symbolic.integration.integral import definite_integral

sage: var("ke, k0, c1, c2, s1, s2, sbc, phi_e, phi_0")
R = RealField(1000)

#assume(sbc > 0)
#assume(k0 > 0)
#assume(ke > 0)

#phi_0 = 2.356194490192345
#k0 = 2
#phi_e = -2.677945044588987
#ke = 0
#sbc = 1.118033988749895

phi_0 = 0

sols = solve([
    ke == k0 + c1*(s1 - s2),
    phi_e == phi_0 + k0*(s1+s2) + 0.5*c1*(s1^2 - s2^2 + 2*s1*s2),
    sbc == s1 + s2,
    c2 == -c1
], s1, s2, c1, c2)

for sol in sols:
    for v in sol:
        print(f"{v.left()} = {v.right()}")

#f = -1/2*(2*ke*sbc + 2*phi_0 - 2*phi_e - sqrt(2*(k0^2 + ke^2)*sbc^2 + 4*phi_0^2 - 8*phi_0*phi_e + 4*phi_e^2 + 4*((k0 + ke)*phi_0 - (k0 + ke)*phi_e)*sbc))/(k0 - ke)
#g = -1/2*(2*ke*sbc + 2*phi_0 - 2*phi_e + sqrt(2*(k0^2 + ke^2)*sbc^2 + 4*phi_0^2 - 8*phi_0*phi_e + 4*phi_e^2 + 4*((k0 + ke)*phi_0 - (k0 + ke)*phi_e)*sbc))/(k0 - ke)

#h = -((k0 + ke)*sbc + 2*phi_0 - 2*phi_e + sqrt(2*(k0^2 + ke^2)*sbc^2 + 4*phi_0^2 - 8*phi_0*phi_e + 4*phi_e^2 + 4*((k0 + ke)*phi_0 - (k0 + ke)*phi_e)*sbc))/sbc^2
#i = -((k0 + ke)*sbc + 2*phi_0 - 2*phi_e - sqrt(2*(k0^2 + ke^2)*sbc^2 + 4*phi_0^2 - 8*phi_0*phi_e + 4*phi_e^2 + 4*((k0 + ke)*phi_0 - (k0 + ke)*phi_e)*sbc))/sbc^2

#limf = limit(f, ke=k0)
#limg = limit(g, ke=k0)
#print()
#print(limf)
#print(limg)

#limh = limit(h, ke=k0)
#limi = limit(i, ke=k0)
#print()
#print(limh)
#print(limi)

print("test:")
f_c1 = -2*(k0*sbc + phi_0 - phi_e - sqrt(k0^2*sbc^2 + phi_0^2 - 2*phi_0*phi_e + phi_e^2 + 2*(k0*phi_0 - k0*phi_e)*sbc))/sbc^2;
print(f_c1.simplify())

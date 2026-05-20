# Data for 'Pushforward Problems and Applications to Isogeny-based Cryptography'

This repository accompanies the paper 'Pushforward Problems and Applications to Isogeny-based Cryptography'.
In particular, it gives an implementation of the algorithms described in Section 5.5.

The file `sos.py` is imported from [sqisign_prism_v2](https://github.com/KULeuven-COSIC/sqisign_prism_v2).


### Example

This example runs on Sagemath version 10.8.
```python
from suborder_generator import (
    parameter_generation, 
    suborder_generator
)

# p = 7 (mod 8)  

# if p = 7 (mod 8), there exist integers (a,b) such that 2a^2 - b^2 = p.
while True:
    p = 8 * randint(0, 2**64) + 7
    if is_prime(p):
        break

N, e, w, d, M = parameter_generation(p)

print(f"N = {factor(N)}")
print(f"e = {factor(e)}")
print(f"w = {factor(w)}")
print(f"d = {factor(d)}")

x = ZZ(sqrt(M))

# We consider the quaternion 
# θ = w + d * x * i 
# of norm N^2 * e

Ends = suborder_generator(p, N, e, d, w, x)
B = QuaternionAlgebra(-1, -p)

θ1, θ2 = Ends.computing_endomorphisms()

basis = (B(1), θ1, θ2, θ1 * θ2)
O_max = B.maximal_order(order_basis = basis)

print(O_max)
```
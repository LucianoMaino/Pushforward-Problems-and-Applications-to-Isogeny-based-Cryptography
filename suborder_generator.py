from sage.all import (
    ZZ,
    QQ,
    ceil, 
    sqrt, 
    randint,
    QuaternionAlgebra,
    vector,
    gcd,
    is_prime,
    BinaryQF, 
    is_square
)

from sos import sum_of_squares

def parameter_generation(p):

    # we need that p = 2 * a^2 - b^2, so we have that 
    # M = (p + b^2)/2 is a square
    # 
    # We can then impose that n = (p - b^2)/2
    # so, we have M^2 - n^2 = p * b^2

    Q = BinaryQF([2, 0, -1])
    _, k = Q.solve_integer(p)
    M = ZZ((p + k**2)/2)

    # We now consider the system of equation
    # x^2 - y^2 - M z^2 = 0 
    # all the solutions are given as 

    # α = 2 * randint(0, sqrt_p) + 1
    # β = 2 * randint(0, sqrt_p) 
    # γ = 2 * randint(0, sqrt_p) + 1

    # x_ = - α**2 - β**2 - M * γ**2 + 2 * α * β
    # y_ = α**2 + β**2 - M * γ**2 - 2 * α * β
    # z_ = 2 * γ  * (β - α )

    sqrt_p = ceil(sqrt(p))

    if M%2 == 1:
        M_ = ZZ(M)
    else:
        assert M%4 == 0
        M_ = ZZ(M//4) 

    while True:
        α = 2 * randint(0, sqrt_p) + 1
        β = 2 * randint(0, sqrt_p) 
        N = - α**2 - β**2 - M_ + 2 * α * β
        w = α**2 + β**2 - M_ - 2 * α * β
        d = 2 * (β - α)

        if M%2 == 1:
            N = ZZ(N/2)
            w = ZZ(w/2)
            d = ZZ(d/2)

        if N < 0:
            N = - N

        if d < 0:
            d = - d

        #removing smooth factors
        B = 19078266889580195013601891820992757757219839668357012055907516904309700014933909014729740190
        e_ = gcd(N, B)
        e = 1
        while e_ != 1:
            N = N // e_
            e *= e_ ** 2 
            e_ = gcd(N, B)

        if M%4 == 0:
            assert d%2 == 0
            d = d // 2

        if is_prime(N) and is_prime(d) and e < sqrt_p:
            break
    
    assert N**2 * e == w**2 + d**2 * M

    return N, e, w, d, M


class suborder_generator:
    """
    This class computes two quaternions θ1, θ2 such that
        - nrd(θ1) = nrd(θ2) = N**2 * e
        - θi = w (mod d) 
        - <1, θ1, θ2, θ1*θ2> is a full-rank suborder
    """

    def __init__(self, p, N, e, d, w, x):
        #forcing everything to be a rational
        self.p = QQ(p) 
        self.N = QQ(N) 
        self.e = QQ(e) 
        self.d = QQ(d)
        self.w = QQ(w)
        self.x = QQ(x) 
        self.y = QQ(0)
        self.z = QQ(0) 
        self.M = (N ** 2 * e - w ** 2) // (d**2)


        # we know that 2M - p is a square and M is a square
        self.sqrt_M = ZZ(sqrt(self.M))
        self.k = ZZ(sqrt(2*self.M - self.p))
        self.sqrt_p = ZZ(ceil(sqrt(p)))
        self.p_one_forth = ZZ(ceil(sqrt(self.sqrt_p)))

        self.quaternion_algebra = QuaternionAlgebra(-1, -p)
        assert N ** 2 * e == w ** 2 +  d ** 2 * x**2 

    def _generate_random_x1_y1_z1(self):
        p = self.p
        sqrt_p = self.sqrt_p
        p_one_forth = self.p_one_forth
        x = self.x

        # some bounds so the factorisation is easier
        α1 = randint(0, sqrt_p)
        β1 = randint(0, p_one_forth)
        γ1 = randint(0, p_one_forth)


        t1_num = - 2 * α1 * x 
        t1_den = α1**2 + β1**2 * p + γ1**2 * p

        x1 = QQ(x + α1 * t1_num/t1_den)
        y1 = QQ(β1 * t1_num/t1_den) 
        z1 = QQ(γ1 * t1_num/t1_den)  

        return x1, y1, z1
    
    def _scalar_product(self, v1, v2):
        p = self.p
        a1, b1, c1 = v1
        a2, b2, c2 = v2
        return a1 * a2 + p * (b1 * b2 + c1 * c2)
    
    def _generate_orthogonal_basis(self, x1_y1_z1):
        p = self.p
        x1, y1, z1 = x1_y1_z1
        w1 = vector([x1, y1, z1])
        w2 = vector([-p * y1, x1, 0])
        w3_prime = vector([-p * z1, 0, x1])

        scalar_w3_prime_w2 = self._scalar_product(w3_prime, w2)
        scalar_w2_w2 = self._scalar_product(w2, w2)
        w3 = w3_prime - (scalar_w3_prime_w2 / scalar_w2_w2) * w2

        return w1, w2, w3
    
    def _solve_special_ternary_quadratic_form(self, x1_y1_z1):

        # we know that 2M - p is a square and M is a square
        p = self.p

        k = self.k
        sqrt_M = self.sqrt_M 

        x1, y1, z1 = x1_y1_z1
        z1_num = z1.numerator()
        z1_den = z1.denominator()
        # remember x1^2 + p y1^2 = M - p z1^2 = (M * z1_den^2 - p * z1_num^2)/(z1_den^2)

        # also, we have that z1_num/z1_den = - 2 * α1 * sqrt(M)/α1**2 + β1**2 * p + γ1**2 * p
        # hence M * z1_den**2 - p * z1_num**2 is divisible by M, which is a square. 

        # However, some common factors between numerator and denominator can be simplified
        t = gcd(ZZ(self.M * z1_den**2 - p * z1_num**2), self.M)
        assert is_square(t)
        sqrt_t = ZZ(sqrt(t))

        # all this machinery is to remove squares from the factorisation
    
        t2_t3 = sum_of_squares(ZZ(self.M * z1_den**2 - p * z1_num**2)//t)

        if t2_t3 == []:
            return []
        t2, t3 = t2_t3
        t2 *= sqrt_t
        t3 *= sqrt_t

        # M- n = k^2
        # back substituions
        s1 = z1_den / (k * p)
        s2 = t2 / (k**2 * sqrt_M)
        s3 = t3  * x1 / ((x1**2 + p *y1**2) * k**2)

        return s1, s2, s3 

                
    def computing_endomorphisms(self):
        
        w = self.w
        d = self.d
        flag = True

        while flag:
            x1, y1, z1 = self._generate_random_x1_y1_z1()

            s1_s2_s3 = self._solve_special_ternary_quadratic_form([x1, y1, z1])
            if s1_s2_s3:
                s1, s2, s3 = s1_s2_s3
                w1, w2, w3 = self._generate_orthogonal_basis([x1, y1, z1])
                scalar_w1_w1 = self._scalar_product(w1, w1)
                scalar_w2_w2 = self._scalar_product(w2, w2)
                scalar_w3_w3 = self._scalar_product(w3, w3)
                v = (s1 / scalar_w1_w1) * w1 + (s2 / scalar_w2_w2) * w2 + (s3 / scalar_w3_w3) * w3
                α2, β2, γ2 = v
                flag = False
                break


        p = self.p

        t2_num = - 2 * α2 * x1 - 2 * β2 * y1* p - 2 * γ2 * z1 * p
        t2_den = α2**2 + β2**2 * p + γ2**2 * p

        x2 = QQ(x1 + α2 * t2_num/t2_den)
        y2 = QQ(y1 + β2 * t2_num/t2_den) 
        z2 = QQ(z1 + γ2 * t2_num/t2_den)  
        
        B = self.quaternion_algebra
        i, j, k = B.gens()

        θ1 = w + d * (x1 * i + y1 * j  + z1 * k)
        θ2 = w + d * (x2 * i + y2 * j  + z2 * k)

        assert θ1 * θ2 != θ2 * θ1

        assert (θ1 * θ2).reduced_trace() in ZZ
        # sanity checks, always good
        N = self.N
        e = self.e
        assert θ1.reduced_norm() == N**2 * e
        assert θ2.reduced_norm() == N**2 * e
 
        return θ1, θ2
    
    def generate_maximal_order_random(self):

        θ1, θ2 = self.computing_endomorphisms()
        B = self.quaternion_algebra
        basis = (B(1), θ1, θ2, θ1 * θ2)
        O_max = B.maximal_order(order_basis = basis)

        return O_max
title: Substepping

# Particle pusher substepping

In `ALaDyn` the currently available particle pushers are the standard [_Boris_](https://books.google.it/books?hl=en&lr=&id=dQRBAQAAIAAJ&oi=fnd&pg=PA3&dq=Relativistic+plasma+simulation-optimization+of+a+hybrid+code&ots=OTEPPhsXrf&sig=nMogaw16LP2bAGu1sitHuW8pOiY#v=onepage&q=Relativistic%20plasma%20simulation-optimization%20of%20a%20hybrid%20code&f=false) (BP) or the [_Higuera-Cary_](https://aip.scitation.org/doi/abs/10.1063/1.4979989) (HC). Each of them can be chosen in the input file by selecting

```fortran
pusher = 1
```

for the HC or

```fortran
pusher = 2
```

for the BP.

## Brief description of the algorithms

The usual implementation of the particle pusher divides the particle motion in four main steps, the first three concerning the momentum update and finally the position update:

1. Particles are accelerated for an half timestep using only the electric field \(\mathbf{E}\),
2. Then, they are rotated under the influence of the \(\mathbf{v} \times \mathbf{B}\) component of the Lorentz force
3. The second half timestep only retaining the \(\mathbf{E}\) field is performed.
4. Once the new momentum has been computed, the particle position is advanced via a standard Leap-frog timestep.

This algorithm, called the _Boris-Buneman_ (BB) rotation, is very powerful since it guarantees a second order accuracy that is consistent with the _Yee_ staggered Electromagnetic field integrator usually implemented in PIC schemes.
In addition, the _Boris-Buneman_ rotation preserves the particle energy when the particle rotates in the magnetic field. In fact, it computes the approximated rotation angle (up to \(\mathcal{O}(\Delta t^3)\)) and then it performes that rotation analitically.
However, laser pulses employed in LWFA schemes, can reach very high intensities, thus requiring particles to rotate for very large angles in a single timestep, as pointed out by [Arefiev](https://aip.scitation.org/doi/pdf/10.1063/1.4965624).

## Sub-cycling

The particle pusher sub-cycling may be a clever workaround to correct the particles trajectory when increasing the whole simulation resolution is not a viable option.

In `ALaDyn` via the substepping option in the input file

```fortran
n_substeps       = 1,
```

it is possible to subdivide the rotation under the magnetic field when the particles are moved to obtain a corrected final momentum.
The BB scheme is therefore modified as following:

1. Particles are accelerated for an half timestep using only the electric field \(\mathbf{E}\),
2. Then, the magnetic field is divided by \(N=n\_substeps\),
effectively reducing the rotation angle per substep as \(\theta_i = \theta/N\),
where \(\theta\) is the rotation angle \(\theta \simeq q \mathbf{B} \Delta t/(m \gamma c)\).
The analytical rotation is performed \(N\) times taking into account the reduced angle and for every substep the intermediate step is cached in a temporary variable that then defines the new momentum to rotate,
3. After \(N\) sub-rotations have been performed, momentum is updated with the second half push under \(\mathbf{E}\).
4. Particle position is advanced via a standard Leap-frog timestep.

The particle substepping is implemented in the `boris_push` module as follows

```fortran
do nss = 1, nstep
    pxp(1:np) = ((gam02(1:np) - b2(1:np))*pxpo(1:np) + 2*(bv(1:np)*bb(1:np, 1) + &
     gam(1:np)*(pypo(1:np)*bb(1:np, 3) - pzpo(1:np)*bb(1:np, 2))))/&
     (gam02(1:np) + b2(1:np))
    pyp(1:np) = ((gam02(1:np) - b2(1:np))*pypo(1:np) + 2*(bv(1:np)*bb(1:np, 2) + &
     gam(1:np)*(pzpo(1:np)*bb(1:np, 1) - pxpo(1:np)*bb(1:np, 3))))/&
     (gam02(1:np) + b2(1:np))
    pzp(1:np) = ((gam02(1:np) - b2(1:np))*pzpo(1:np) + 2*(bv(1:np)*bb(1:np, 3) + &
     gam(1:np)*(pxpo(1:np)*bb(1:np, 2) - pypo(1:np)*bb(1:np, 1))))/&
     (gam02(1:np) + b2(1:np))
    pxpo(1:np) = pxp(1:np)
    pypo(1:np) = pyp(1:np)
    pzpo(1:np) = pzp(1:np)
end do
```

where for every subcycle `pxpo`, `pypo` and `pzpo` are the initial momenta before the rotation,
the `b2` vector caches \(|\mathbf{B}|^2\),
`bb` is the magnetic field \(\mathbf{B}\),
`bv` contains the scalar product \(\mathbf{B}\cdot \mathbf{u}^-\), with \(\mathbf{u}^-\) the updated momentum after the step 1. of the pusher scheme.

The `gam` and `gam02` vectors contain the Lorenz factor and the squared Lorentz factor respectively, where that could be defined as in the Boris push
$$\gamma_B^2=1+|\mathbf{u}^-|^2,$$
or as in the Higuera-Cary one
$$\gamma_H^2=\frac{1}{2}\left(\gamma_B^2 - \beta^2 + \sqrt{\left(\gamma_B^2 - \beta^2\right)^2+4\left(\beta^2+|\beta \cdot \mathbf{u}^-|\right)}\right),$$
where \(\beta=\gamma \theta/2 = q \mathbf{B}\Delta t/2mc\).

We note that since the BB rotation preserves the norm of the momentum and the scalar product
\(|\beta \cdot \mathbf{u}|\), the Lorentz factor can be computed just once before the rotation regardless te number of timesteps employed.
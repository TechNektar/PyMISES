<a name="br1"></a> 

A User’s Guide to MISES 2.63

Mark Drela, Harold Youngren

MIT Aerospace Computational Design Laboratory

February 2008

This is a user manual for the MISES viscous/inviscid cascade analysis and design system.



<a name="br2"></a> 

Contents

1 Overview

4

4

5

2 Internal Reference Quantities

3 Streamsurface and Blade geometry deﬁnition

4 Input Files

6

7

9

4\.1 Blade coordinate ﬁle blade.xxx . . . . . . . . . . . . . . . . . . . . . . . . . . .

4\.2 Geometry parameter ﬁle bparm.xxx . . . . . . . . . . . . . . . . . . . . . . . .

4\.3 Modiﬁed-geometry parameter ﬁle bspec.xxx . . . . . . . . . . . . . . . . . . . 10

4\.4 Geometry parameter speciﬁcation ﬁle bplist.xxx . . . . . . . . . . . . . . . . 10

4\.5 Stream surface ﬁle stream.xxx . . . . . . . . . . . . . . . . . . . . . . . . . . . 11

4\.6 Prescribed-loss ﬁle loss.xxx . . . . . . . . . . . . . . . . . . . . . . . . . . . . 12

4\.7 Wall-suction speciﬁcation ﬁle suct.xxx . . . . . . . . . . . . . . . . . . . . . . . 13

4\.8 Flow condition ﬁle ises.xxx . . . . . . . . . . . . . . . . . . . . . . . . . . . . 14

4\.8.1 Variable,Constraint indices . . . . . . . . . . . . . . . . . . . . . . . . . . 21

4\.8.2 Pressure-correction term . . . . . . . . . . . . . . . . . . . . . . . . . . . . 22

4\.8.3 Momentum/Entropy conservation . . . . . . . . . . . . . . . . . . . . . . 22

4\.8.4 Artiﬁcial dissipation . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 23

4\.8.5 Artiﬁcial dissipation level selection . . . . . . . . . . . . . . . . . . . . . . 25

4\.8.6 Dissipation enhancement during convergence . . . . . . . . . . . . . . . . 27

4\.9 Example ises.xxx input-ﬁle lines . . . . . . . . . . . . . . . . . . . . . . . . . 28

4\.9.1 Lines 1–4. Variables, constraints, ﬂow conditions. . . . . . . . . . . . . . . 28

4\.9.2 Lines 6–7. Viscous ﬂow parameters. . . . . . . . . . . . . . . . . . . . . . 30

4\.9.3 Line 8. Isentropy and dissipation . . . . . . . . . . . . . . . . . . . . . . . 31

4\.9.4 Line 9. Streamtube thickness mode amplitudes . . . . . . . . . . . . . . . 32

4\.10 Geometry perturbation mode speciﬁcation ﬁle modes.xxx . . . . . . . . . . . . 32

4\.11 Design-parameter speciﬁcation ﬁle params.xxx . . . . . . . . . . . . . . . . . . 33

1



<a name="br3"></a> 

5 Program Descriptions

34

5\.1 ISET . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 35

5\.1.1 Basic Initialization . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 35

5\.1.2 Panel solution . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 35

5\.1.3 Initial surface gridding . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 36

5\.1.4 Grid smoothing . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 38

5\.1.5 Initial solution ﬁle output . . . . . . . . . . . . . . . . . . . . . . . . . . . 38

5\.1.6 Grid parameters . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 38

5\.1.7 Grid parameter saving, recall . . . . . . . . . . . . . . . . . . . . . . . . . 39

5\.1.8 Smoothing and writing the grid . . . . . . . . . . . . . . . . . . . . . . . . 39

5\.2 ISES . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 40

5\.2.1 Inﬂow boundary conditions . . . . . . . . . . . . . . . . . . . . . . . . . . 40

5\.2.2 Outﬂow boundary conditions . . . . . . . . . . . . . . . . . . . . . . . . . 42

5\.3 IPLOT . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 43

5\.3.1 Blade surface plots . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 43

5\.3.2 Suction . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 45

5\.3.3 Streamtube plots . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 45

5\.3.4 Contour/grid plots . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 47

5\.3.5 Wake proﬁle plots . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 47

5\.3.6 r,b vs m<sup>′</sup> stream surface deﬁnition plots . . . . . . . . . . . . . . . . . . . 47

5\.3.7 Wheel view . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 47

5\.4 EDP . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 47

5\.4.1 Surface parameterization . . . . . . . . . . . . . . . . . . . . . . . . . . . 48

5\.4.2 EDP execution . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 49

5\.4.3 Modal-Inverse . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 52

5\.4.4 Parametric-Inverse . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 53

5\.4.5 Blade Translation, Scaling, Rotation . . . . . . . . . . . . . . . . . . . . . 53

5\.4.6 Modiﬁed-Blade Output . . . . . . . . . . . . . . . . . . . . . . . . . . . . 53

5\.4.7 ISES Parameter Changes . . . . . . . . . . . . . . . . . . . . . . . . . . . 54

2



<a name="br4"></a> 

5\.4.8 Inverse Design Session . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 54

5\.4.9 Parameter-Modiﬁcation Design Session . . . . . . . . . . . . . . . . . . . . 55

5\.5 POLAR . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 56

5\.6 BLDSET . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 58

6 Optimization

7 Graphics

59

59

59

8 General Hints

8\.1 Viscous solutions . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 59

8\.2 Inverse solutions . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 60

8\.3 Grid resolution . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 60

8\.4 Execution times . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 61

3



<a name="br5"></a> 

1 Overview

The MISES system is a collection of programs for cascade analysis and design. This includes

programs for grid generation and initialization, ﬂow analysis, plotting and interpretation of

results, and an interactive program to specify design conditions.

The block diagram for these programs is given at the end of this manual. The basic grid and

ﬂow data ﬁle for a case is the so-called state ﬁle named idat.xxx, where “xxx” is an extension

suﬃx used to designate the case being run. The state ﬁle is initialized using ISET from the

blade geometry ﬁle blade.xxx and the optional stream surface geometry ﬁle stream.xxx and

the prescribed loss schedule ﬁle loss.xxx. The ﬂow solver ISES uses the state ﬁle and a ﬂow

condition ﬁle ises.xxx that speciﬁes the ﬂow conditions and program conﬁguration ﬂags. The

POLAR program performs the same calculations as ISES , but for a set of speciﬁed parameters.

Additional design condition information can be interactively added to the state ﬁle using the

EDP pressure edit program. The IPLOT program plots the ﬂow and geometry data from the

state ﬁle in an interactive plotting session.

2 Internal Reference Quantities

All ﬂow variables used by MISES are deﬁned in the relative frame. Internally, MISES employs

rotation-corrected stagnation density and speed of sound, ρ<sub>oa</sub>, a<sub>oa</sub>, as the basic reference ﬂow

variables, so that ρ<sub>oa</sub> = 1 and a<sub>oa</sub> = 1 by deﬁnition. The corresponding rotation-corrected

stagnation pressure p<sub>oa</sub> and enthalpy h<sub>oa</sub> ≡ I (i.e. rothalpy) are then related as follows.

2

γ p<sub>oa</sub> = ρ<sub>oa</sub>

a

= (γ − 1) ρ<sub>oa</sub>

I

oa

The Fortran names and assigned values of the internal reference quantities are ρ<sub>oa</sub> = RSTR0 = 1,

and I = HSTR0 = 1/(γ −1). Normally, these are not of concern for the user, since all input

and output is typically done via ratios and related dimensionless quantities, following common

conventions. For example, the outlet pressure is speciﬁed as p /p , where p is the conventional

2

o1

o1

relative-frame total pressure at the inlet at the speciﬁed radius r<sub>1</sub>.

The ( )<sub>oa</sub> notation means an “absolute” total quantity (not to be confused with an absolute-

frame total quantity), in the sense that it implies an isentropic process where the ﬂuid is brought

to rest in the relative frame, and taken to the rotation center r = 0. Bringing the ﬂuid to rest

in the relative frame at a ﬁxed radius r gives the conventional stagnation quantities ρ , a ,

o

o

etc. The absolute and conventional stagnation quantities are related by the usual isentropic

relations.



!

γ

2

2

I + Ω r /2

γ−1

p<sub>o</sub>

\=

p<sub>oa</sub>

I



!

1

γ−1

2

2

I + Ω r /2

ρ<sub>o</sub>

ρ<sub>oa</sub>

2

\=

\=

I

2

I + Ω r /2

2

a

o

a <sup>2</sup>

oa

I

4

![ref1]![ref2]![ref3]![ref1]![ref2]![ref3]![ref1]![ref2]

<a name="br6"></a> 

The reason ρ<sub>oa</sub> and a<sub>oa</sub> were chosen for the internal reference quantities is precisely because they

are independent of radius, and thus considerably simplify the internal formulation of the code.

But again, they are transparent to the user, and are only described here in case source-code

additions are being contemplated.

3 Streamsurface and Blade geometry deﬁnition

The blade airfoil and grid domain geometry are deﬁned in the standard m<sup>′</sup> − θ streamsurface

coordinate system, shown in Figure 1. With z denoting the cylindrical axis coordinate and r

the local streamsurface radius, the m<sup>′</sup> coordinate is deﬁned by

Z

Z √

dm

r

dr<sup>2</sup> + dz<sup>2</sup>

m<sup>′</sup> =

\=

r

while θ is the usual circumferential angle. The total arc length increment ds in the stream

surface is given by

ds = <sup>q</sup>dr<sup>2</sup> + dz<sup>2</sup> + (r dθ)<sup>2</sup> = <sup>q</sup>dm<sup>2</sup> + (r dθ)<sup>2</sup> = r<sup>p</sup>dm<sup>′2</sup> + dθ<sup>2</sup> = r ds<sup>′</sup>.

√

The normalized arc length ds<sup>′</sup> = dm<sup>′2</sup> + dθ<sup>2</sup> will be used later as a spline parameter to

√

deﬁne the blade shape in the m<sup>′</sup> − θ plane. The intermediate coordinate dm = dr<sup>2</sup> + dz<sup>2</sup> is

the physical arc length increment projected onto the meridional r − z plane, and is not used

explicitly in MISES.

The transformation from physical space to the m<sup>′</sup> − θ plane is angle-preserving. Hence,

no transformation is required for ﬂow angles or surface normal directions. This simpliﬁes

imposition of boundary conditions such as a speciﬁed inlet ﬂow angle, or the normal-oﬀsetting

of the inviscid ﬂow by the viscous displacement thickness.

For 2-D cascades, r becomes an arbitrary constant scaling length, and hence

z

m<sup>′</sup> =

(2-D cascade).

r

For a purely radial cascade (e.g. squirrel cage fan) z is a constant, and in this case

Z

dr

m<sup>′</sup>

m<sup>′</sup>

\=

\=

= ln r

(radial outﬂow cascade)

(radial inﬂow cascade)

r

Z

−dr

= ln 1/r

r

For some other analytically-deﬁned r(z) distributions, the general m<sup>′</sup> deﬁnition might also

be integrable in closed form. However, for a general blade section “slice” deﬁned in discrete

Cartesian x y z coordinates, numerical integration for the corresponding discrete m<sup>′</sup> θ is

i

i

i

i

i

necessary. A simple trapezoidal integration is appropriate.

θ<sub>i</sub> = arctan <sup>ꢀ</sup>

ꢁ

y<sub>i</sub>

x<sub>i</sub>

5

![ref1]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAD8DASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAf/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAEsDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAf/Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAHoDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAAAAAAAf/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAFkDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAAAf/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAEcDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAf/Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAEcDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAf/Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAD8DASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAf/9k=)![ref4]![ref5]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAB0DASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAf/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADABADASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAf/Z)

<a name="br7"></a> 

*r*

*V*

*2*

<i>S<sub>2</sub></i>

*b*

*m*

*V*

*1*

*r*

<i>S<sub>1</sub></i>

*m’*

*z*

*m’*

*1*

*m’*

*2*

Figure 1: Streamsurface deﬁnition.

q

2

2

r<sub>i</sub>

\=

x + y

i

i

<sup>q</sup>(r<sub>i</sub> − r<sub>i−1</sub>)<sup>2</sup> + (z<sub>i</sub> − z<sub>i−1</sub>)<sup>2</sup>

2

m<sup>′</sup><sub>i</sub> = m<sup>′</sup>

\+

i−1

<sup>r</sup>i <sup>+ r</sup>

i−1

If the discrete integration method is used, it may be necessary to slightly adjust m<sup>′</sup> θ at the

i

i

i = 1 and i = N points to obtain the exact trailing edge gap desired. Even if these two points

have identical x, y, z coordinates, as in the case of a sharp trailing edge, they will not generally

have identical m<sup>′</sup>, θ coordinates due to numerical integration errors. In general, these errors

will be quite small if the original x, y, z coordinates are reasonably dense. The initial m<sup>′</sup>

i=1

coordinate is arbitrary, and merely shifts all the succeeding m<sup>′</sup><sub>i</sub> values. Likewise, an arbitrary

constant can be added to all θ values. These shifts are best done after the integration to position

the blade “in space” where desired. Placing the leading edge near the origin is convenient.

4 Input Files

The principal input ﬁles needed to set up a MISES solution case are listed below. They are

described in more detail in subsequent sections. The “xxx” is an arbitrary extension suﬃx used

to designate the case being run.

blade.xxx Required unless bparm.xxx is used. Deﬁnes blade shape via a relatively large

number of coordinate pairs.

bparm.xxx Required unless blade.xxx is used. Alternative way to deﬁne a blade shape via

an arbitrary set of parameters, interpreted by user-supplied routines.

bspec.xxx Optional — used for redesign cases only. Deﬁnes a new blade shape to be

imposed on an existing ﬂow solution after modest geometry-parameter changes

(this is much more eﬃcient than starting a new solution). This ﬁle must have

the same format as the bparm.xxx ﬁle.

6

![](data:image/jpeg;base64,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)![](data:image/jpeg;base64,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)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAOAAkDASIAAhEBAxEB/8QAFwAAAwEAAAAAAAAAAAAAAAAABQgJCv/EACcQAAEEAgADCQEAAAAAAAAAAAYDBAUHAggAAQkTFRYZVVmVl9XW/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL1bx9RSDr/as91zINuQ/TZjSVegVsBhETxsgUN75ss8cEiMKFSGTF1DtWQIAZDKzM5B3D1dUudHIy4ylB/kP9nLL35xfV79kO1vtFL+B4GWlJxj/c/ZTo7OWHOVidprzF9v7LstwkjFr5UdsC8nVLYpwezRyk5EfMRl7XoWoEmCGDtkRISU/hLMhPnHM++NJngiJ9MS+SW/O4D/2Q==)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAANABUDASIAAhEBAxEB/8QAGAAAAwEBAAAAAAAAAAAAAAAAAAgJBgr/xAAgEAACAgICAwEBAAAAAAAAAAAEBgUHAgMBCAAJFhUS/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AOm64fbvQ1L91K26lSYrmSBOJl4SVhswtQ3/ACJie5VpL1dGq68rREPVMiA/wbVg7Tu49rWZKQgFfmADxkD9X7ofOzQux1gdwO2XYHq5nbduUXV3WNQo1oIKoRvNQrJtBnvSEbJuJkpexAdg8mrKtfhJ8vFkp44snpdimePmj9kYQpDDmUYmaVq9gtdGvGZUgzrWrVXe0tHcdhkpgavrFmEqpjzEDB6T9cURoYCUhW279pwBRIvMRq4j94mO8zEicTjXr5t9phdxgvi5CKSZ06sBJEruHRJYCQb2hwn60YI93tBx0WHiA7/C/FkxCTD8JEMdBQ7YyhjM2OMhuy2AlV3dnvYN0ktGe68dA6Be/awAsxy9OW3ldF7Br9lddJ5rjcZhZXZm0nnkcu24yx4zbIsi/wAh6cNSJGxGMBu4/snDw8cD1IKNpLKv2SNuB+RLRZrAvVktbNnWKvna8kws7JmWNrKTztzFbFsGTCqmbpDiCrkDg+N+WWB8IjP9PnngrA8D/9k=)![ref6]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADADkDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAf/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAKoDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAAAAAAAAAAAAAf/Z)

<a name="br8"></a> 

bplist.xxx Optional — used in conjunction with bspec.xxx to specify which parameters

are to be modiﬁed.

stream.xxx Required except for strictly-2D cases. Speciﬁes streamsurface radius and thick-

ness distributions.

loss.xxx Optional. Deﬁnes an additional loss distribution imposed on the solution.

suct.xxx Optional. Speciﬁes one or more wall-suction locations.

ises.xxx Required. Speciﬁes the ﬂow conditions and solver control ﬂags.

modes.xxx Optional. Speciﬁes the blade geometry perturbation mode shapes.

params.xxx Optional. Speciﬁes new design-parameter values which overrule those speciﬁed

in the ises.xxx ﬁle.

4\.1 Blade coordinate ﬁle blade.xxx

The discrete blade airfoil geometry coordinates, however generated, are speciﬁed in the for-

matted ﬁle blade.xxx, which is used by the initialization program ISET to deﬁne the initial

streamline grid. This ﬁle has the following structure.

NAME

SINL SOUT CHINL CHOUT PITCH

X(1) Y(1)

X(2) Y(2)

←− Blade 1

.

.

.

.

X(I) Y(I)

999\.

←− Start new blade deﬁnition

X(1) Y(1)

X(2) Y(2)

X(3) Y(3)

←− Blade 2

.

.

.

.

X(I) Y(I)

NAME

SINL

is the name of the case, not more than 32 characters.

is the initial S = tan(β ), the tangent of the inlet ﬂow angle relative to the axial

1

1

direction. This is the default inlet ﬂow slope at ISET startup, and can be changed

interactively.

SOUT

is the initial S = tan(β ). With the new panel-solver grid generator in ISET ,

2

2

SOUT is no longer used.

7



<a name="br9"></a> 

CHINL is the distance in m<sup>′</sup> from the blade 1 leading edge to the grid inﬂow plane.

CHOUT is the distance in m<sup>′</sup> from the blade 1 trailing edge to the grid outﬂow plane.

PITCH is the circumferential pitch of the cascade in radians = 2π/(number of blades)

The speciﬁed inlet ﬂow slope SINL and the outlet ﬂow slope SOUT calculated in ISET are

merely the initial values for the ﬂow slope variables S S , which will usually change during

1

2

the ﬂow calculation. For a well-behaved solution process, SINL should be comparable to the

ﬁnal inlet ﬂow slope S<sub>1</sub> expected in the ﬂow, although in practice there is quite a lot of room

for “error” in the initial value, especially for subsonic ﬂows. For supersonic inlet ﬂows, SINL

is preferably quite close to the ﬁnal expected S<sub>1</sub>, since the ﬂow is then very sensitive to small

inlet ﬂow angle changes.

The grid inﬂow and outﬂow plane locations speciﬁed by CHINL and CHOUT are not the locations

of the inlet and outlet condition-deﬁning planes m and m shown in Figure 1. The latter are

1

2

deﬁned separately in the ﬂow condition input ﬁle to be described later. It is only necessary that

CHINL and CHOUT be large enough so that m and m fall inside the grid. The only drawback

1

2

to making them too large is that a larger grid results, which produces correspondingly longer

run times. The run time scales very nearly linearly with the number of streamwise grid points.

The blade coordinates X(1),Y(1) through to X(I),Y(I) are the m<sup>′</sup>, θ coordinates of the

blade surface, starting at the trailing edge, going round the leading edge in either direction,

then going back to the trailing edge. For multiple blades the ﬁrst blade deﬁnition is ended

with a 999.0 coordinate. The second blade deﬁnition follows the same format as the ﬁrst and

is assumed to be in the same coordinates matching the main blade. If necessary, it will be

relocated so that it lies within the passage formed by the main blade and its periodic image

at Y - PITCH. Note that this diﬀers from previous MISES conventions. The new convention is

necessary for the new ISET grid generator.

A blunt trailing edge is speciﬁed by leaving the blade “open”, so that the ﬁrst and last coor-

dinate points do not cooincide. If the actual blade has the common semi-circular trailing edge,

it must be “cut oﬀ” near the tangency points. The Kutta condition is applied between these

two points. ISES incorporates a blunt trailing edge model which accounts for the additional

losses associated with a blunt trailing edge. Recent investigations indicate that this model

underpredicts the losses if signiﬁcant vortex shedding is present, since the additional Reynolds

stresses associated with the vortex motion are not represented by the standard turbulent-wake

formulation. A future MISES version may address this deﬁciency.

A sharp leading edge may be speciﬁed for either or both blades by repeating the leading edge

coordinates. However, sharp leading edges may lead to problems with singular velocities for

subsonic or transonic cascades. This is not a problem at supersonic Mach numbers where the

“unique incidence” eﬀect aligns the ﬂow with the leading edge.

The ﬁrst two lines in blade.xxx can be omitted, in which case the user will be prompted to

enter the missing information from the keyboard.

8



<a name="br10"></a> 

For 2-D cases, X, Y are just the Cartesian coordinates, in units of some arbitrary reference

length L<sub>ref</sub> which is also used to deﬁne the Reynolds number and the rotational wheel speed,

discussed below.

Inside all the programs, the blade shape is deﬁned analytically as a parametric cubic spline

of the form m<sup>′</sup>(s<sup>′</sup>), θ(s<sup>′</sup>), where the spline parameter s<sup>′</sup> is the arc length in the m<sup>′</sup>–θ plane.

Z

Z

<sup>p</sup>dm<sup>′2</sup> + dθ<sup>2</sup>

ds

r

s<sup>′</sup> =

\=

The spline parameter s<sup>′</sup> and the necessary spline derivatives dm<sup>′</sup>/ds<sup>′</sup> and dθ/ds<sup>′</sup> are calculated

directly from the input X, Y coordinates, which are used as the spline knots where continuity

of second derivatives is enforced. Constant-curvature end conditions (zero third derivative) are

used.

4\.2 Geometry parameter ﬁle bparm.xxx

MISES supports the ability to deﬁne the blade geometry via an arbitrary set of parameters,

G , k = 1 . . . K. Examples of G might be Bezier coeﬃcients, complex-mapping coeﬃcients, or

k

k

some set of camber and thickness modes. These parameters are speciﬁed in the ﬁle bparm.xxx,

which replaces the raw-coordinate blade.xxx ﬁle described above.

The following routines, located in the src/geo/ directory, are programmed by the user to

read, write, and process the information in the bparm.xxx ﬁle:

Routine source ﬁle Purpose

BPREAD bpario.f Reads G<sub>k</sub> from ﬁle bparm.xxx

BPWRIT bpario.f Writes G<sub>k</sub> to ﬁle bparm.xxx

BLDGEN bldgen.f Generates m<sup>′</sup>, θ points, and optionally ∂m<sup>′</sup>/∂G , ∂θ/∂G , from G

k

k

k

BPCON bpcon.f Sets constraint residuals R(G ) and also ∂R/∂G .

k

k

The geometry sensitivities ∂m<sup>′</sup>/∂G , ∂θ/∂G allow Parametric-Inverse calculations to be per-

k

k

formed, with ISES determining the combination of G<sub>k</sub> which best matches a speciﬁed pressure

distribution in a least-squares sense. This is a generalization of the Modal-Inverse method em-

ployed by all previous MISES versions. If Parametric-Inverse calculations are not required, then

all ∂( )/∂G values can be returned as zero. The residuals R(G ) embody constraints which can

k

k

be imposed on the geometry in Parametric-Inverse calculations. Neither R(G<sub>k</sub>) nor ∂R/∂G<sub>k</sub>

are needed if the constraints are not used.

As can be seen from the call lists of BPREAD and BPWRIT, these subroutines must return

and accept the case name and grid parameters — the same information which is speciﬁed in

the ﬁrst two lines of the blade.xxx ﬁle. Hence, it is natural for the bparm.xxx ﬁle to have the

form shown below, although any format is admissible as long as it is understood by BPREAD

and BPWRIT.

NAME

9

![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADABEDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAf/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAEcDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAf/Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAlsDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAf/9k=)

<a name="br11"></a> 

SINL SOUT CHINL CHOUT PITCH

< geometry parameters >

The geometry parameters read and returned by BPREAD will be passed to BLDGEN for

conversion to m<sup>′</sup>, θ coordinates. Subsequent processing is then essentially the same as though

these m<sup>′</sup>, θ coordinates were speciﬁed via blade.xxx.

BLDGEN itself can of course call other private routines as needed; its internal operations

are of no consequence to MISES. For normal analysis calculations, BLDGEN will be called

only once by ISET . For Parametric-Inverse calculations, it will also be called repeatedly by

ISES to regenerate the blade geometry if the G<sub>k</sub> parameters have been modiﬁed during a

Newton iteration.

The BPCON routine deﬁnes user-supplied constraints on the geometry parameters. Examples

might be constraints on the blade cross-sectional area, bending inertia, trailing-edge angle, etc.

These can be imposed in Parametric-Inverse calculations (described later).

4\.3 Modiﬁed-geometry parameter ﬁle bspec.xxx

This has exactly the same format as ﬁle bparm.xxx, and contains modiﬁed parameter values

which are to be imposed on a solution in redesign cases.

4\.4 Geometry parameter speciﬁcation ﬁle bplist.xxx

This ﬁle lists the parameters which are to be modiﬁed, either by driving them to the modiﬁed

values in bspec.xxx, or by performing a Parametric-Inverse calculation. It has the following

format.

1

! Parameter\_1\_name

-3 ! Parameter\_3\_name

4

8

.

.

K

! Parameter\_4\_name

! Parameter\_8\_name

! Parameter\_K\_name

Only parameters whose indices appear in this list are treated as global variables to be modiﬁed.

All other parameters are held ﬁxed at their current values. The parameter indices correspond

to the G(k) array indices “k” in the SUBROUTINE BPREAD, BPWRIT, BLDGEN, BPCON

call lists. A negative index is ignored, so that parameter 3 in the example above is eﬀectively

absent. Also, all input after the “!” is currently ignored, and the parameter-name strings

appear only for the user’s convenience. Future MISES versions will likely read the parameter

names for use in interactive menus.

10



<a name="br12"></a> 

4\.5 Stream surface ﬁle stream.xxx

This is an optional formatted ﬁle which speciﬁes the radius and thickness r(m<sup>′</sup>) and b(m<sup>′</sup>) of

the stream surface on which the blade-to-blade ﬂow is calculated (see Figure 1). It is used by

the initialization program ISET to deﬁne the ﬂow domain. The actual streamtube thickness

used in MISES is deﬁned as

b(m<sup>′</sup>) = b (m<sup>′</sup>) + B b (m<sup>′</sup>) + B b (m<sup>′</sup>)

0

1

1

2

2

where b is the baseline thickness distribution, and b b are optional modiﬁcation modes con-

0

1

2

trolled by the mode amplitudes B<sub>1</sub> B<sub>2</sub> (Fortran names: BVRN(1) BVRN(2) ). The purpose of

the two modiﬁcation terms is to more easily allow adjustment of the streamtube thickness to

account for eﬀects such as endwall losses and/or cooling mass ﬂow addition. The mode ampli-

tudes B B can be set implicitly by the ﬂow solver to attain a speciﬁed mass ﬂow or outlet

1

2

pressure, for example. This will be described later. Note: MISES 2.5 is actually coded for an

arbitrary number of b<sub>k</sub>(m<sup>′</sup>) modes, although the current array dimension is set as PARAMETER

(NBVRX = 2) in the src/INC/STATE.INC ﬁle.

The stream.xxx ﬁle has the structure

ROTREL

X(1) R(1) B 0(1) [ B 1(1) B 2(1) ] ←− optional

X(2) R(2) B 0(2) [ B 1(2) B 2(2) ]

.

.

.

.

.

.

.

.

.

.

X(I) R(I) B 0(I) [ B 1(I) B 2(I) ]

which speciﬁes the m<sup>′</sup>, r, b , b , b values. The last two columns can be omitted.

0

1

2

ROTREL is the non-dimensionalized wheel speed

ROTREL = Ω L<sub>ref</sub>/a<sub>o1</sub>

with Ω positive for the blade row moving “up” in the positive θ = Y direction as shown in

Figure 1. The normalizing quantities are the previously-described reference length L<sub>ref</sub> and the

relative-frame inlet stagnation speed of sound a<sub>o1</sub>.

The X values (m<sup>′</sup> coordinates) are in the same set of axes as those used to deﬁne the blade

airfoil geometry X values in blade.xxx. The stream surface radii R are in units of the same

reference length L<sub>ref</sub> used to deﬁne the wheel speed ROTREL:

r

R =

<sup>L</sup>ref

The streamtube thickness modes B 0, B 1, B 2, can have any length units (only the percentage-

wise changes in b are signiﬁcant). If B 1 and/or B 2 are omitted or all zeros, then default b<sub>1</sub>

and/or b<sub>2</sub> distributions will be calculated using SUBROUTINE BHDEF (in src/ises/rbcalc.f).

11

![ref7]![ref7]![ref7]![ref7]![ref7]![ref7]![ref7]![ref7]![ref7]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADABsDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAf/9k=)![ref7]![ref7]![ref7]![ref7]![ref7]

<a name="br13"></a> 

The current default distributions are

1

m<sup>′</sup> − m<sup>′</sup>

1

b<sub>1</sub>

\=

[1 − cos(πt)]

;

;

t =

2

m<sup>′</sup> − m<sup>′</sup>

2

1

1

1

2

m<sup>′</sup> − (m<sup>′</sup> + m<sup>′</sup> )

2

1

2

b<sub>2</sub>

\=

[1 − cos(πt)]

t =

1

2

m<sup>′</sup> − (m<sup>′</sup> + m<sup>′</sup> )

2

1

2

where m<sup>′</sup> and m<sup>′</sup> are the inlet and outlet condition speciﬁcation locations as described below.

1

2

If t falls outside the range 0 ≤ t ≤ 1, it is reset to the endpoint value, so that b and b are

1

2

conveniently either zero or unity outside this range, although in practice the mode shapes are

quite arbitrary. Care must be taken in deﬁning and using these modes so that the total b(m<sup>′</sup>)

does not become negative anywhere in the grid domain!

The stream surface R and B deﬁnitions are splined in X to generate intermediate values and

derivatives. A repeated coordinate pair can be used to specify a slope discontinuity in the

deﬁnitions, although this is not realistic and not generally advised. To avoid extrapolating past

the spline endpoints, the stream surface should be deﬁned so that the domain inﬂow boundaries

(typically Xle-CHINL to Xte+CHOUT) fall well inside the spline X limits. The splined r(m<sup>′</sup>)

amd b(m<sup>′</sup>) distributions can be plotted in IPLOT after the initial grid is generated. Any b<sub>1</sub>(m<sup>′</sup>)

and/or b<sub>2</sub>(m<sup>′</sup>) mode contributions are shown as well.

4\.6 Prescribed-loss ﬁle loss.xxx

This is an optional formatted ﬁle which speciﬁes the total pressure losses added to the ﬂowﬁeld

conservation equations. The resulting streamwise momentum equation has the form

∆ (ln p<sub>oa</sub>) = ∆P

where ∆( ) implies a change over a cell along a streamtube. With a constant rothalpy being

assumed, the righthand side forcing term is in eﬀect a prescribed entropy change

∆S

∆P = ∆ (ln p<sub>oa</sub>)<sub>prescribed</sub> = −

R

which is intended as a means of modeling endwall losses, for example. Naturally, P should

monotonically decrease for any physical loss-generating process.

The prescribed distribution of P(m<sup>′</sup>) is given in the optional ﬁle loss.xxx, which has the

following structure.

X(1) P(1)

X(2) P(2)

.

.

.

.

X(I) P(I)

The X values are m<sup>′</sup> coordinates in the same set of axes as those used to deﬁne the blade airfoil

geometry X values in blade.xxx. The P values are P as deﬁned above. Because only changes

in P are used in the solution, it can have an arbitrary additive constant.

12

![ref4]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADADsDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAf/9k=)![ref4]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAAgDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAf/Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAHQDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAAAAAAAf/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAAgDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAf/Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADABkDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAf/9k=)

<a name="br14"></a> 

As with the stream.xxx ﬁle, the X value range must contain the entire grid m<sup>′</sup> range. The

prescribed loss at any m<sup>′</sup> location is obtained from the spline representation P(m<sup>′</sup>) just like the

other streamsurface quantities r(m<sup>′</sup>) and b(m<sup>′</sup>).

4\.7 Wall-suction speciﬁcation ﬁle suct.xxx

The ﬁle deﬁnes the location and strength of one or more suction slots on the blade surface(s).

CQsuct(1) Sbeg(1) Send(1) Iside(1)

CQsuct(2) Sbeg(2) Send(2) Iside(2)

.

.

Each line speciﬁes the parameters for one suction slot:

CQsuct = suction mass-ﬂow coeﬃcient C<sub>Q</sub>

Sbeg = fractional arc length location s<sup>′</sup><sub>beg</sub>/s<sub>s</sub><sup>′</sup><sub>ide</sub> of front of slot

Send = fractional arc length location s<sup>′</sup><sub>end</sub>/s<sub>s</sub><sup>′</sup><sub>ide</sub> of rear of slot

Iside = blade side containing suction slot

The suction coeﬃcient is the ratio of the suction to total mass ﬂows. These are deﬁned as

Z

m˙ <sub>suct</sub>

s<sup>′</sup>

end

C<sub>Q</sub>

\=

m˙ <sub>suct</sub>

\=

−ρ v b r ds<sup>′</sup>

m˙ = ρ u b r P

1 1 1

w

w

1

m˙

s<sup>′</sup>

beg

where P is the pitch in radians, and ()<sub>1</sub> is a reference inlet quantity. The distribution of the

mass ﬂux ρ v over the extent of the slot is assumed to be parabolic in s<sup>′</sup>, with a maximum in

w

w

the middle of the slot. A uniform distribution is also implemented, but is currently commented

out in the deﬁning SUBROUTINE SETSUCTION (in isesinit.f).

The slot extends over the fractional arc length region Sbeg < s<sup>′</sup>/s<sup>′</sup> < Send, where s<sup>′</sup> is

side

side

the distance in the m<sup>′</sup>-θ plane from the stagnation point to the trailing edge. Note that this

will make the slot move slightly along the surface as the stagnation point moves. This was done

for implementation simplicity, since the slot is then “ﬁxed” on the sliding grid. Future versions

may allow speciﬁcation of the slot at a ﬁxed physical location.

The slot is located on blade side Iside. The sides are numbered as follows.

Side 1: Blade 1 top

Side 2: Blade 1 bottom

Side 3: Blade 2 top

Side 4: Blade 2 bottom

13

![ref8]

<a name="br15"></a> 

4\.8 Flow condition ﬁle ises.xxx

The ﬁle deﬁnes the ﬂow conditions and boundary conditions used by the solver program ISES .

It also conﬁgures the program into its analysis and design modes by specifying appropriate

global variables and constraints.

GVAR(1) GVAR(2) ... GVAR(N)

GCON(1) GCON(2) ... GCON(N)

MINLin P1PTin SINLin XINL [ V1ATin ] <-- optional

MOUTin P2PTin SOUTin XOUT [ V2ATin ] <-- optional

MFRin HWRATin [ XSHKin MSHKin ] <-- optional

REYNin NCRIT

TRANS1 TRANS2 (TRANS1 TRANS2 for blade 2) ...

ISMOM MCRIT MUCON

BVR1in BVR2in

<-- optional (see below)

MOVX MOVY SCAL ROTA (MOVX MOVY ... for blade 2)... <-- optional (see below)

KMOD GMODin

<-- optional (see below)

<-- optional (see below)

<-- optional (see below)

KMOD GMODin

KMOD GMODin

.

.

Global variables and constraints. Lines 1,2

The list of integers GVAR(1) ... GVAR(N) , given in any order, speciﬁes the global variables

to be treated as formal unknowns in the overall Newton equation system. If any variable is to

take on a new value, it must be speciﬁed in the list — otherwise it will retain its current value

in idat.xxx. The global variables selected also indicate to ISES which mode it should operate

in (direct, inverse, etc).

The list of possible global variables is,

1 SINL inlet flow slope (S1)

2 SLEX grid exit slope

5 SBLE LE stagnation point (for each non-sharp LE blade)

6 PREX grid exit static pressure

7 BVR1 streamtube thickness mode 1 DOF

8 BVR2 streamtube thickness mode 2 DOF

10 REYN stagnation Reynolds number

11 PDX0 zeroth mixed inverse prescribed Pi DOF

12 PDX1 first mixed inverse prescribed Pi DOF

13 PDD0 second mixed inverse prescribed Pi DOF

14 PDD1 third mixed inverse prescribed Pi DOF

14

![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAUADASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAf/9k=)

<a name="br16"></a> 

15 MINL inlet Mach number

16 MAS1 differential mass fraction DOF

20 GMODn modal geometry DOF flag n = 1, 2, ... NGMOD

31 MOVX x-movement DOF for each blade

32 MOVY y-movement DOF for each blade

33 SCAL scaling DOF for each blade

34 ROTA rotation DOF for each blade (in degrees, CCW)

40 GPARk geometry parameter DOF flag k = 1, 2, ... NGPAR

The list of integers GCON(1) ... GCON(N) , in any order, speciﬁes the global constraints to

be used to close the Newton equation system. In eﬀect, these constrain the speciﬁed global

variables.

The list of possible global constraints is,

1 Drive inlet slope S1 to SINLin

2 Drive outlet slope S2 to SOUTin

3 Set LE Kutta condition (for all non-sharp LE blades)

4 Set TE Kutta condition (for all blades)

5 Drive over/under splitter mass fraction ratio to MFRin

6 Drive inlet P0a to PSTr0 ( = 1/gamma )

7 Drive streamtube thickness mode 1 amplitude to BVR1IN

8 Drive streamtube thickness mode 2 amplitude to BVR2IN

9 Drive inlet v1/ao1 to V1ATin

10 Drive outlet v2/ao1 to V2ATin

11 Fix left endpoint of freewall segment

12 Fix right endpoint of freewall segment

13 Fix dP/ds2 at left endpoint of freewall segment

14 Fix dP/ds2 at right endpoint of freewall segment

15 Drive inlet Mach M1 to MINLin

16 Drive inlet pressure P1/Po1 to P1PTin

17 Drive outlet Mach M2 to MOUTin

18 Drive outlet pressure P2/Po1 to P2PTin

19 Drive inlet Reynolds number to REYNIN

20 Drive GMODn to GMODnin n = 1, 2, ... NGMOD

21 Set Xshock from XSHK to XSHKIN

31 Drive MOVX to MOVXin

32 Drive MOVY to MOVYin

33 Drive SCAL to SCALin

15



<a name="br17"></a> 

34 Drive ROTA to ROTAin

40 Drive GPARk to GPARkin k = 1, 2, ... NGPAR

41 Set Geometry-Parameter Constraint 1

42 Set Geometry-Parameter Constraint 2

.

.

Inlet, outlet conditions. Lines 3,4

¯

MINLin = inlet relative Mach number M<sub>1</sub>

P1PTin = inlet static/inlet-total pressure ratio p¯<sub>1</sub>/p<sub>o1</sub>

¯

SINLin = inlet relative ﬂow slope S = tan(β ) = v¯ /u¯

1

1

1

1

XINL = inlet-condition location m<sup>′</sup>

1

V1ATin = inlet relative tangential velocity ratio v¯<sub>1</sub>/a<sub>o1</sub>

¯

MOUTin = outlet relative Mach number M<sub>2</sub>

P2PTin = outlet static/inlet-total pressure ratio p¯<sub>2</sub>/p<sub>o1</sub>

¯

SOUTin = outlet relative ﬂow slope S = tan(β ) = v¯ /u¯

2

2

2

2

XOUT = outlet-condition location m<sup>′</sup>

2

V1ATin = outlet relative tangential velocity ratio v¯<sub>2</sub>/a<sub>o1</sub>

The V1ATin or V2ATin values can be omitted if corresponding constraints (9) or (10) are not

speciﬁed. This allows the use of MISES 2.4 ﬁle format which did not allow speciﬁcation of

v¯<sub>1</sub>/a<sub>o1</sub> or v¯ /a . It is useful to note that the inlet Mach number, ﬂow slope, and tangential

2

o1

velocity are related by

v¯<sub>1</sub>

M<sub>1</sub>

S<sub>1</sub>

q

q

\=

a<sub>o1</sub>

γ

1

1 + − M<sup>2</sup> <sub>1 + S</sub>2

2

1

1

so that setting any two ﬁxes the third.

The inlet and outlet conditions above are in terms of hypothetical mixed-out states at the

speciﬁed locations m<sup>′</sup> and m<sup>′</sup> (XINL, XOUT) shown in Figure 1. The mixed-out state ρ¯, u¯, v¯, p¯

1

2

is deﬁned to have the same total mass ﬂow, m<sup>′</sup>-momentum, θ-momentum, and total enthalpy

as the actual ﬂow at that same m<sup>′</sup> station. At the inlet station m<sub>1</sub><sup>′</sup> , the ﬂow is assumed to be

isentropic, so that the “mixed-out” state is obtained analytically from the known stagnation

conditions, mass ﬂux, and angular momentum. At the outlet station m<sup>′</sup><sub>2</sub>, the mixed-out state

is deﬁned implicitly by

R

R

ρ¯u¯ Pbr = ρu br dθ

ꢃ

\=

dm˙

− m˙ <sub>suct</sub>

ꢂ

R ꢂ

ꢃ

R

R

2

2

ρ¯u¯ + p¯ Pbr = ρu + p br dθ = u dm˙ + p br dθ − u<sub>e</sub>m˙ <sub>suct</sub> − ρ V u Θ b

e

e

e

R

R

R

ρ¯u¯v¯Pbr = ρuv br dθ

R

\=

v dm˙

− v<sub>e</sub>m˙ <sub>suct</sub> − ρ V v Θ b

− h<sub>oe</sub>m˙ <sub>suct</sub> − ρ V h δ b

e

e

e

¯

ρ¯u¯h Pbr = ρuh<sub>o</sub> br dθ

\=

h<sub>o</sub> dm˙

o

e

e

oe

h

16

![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAPgDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAf/9k=)![ref1]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAGAFcDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAUHCAr/xAAkEAAABAUFAQEBAAAAAAAAAAABAgMEAAUGBwkIFVeX2BETFP/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDv4ivLmluKSjanWtM1ot/cYsnWVpBhcKZTuSUU7naIl/lQqmbU1KJ/PWsnMUyorqy2TTFyQwJ/m1OAmEEIDBrKp8tLxRdBrbXHYDhum1Xc/tezUokiBXxFFEBbnRsAoodQ5Uji6KqQhSHAgJGUATCEhuuXjjfHJ3rqd89QhAN1y8cb45O9dTvnqG65eON8cneup3z1CEA3XLxxvjk711O+eobrl443xyd66nfPUIQDdcvHG+OTvXU756huuXjjfHJ3rqd89QhAaA0/vNabmcVAXVJS+mCQSAstajSy1hLhXZrObuJwLoQeJT9rcW2lCsmUtKy+GbLy52/dKOvqardNL4qKEID/2Q==)![ref5]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAGADkDASIAAhEBAxEB/8QAGQABAAIDAAAAAAAAAAAAAAAAAAIIBQcK/8QAKRAAAQMDAgUDBQAAAAAAAAAAAwECBAAFBwYJCBIVV5cRE9gUFlST0//EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDv4rX2Sm5BZo7U5cUxdHT8ipZpJNHwdf3K9WXRcu/Bb6wQapu2m7Tfr7FspHOckstts9xlDRGqKKVVVEUoKEw9Ubs8w540bGu3ckoIY8qR72bOJMUdBS1KgXAIHABCuM9Qk+pY8bBsRBe28nM7kyXVd3jtvtyedeJ349UpQOq7vHbfbk868Tvx6p1Xd47b7cnnXid+PVKUEhXndqZKhdSx3t3igPlgZMJBzZxKypjIyu9TujBl4Ajx3nQSPUSEMNik5Uc5rVVyXp+4ZH4MX9r/AOVKUH//2Q==)

<a name="br18"></a> 

where P is the angular pitch, u, v are the m<sup>′</sup>, θ velocity components, and V =

<sup>√</sup>u<sup>2</sup> + v<sup>2</sup> is

R

the speed. The integrals ( ) dm˙ on the righthand side are summed over all the inviscid

streamtubes at the m<sup>′</sup> station (in SUBROUTINE MIXOUT), which then requires including the

2

additional terms involving the momentum thickness Θ and total-enthalpy thickness δ<sub>h</sub>. These

inviscid streamtubes include the “removed” suction ﬂow (this is discussed later), which then

requires subtracting these ﬁctitious contributions explicitly via the m˙ <sub>suct</sub> terms.

Downstream of the blade row, the momentum and displacement thicknesses are extracted

directly from the solution at the sampled location. The total enthalpy thickness is not directly

available, since the thermal energy equation is not solved — the Reynolds analogy is used

locally at every surface location instead. Speciﬁcally, the surface heat ﬂux q<sub>w</sub> into the wall is

determined from the relation

q<sub>w</sub> = ρ u (h − h<sub>w</sub>) C<sub>h</sub>

e

e

aw

1

2

2

1

2

2

1

2

2

e

h<sub>aw</sub> = h + f u = I + (Ωr) + (f − 1) u

e

r

e

r

1

−2/3

C<sub>h</sub> = 0.22 Re− Pr

θ

where h<sub>w</sub> is the speciﬁed wall enthalpy, h<sub>aw</sub> is the usual adiabatic-wall enthalpy, and f<sub>r</sub> is the

temperature recovery factor. The correlation for the Stanton number C<sub>h</sub> is strictly correct only

for zero pressure gradients, although comparisons with ﬁnite-diﬀerence BL calculations indicate

that it is rearely more than 20% in error even for severe pressure gradient cases like high-work

turbines.

Since the evolution of the total enthalpy thickness δ<sub>h</sub> is not tracked, its ﬁnal value downstream

of the blade row must therefore be obtained from the thermal energy equation via the integrated

heat ﬂux over the blade surface.

Z

ρ V h δ<sub>h</sub>

≡

(h<sub>oe</sub> − h<sub>o</sub>) ρu r dθ

e

e

oe

d

(ρ V h δ b) = q br

e

e

oe

h

w

ds<sup>′</sup>

Z

˙

(ρ V h δ<sub>h</sub> b)<sub>exit</sub>

\=

q br ds<sup>′</sup> ≡ H

e

e

oe

w

w

Note that δ<sub>h</sub> is zero for an adiabatic-wall blade.

Using the state equation, the mixed-out total enthalpy is expressed as

γ

1

2

¯

¯ <sup>2</sup>

h

\=

p¯/ρ¯ +

V

o

γ−1

and the four relations above then form a closed system for the mixed-out state. Note that the

mixed-out state is assumed to have the same radius r and streamtube thickness b as the station

where the integration is performed.

Choosing the deﬁning station to be m<sup>′</sup> determines the mixed-out ﬂow quantities ρ¯ , u¯ , etc,

2

which are then used to impose global constraints on the outlet ﬂow angles, Mach numbers

2

2

and/or pressure, and also to compute the mixed-out losses at the exit:

p¯<sub>o2</sub> = p¯<sub>2</sub> <sup>ꢀ</sup>1 +

ꢁ

ꢄ

ꢅ

γ

<sup>v¯</sup>2

<sup>ρ¯</sup>2

γ−1

γ

−

1

¯

S =

2

¯ <sup>2</sup>

2

2

¯ <sup>2</sup>

M

2

M

\=

u¯ + v¯

2

2

2

u¯<sub>2</sub>

γp¯<sub>2</sub>

2

17

![ref9]![ref10]![ref10]![ref10]![ref3]![ref11]![ref4]![ref12]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADABkDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAf/9k=)![ref11]![ref3]

<a name="br19"></a> 

It is reassuring to note that for a viscous calculation, the mixed-out exit state and the loss in

particular (discussed later) are essentially independent of where it is calculated downstream

of the trailing edge. This veriﬁes that the quasi-3D coupled boundary layer formulation is

mass-,momentum-, and energy-conserving.

The inlet Mach number constraint (15) using MINLin and inlet pressure constraint (16) using

P1PTin are essentially equivalent, since

p¯<sub>1</sub>/p<sub>o1</sub>

\=

<sup>ꢀ</sup>1 +

ꢁ

−γ

γ

−1

γ−1

¯ <sup>2</sup>

M

1

2

Both are provided only for convenience.

Normally, the input parameters MINLin, P1PTin, MOUTin, P2PTin, are used only if the cor-

responding constraints (15),(16),(17),(18), are selected, with the inlet Mach number MINL (15)

being the appropriate global variable. If none of these constraints are selected, then MINLin

and SINLin will be used to simply set the inlet Mach and the mass ﬂow, which is then held

ﬁxed in the calculation. Because of this, it is normally not necessary to specify the inlet Mach

constraint (15) for subsonic ﬂows, since it will be satisﬁed anyway if constraints (15). . . (18) are

all omitted.

It is important to note that choked cases with speciﬁed subsonic inlet Mach number or inlet

pressure (i.e. a speciﬁed mass ﬂow) are ill-posed — physical considerations require that one

of the outlet constraints (17),(18) be used instead. Computationally, this does not suﬃce,

however. Normally, the ISES Newton procedure adjusts the inlet Mach variable MINL to meet

any speciﬁed conditions, but the inlet Mach and (and all other global variables) are treated

as being temporarily ﬁxed when the Newton matrix for the ﬂowﬁeld variables is set up for

each iteration. Each such iteration is therefore ill-posed, and ISES will complain either with

enormous residuals and/or an arithmetic fault due to a nearly-singular Newton matrix. A

similar problem occurs in axially-supersonic ﬂows, where the inlet Mach number cannot be

inﬂuenced by iterating on the downstream conditions.

To allow calculation of choked and/or axially-supersonic cases, it is necessary to select the

grid-exit pressure PREX (6) as a global variable. The corresponding global constraint is the

grid-inlet stagnation pressure constraint (6). This combination enables a computational “trick”

by which ISES can alter the inlet Mach number simultaneously with the ﬂowﬁeld variables, by

temporarily allowing variation of the inlet stagnation pressure. This stagnation pressure is then

driven back to its correct value by constraint (6). The overall procedure is merely the equivalent

of performing partial pivoting in the overall Newton matrix solution process. This produces a

diﬀerent iteration history (without arithmetic faults in particular!), but there is no eﬀect on the

ﬁnal ﬂow solution. The (6),(6) combination can of course be speciﬁed for unchoked subsonic

ﬂows, but it does produce a ∼ 10% CPU penalty and therefore should be omitted if no choking

is expected.

Over/Under splitter mass ﬂow ratio, wall temperature ratio. Line 5

MFRin = Speciﬁed mass ﬂow ratio above/below splitter blade

18

![ref11]![ref3]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAfQDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAf/9k=)

<a name="br20"></a> 

HWRATin = Speciﬁed h<sub>wall</sub>/h<sub>oa</sub>

= 0.0 for adiabatic cases

XSHKin = Speciﬁed shock x location for constraint (21) (optional)

MSHKin = Speciﬁed shock Mach number for constraint (21) (optional)

MFRin is used by constraint (5), typically in conjunction with variable MOVY (32), to position

the splitter blade 2 in the passage.

The wall temperature HWRATin is speciﬁed as a ratio relative to the rothalpy, h /h . If this

w

oa

ratio is speciﬁed as zero, the blade surface is taken to be adiabatic, as was assumed in MISES

2\.4 and all earlier versions.

XSHKin and MSHKin are optional, and needed only if constraint (21) is chosen. This can

be used to ﬁnd the exit pressure or inlet Mach which puts the shock at a speciﬁed location.

Typically, MSHKin = 1.0 would be speciﬁed.

Viscous ﬂow parameters. Lines 6,7

REYNin = Reynolds number

= 0.0 → inviscid calculation

(restarting a viscous case with REYNin = 0 “freezes” the boundary layers)

NCRIT = (+) critical ampliﬁcation factor “n<sub>crit</sub>” for e<sup>n</sup> transition model

= (−) freestream turbulence level (τ = -NCRIT, in %) for modiﬁed

Abu-Ghannam–Shaw bypass transition model

TRANS1 = side 1 surface transition trip m<sup>′</sup>/chord location

TRANS2 = side 2 surface transition trip m<sup>′</sup>/chord location

The input Reynolds number REYNin is based on the mixed-out static density, viscosity, and

relative speed at m<sub>1</sub><sup>′</sup> , and the reference length L<sub>ref</sub>.

¯

ρ¯ V L

1

1

ref

REYNin =

µ¯<sub>1</sub>

The reference length L<sub>ref</sub> is the same as was used used to deﬁne the streamsurface radii R in

the stream.xxx ﬁle described earlier. If a constant R=1 is speciﬁed (the default case for 2-D

cascades), then L<sub>ref</sub> becomes the length unit of the blade coordinates (X,Y). If (X,Y) are also

deﬁned so that the blade chord is unity, REYNin is then the usual chord-based Reynolds number.

The Reynolds number is used in conjunction with Sutherland’s formula to set the local viscos-

ity as a function of the local temperature. This requires another parameter, namely Sutherland’s

constant (T<sub>S</sub> = 110 K◦). This is stored internally in variable TSRAT as a temperature ratio:

T<sub>S</sub>

TSRAT =

<sup>T</sup>o1

where T<sub>o1</sub> is the relative-frame total temperature at m<sup>′</sup><sub>1</sub>. Currently, TSRAT is hard-wired to

0\.35 in src/iset/iset.f, although it can be read as an input parameter via the ises.xxx

19

![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAQIDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAf/Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADADgDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAf/9k=)![ref5]

<a name="br21"></a> 

input ﬁle if desired. The ratio of speciﬁc heats c /c = γ (GAM) can likewise be read in as

p

v

an input parameter. The read statement for both parameters is currently commented out in

SUBROUTINE ISESINIT (in src/ises/isesinit.f).

MISES incorporates a modiﬁed version of the Abu-Ghannam–Shaw (AGS) bypass transition

model (see separate document). Both the e<sup>n</sup> and the AGS model are active all the time and

either one may be decisive in inducing transition. Their respective input parameters n<sub>crit</sub> and

τ (Fortran names: ACRIT, FTURB) are always related through a modiﬁed Mack’s correlation,

and can be input either way. If a positive NCRIT is input, then this is taken as n<sub>crit</sub>,and the %

turbulence level τ for the AGS model is calculated from the modiﬁed Mack’s correlation:

<sup>n</sup>crit = NCRIT

ꢀ

τ<sup>′</sup> = 100 exp −

ꢀ

ꢁ

8\.43 + n<sub>crit</sub>

2\.4

ꢁ

2\.7 1 + τ<sup>′</sup>/2.7

τ

\=

ln

2

1 − τ<sup>′</sup>/2.7

If a negative NCRIT is input, then this is taken as the % turbulence level, and n<sub>crit</sub> is calculated

instead.

τ

= -NCRIT

τ<sup>′</sup> = 2.7 tanh(τ/2.7)

n<sub>crit</sub> = −8.43 − 2.4 ln <sup>ꢀ</sup>

ꢁ

τ<sup>′</sup>

100

The Mack modiﬁcation function τ<sup>′</sup>(τ) prevents negative n<sub>crit</sub> values for large τ values, and is

deemed reasonable given that Mack’s original correlation was developed for small τ levels.

The transition trip locations TRANS are deﬁned in terms of the fractional m<sup>′</sup> position on the

blade (which is not necessarily the fractional arc length position if r varies along the blade).

If an additional blade is present, simply add another TRANS3, TRANS4 pair on the same line.

n

Setting TRANSx ≥ 1.0 implies there is no transition strip on that side. Note that the e and

AGS criteria are always active, and free transition by either criterion may occur upstream of

TRANSx.

Isentropy and dissipation. Line 8.

ISMOM = 1 → use S-momentum equation

2 → use isentropic condition

3 → use S-momentum equation, replaced by an isentropic condition

only near the leading edge to minimize truncation errors there.

4 → use isentropic condition, replaced S-momentum equation only at shocks

where dissipation is active

MCRIT = critical Mach number in the deﬁnition of bulk viscosity

= 0.98 usually for weak shocks

= 0.85 for exceptionally strong shocks

20

![ref2]![ref3]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAEEDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAf/Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADABkDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAf/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAPoDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAf/9k=)

<a name="br22"></a> 

MUCON = artiﬁcial dissipation coeﬃcient (= 1.0 normally)

A negative MUCON value disables 2nd-order dissipation.

The value of MUCON may need to be increased to higher values (up to 1.5 or so) for strong shocks

and/or highly sheared grids. It may be appropriate to at the same time reduce MCRIT. Second

order dissipation is not recommended when when strong shocks traverse quasi-normal grid lines

— e.g. normal shocks on sheared grids, or strong oblique shocks on orthogonal grids. Note

that the sign convention on MUCON is reversed from previous MISES versions. This is to make

it consistent with the multielement MSES code convention from which the new dissipation

triggers were taken.

Streamtube thickness mode amplitudes. Line 9.

BVR1in = Speciﬁed streamtube thickness mode 1 amplitude BVR1

BVR2in = Speciﬁed streamtube thickness mode 2 amplitude BVR2

Geometry movement,scaling,rotation mode amplitudes. Line 10.

MOVXin = Speciﬁed x-displacement mode MOVX

MOVYin = Speciﬁed y-displacement mode MOVY

SCALin = Speciﬁed scaling mode SCAL

ROTAin = Speciﬁed rotation mode ROTA

Geometry shape mode amplitudes. Lines 11...end

KMOD = Speciﬁes geometric design mode

GMOD = Speciﬁed magnitude of geometric design mode MOD(KMOD)

4\.8.1 Variable,Constraint indices

The use of global variables and constraints gives the user a very ﬂexible means to apply boundary

conditions that specify the ﬂow condition or design conditions. In general, any number of global

variables can be speciﬁed, as long as the same number of global constraints are also speciﬁed

(except if a Modal- or Parametric-Inverse calculation is to be performed, as explained later). It

is only necessary that all the constraints properly constrain all the variables, and that the ﬂow

does not admit any non-physical situations. For example, the grid-exit slope variable DSLEX

(variable 2) can be constrained either by specifying it directly (constraint 2), or by specifying

the trailing edge Kutta condition (constraint 4) instead. Not specifying the Kutta condition is

not physical, however, and may produce strong pressure spikes at the trailing edge.

Note that for a multiple-blade case some of the variable and constraint options add a global

variable and a global constraint for each blade. For example, the speciﬁcation of trailing edge

21

![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAWEDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAf/Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAdsDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAf/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAXQDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAf/9k=)

<a name="br23"></a> 

Kutta condition (constraint 4) adds two constraints to the system — one for each blade trailing

edge. Speciﬁc examples of global variables and constraints are provided in a following section.

The GMODn (20) ﬂag indicates that some number of geometry design mode amplitudes are

selected, with their amplitudes speciﬁed in lines 11–EOF. The mode shapes associated with

each of these amplitude variables are deﬁned in the modes.xxx ﬁle, described later.

4\.8.2 Pressure-correction term

The pressure-correction term P<sub>corr</sub> (described in H. Youngren’s thesis) eﬀectively adds tension

to all the streamlines, and thus suppresses sawtooth modes in the grid. It is really only required

in inverse cases and in viscous cases with boundary layer separation, since the sawtooth modes

are adequately constrained by solid-wall and δ∗-oﬀset boundary conditions. It is also helpful in

cases with strong shocks traversing a sheared grid which has large streamtube widths compared

to cell lengths. Since this term is not “smoothing” in the normal sense, and is not dissipative,

it is simplest to use it whether it is necessary or not.

Since MISES v 1.4, the P<sub>corr</sub> term has been reformulated slightly from H. Youngren’s form.

The dependence on the local streamtube area has also been eliminated, making the term have

more or less equal inﬂuence throughout the ﬂowﬁeld. Previously, it was often too strong in the

thin streamtubes adjacent to the blade, and too weak in the larger interior streamtubes. The

result of the reformulation is that a larger P<sub>corr</sub> weighting factor PCWT can now be safely used

for diﬃcult cases, particularly those with strong shocks traversing sheared grids. For MISES v

2\.5, PCWT is hard-wired in SUBROUTINE ISESINIT for simplicity, since there is little reason

to treat it as an input parameter.

4\.8.3 Momentum/Entropy conservation

The ISMOM ﬂag controls whether and where S-momentum (streamwise momentum) or stream-

line total pressure (entropy) are conserved. The streamtube-cell equation residuals for each case

are:

m

ISMOM=1 :

ISMOM=2 :

R<sub>1</sub> ≡ ∆p + ∆q˜ − P<sub>s</sub> − p ∆P = 0

A

R<sub>2</sub> ≡ p ∆ (ln p˜<sub>oa</sub>) − p ∆P = 0

where m is the streamtube’s mass ﬂow, A is the streamtube’s cross-sectional area, and P is the

prescribed loss described earlier (R is slightly modiﬁed for the case of a curved streamtube).

1

The changes ∆( ) are along the streamwise direction in the cell. P is the streamwise centrifugal

s

force, plus the streamwise pressure force contribution to the cell from streamtubes above and

below, arranged so that the net equation is strongly conservative. Dissipation is introduced in

the form of an upwinded speed q˜ described in the next section. The upwinded “absolute” total

pressure (described later) is deﬁned in the standard manner, but using the upwinded speed.

p˜<sub>oa</sub> = p <sup>ꢀ</sup>

ꢁ

γ

I

γ

−1

I − q˜<sup>2</sup>/2 + Ω<sup>2</sup>r<sup>2</sup>/2

22

![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAA8DASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAf/Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAHYDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAAAAAAAf/9k=)![ref3]

<a name="br24"></a> 

ISMOM=1 gives a standard momentum-conserving Euler solver, while ISMOM=2 gives the

equivalent of a standard entropy-conserving Full-Potential solver. ISMOM=3,4 are hybrids

which attempt to make the use of the best features of the two types of solvers: correct Rankine-

Hugoniot shock jumps of the Euler solver, and zero total pressure loss of the Full-Potential

solver.

For ISMOM=3, the S-momentum equation R = 0 is used everywhere except in a region

1

near the leading edge, where the isentropic equation R = 0 is used instead. The size of the

2

isentropic region is hard-wired into a data statement at the top of SUBROUTINE SETUP

(in src/ises/setup.f), and typically extends from the inlet plane to ∼ 10 cells downstream

of the stagnation point, and 4 cells below and above the stagnation streamline. This should

cover most typical airfoils, but can be changed if appropriate for an unusual case. The only

exception is supersonic inlet casades, in which the bow shock will traverse the isentropic region.

IPLOT allows the display of the entropy-conserving region on the grid and ﬂow contour plots,

the remainder being the momentum-conserving region.

MISES v 2.0 is the ﬁrst to incorporate a new ISMOM=4 option, which makes all cells isen-

tropic except those where artiﬁcial dissipation contributes signiﬁcantly to the momentum ﬂux

(described later). In practice, this is nearly the same as ISMOM=3, but it will always give

momentum conservation at shocks, and hence is “safer” than ISMOM=3. It also doesn’t rely

on the ad-hoc hard-wired isentropic region, and hence is more automatic. Note also that IS-

MOM=4 is equivalent to ISMOM=2 for subcritical cases. The only possible problem with

ISMOM=4 is that because the switching between equations is ﬂow-dependent (as opposed to

being hard-wired), there is the possibility of the Newton iteration process stagnating at a limit

cycle. This has been observed occasionally in supersonic fans, so ISMOM=3 might be better

for these cases.

Because minimizing total pressure losses is much more important than conserving momentum

(except at shocks, of course), it is strongly recommended that ISMOM=3 or 4 be used for all

ﬂows. Even very small total pressure errors in the vicinity of the leading edge (where they are

most likely to occur) can cause dramatic errors in lift and drag of a viscous case at near-stall

conditions. This is because for a given imposed streamwise pressure gradient dp/ds, the edge

velocity gradient du<sub>e</sub>/ds which drives the viscous layers depends on the total pressure. Enforcing

isentropy (even if only at the leading edge for ISMOM=3) avoids most such problems, and hence

signiﬁcantly reduces grid density requirements. With ISMOM=3, it is only necessary to ensure

that no shocks traverse the isentropic region, otherwise incorrect Rankine-Hugoniot jumps will

result and the wave drag will not be properly predicted. ISMOM=4 does not have this possible

problem, and is recommended for general use.

4\.8.4 Artiﬁcial dissipation

The artiﬁcal dissipation in MISES is a speed-upwinding formulation analogous to bulk viscosity.

Instead of the actual speed q, the momentum and/or isentropy equations are deﬁned using an

23



<a name="br25"></a> 

upwinded speed q˜ deﬁned by

q˜ = q − µ

<sup>ꢀq</sup>i <sup>− q</sup>

i−1

ꢁ

i−1

<sup>+ µ</sup>i <sup>ꢀ</sup>q <sup>− q</sup>i−2<sup>ꢁ</sup>

(1)

i

(2)

i

i

where i is the grid node index along a streamtube, and µ<sup>(</sup><sub>i</sub><sup>1)</sup>, µ<sub>i</sub><sup>(2)</sup> are the ﬁrst- and second-

order dissipation coeﬃcients. To maintain numerical stability and allow shock capturing, the

following formulas for the dissipation coeﬃcients used, as suggested by a stability analysis.

!#

1 − 1/M<sup>2</sup>

(1)

i

<sup>C</sup>µ

1

2

2

2

2

i−1

µ

\=

(1 − M<sub>crit</sub>) log <sup>"</sup>1 + exp</sup>  

,

with M = (M + M

)

i

γ

1 − M

crit







(1)

i

µ

;

;

2nd-order dissipation

1st-order dissipation

(2)

i

µ

\=

0

In the limit M<sub>crit</sub> → 1, the above formula asymptotes to

µ → max <sup>ꢀ</sup> 0 , (1 − 1/M<sup>2</sup>)<sup>ꢁ</sup>

C<sub>µ</sub>

γ

which is the form indicated by a formal stability analysis, and is the dotted curve in Figure 2.

For M<sub>crit</sub> < 1, the full formula produces somewhat larger µ values, which asymptote rapidly to

zero below M = M<sub>crit</sub>. The intent of introducing M<sub>crit</sub> is to provide a user-adjustable margin

of safety. Figure 2 shows the variation of µ versus M<sub>crit</sub> and the overall scaling constant C<sub>µ</sub>.

Previous MISES versions used the simpler form µ ∼ max(0, 1 − M /M<sup>2</sup>), which has a slope

2

crit

discontinuity at M = M<sub>crit</sub>. The new form appears to have better shock propagating properties,

most likely due to its exponential “tail” for subsonic Mach numbers.

The second-order term in the formula for q˜<sub>i</sub> above is constructed to cancel the ﬁrst-order term

in the case of a linear variation of q(s). Analytically, the net result of either the ﬁrst-order or

second-order upwinding is to add a term to the streamwise momentum equation

∂p

∂s

∂q˜

\+ ρ q − ρΩ<sup>2</sup> r = 0 −→

∂s ∂s

∂r

∂p

∂s

∂q

\+ ρ q − ρΩ<sup>2</sup> r = −ρ q

∂s ∂s

∂r

∂(δq)

∂s

where δq = q˜− q is the upwinding modiﬁcation to the speed. This added term is dissipative,

producing streamwise changes in the absolute total pressure

1 ∂p<sub>oa</sub>

p<sub>oa</sub> ∂s

1 ∂(δq)

= −γM<sup>2</sup>

q ∂s

and is the mechanism which allows captured shocks to generate total pressure loss. In the

ISMOM=4 option, the magnitude of the fractional upwinded total pressure change over a cell

implied by δq is monitored:

∆ (ln p˜ ) ≃ δq ∆ <sup>ꢀ</sup>

ꢁ

ρq

oa

p

If the righthand side term is suﬃciently negative, below some small tolerance −ǫ (typically at

p

a shock), then R = 0 is used instead of R = 0, so that this total pressure loss is realized

1

2

in the solution. The local velocity gradient is also examined, with acceleration favoring R<sub>2</sub>

24

![ref3]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAD8DASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAf/9k=)![ref10]![ref3]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADABMDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAf/9k=)![ref13]![ref13]![ref14]![ref13]![ref13]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADACUDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAf/Z)![ref1]![ref11]![ref4]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADACUDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAf/Z)![ref12]

<a name="br26"></a> 

*M*

crit

*C*

**reduction of**

**increase of**

*M*

*M*

*M 1*

crit

*M*<sub>crit</sub>

Figure 2: Eﬀect of dissipation parameters M<sub>crit</sub> and C<sub>µ</sub> on dissipation level.

and deceleration favoring R . To minimize the tendency for non-convergent limit cycles, the

1

equation switch is implemented as a continuous blend. The actual general equation solved with

the ISMOM=4 option is

f R<sub>1</sub> + (1 − f) R<sub>2</sub> = 0

where 0 ≤ f ≤ 1 is a blending fraction which depends on the local magnitude of the total

pressure loss, as well as the local speed gradient.

It should be pointed out that this seemingly ad-hoc blending of R and R is perfectly

1

2

legitimate within the accuracy of the numerical scheme, since in smooth ﬂow these two equations

2

are equivalent to within O(∆s ). Hence, away from shocks, any change in f will produce a

2

solution change of at most O(∆s ) — the same as the truncation error of the discretization

scheme.

4\.8.5 Artiﬁcial dissipation level selection

The magnitude of the upwinding (e.g. magnitude of δq) is controlled by the approximate

threshold M<sub>crit</sub> (MCRIT), and the weighting factor C<sub>µ</sub> (MUCON). If MUCON is speciﬁed as

negative, then C<sub>µ</sub> = |MUCON|, and µ<sup>(2)</sup> = 0.

Lowering M<sub>crit</sub> and increasing C<sub>µ</sub> both increase the amount of upwinding, but in diﬀerent

ways as shown in Figure 2. The eﬀect of C<sub>µ</sub> on the numerical normal-shock proﬁle is shown in

Figure 3. In general, C ≃ 1 gives the cleanest normal shocks. For oblique shocks, the eﬀective

µ

2

coeﬃcient is C / cos θ, with θ being the shock angle. This favors somewhat smaller values of

µ

C<sub>µ</sub>. The stability analysis indicates that in general the minimum allowable coeﬃcients are

1st-order dissipation: C<sub>µmin</sub> = 1/2

2nd-order dissipation: C<sub>µmin</sub> = 1/4

with M<sub>crit</sub> ≤ 1 being required in all cases. Violation of these thresholds will produce numerical

instability and a nearly-singular Newton matrix.

25

![](data:image/jpeg;base64,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)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAANAAsDASIAAhEBAxEB/8QAFgABAQEAAAAAAAAAAAAAAAAACQgK/8QAHBAAAgMBAQEBAAAAAAAAAAAABQYDBAcIAgET/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/ANOPVfVDfn3Rb0WFOHmviHEfJep9B7iL+NkWfi7PQB2kLpc05hpbDf8AcHlhUNITb+nkRgEfAWoenEArkbktZhHABRJB+VW1z0DnDFtD0E5TYW7Rc+X9DKXh4Sqv0q0b5U8towDWG07d2D81YQYoLPwj+3mY19EfTVqtStEJqdeMlfinH9audJMPRIejuZHVOinf0c8sERsCOkzlWXnNByvNSgQAyVwTHQzBMe2gQDMGBk89q8QrNE9GJsBrhsOk6WrgEhPVkxUGQhVdRXgywth68k8sAoCBH1xYgdFNalntTeKQ+rXreZrM81iX5H+k8skvr17+h//Z)![](data:image/jpeg;base64,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)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAANAAsDASIAAhEBAxEB/8QAFgABAQEAAAAAAAAAAAAAAAAACQgK/8QAHBAAAgMBAQEBAAAAAAAAAAAABQYDBAcIAgET/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/ANOPVfVDfn3Rb0WFOHmviHEfJep9B7iL+NkWfi7PQB2kLpc05hpbDf8AcHlhUNITb+nkRgEfAWoenEArkbktZhHABRJB+VW1z0DnDFtD0E5TYW7Rc+X9DKXh4Sqv0q0b5U8towDWG07d2D81YQYoLPwj+3mY19EfTVqtStEJqdeMlfinH9audJMPRIejuZHVOinf0c8sERsCOkzlWXnNByvNSgQAyVwTHQzBMe2gQDMGBk89q8QrNE9GJsBrhsOk6WrgEhPVkxUGQhVdRXgywth68k8sAoCBH1xYgdFNalntTeKQ+rXreZrM81iX5H+k8skvr17+h//Z)

<a name="br27"></a> 

*C* < *C*

*C* < 1

*C* = 1

*C* > 1

min

(unstable)

Figure 3: Eﬀect of the dissipation weight C<sub>µ</sub> on the numerical structure of a captured shock.

Assumes M<sub>crit</sub> ≃ 1.

*M*

*M*

1

*M*<sub>crit</sub>

*M*<sub>crit</sub>

Figure 4: Eﬀect of the dissipation threshold M<sub>crit</sub> on the numerical structure of weak and strong

captured shocks. Assumes C<sub>µ</sub> ≃ 1.

Reduction of M<sub>crit</sub> mainly has an eﬀect where locally M ≃ 1. Speciﬁcally, the sharp subsonic

recovery downstream of the shock is smeared if the post-shock Mach number falls signiﬁcantly

above M<sub>crit</sub>, as shown in Figure 4. Hence, reduction of M<sub>crit</sub> tends to smear weak shocks, but

has little eﬀect on strong shocks – a rather undesirable situation. It is therefore desirable to set

M<sub>crit</sub> just below 1.0 to give a small margin of safety against numerical instability.

The second-order dissipation substantially reduces spurious total pressure errors wherever the

S-momentum equation is used, compared to the original ﬁrst-order dissipation. It also gives

crisper shocks for a given value of C<sub>µ</sub>, gives more reliable wave drag results, and is nice in

general. The beneﬁts are greatest for oblique shocks, which tend to be quite heavily smeared

with ﬁrst-order dissipation.

The major possible drawback of second-order dissipation is that it has been found to induce

more dispersion noise near a normal or near-normal shock which traverses a sheared grid. This

can be suppressed to some extent by increasing PCWT to 2.0 or so, but only up to a point.

One might also try a somewhat larger positive value of MUCON (= +1.5, say), which will still

give less smearing than a modest 1st-order MUCON (= -1.0, say). In any case, the option of

26

![](data:image/jpeg;base64,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)![](data:image/jpeg;base64,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)![](data:image/jpeg;base64,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)![](data:image/jpeg;base64,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)![](data:image/jpeg;base64,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)![](data:image/jpeg;base64,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)

<a name="br28"></a> 

disabling the second-order dissipation (which then reverts to the old ﬁrst-order scheme) has

been provided. The requirement for extra dissipation is obviously alleviated by using oﬀset

grids, which are much less sheared than periodic grids.

At the grid inﬂow plane, the upstream speeds q<sub>i−1</sub> and q<sub>i−2</sub> in the deﬁnition of q˜ above are

not available. Here, the upwinded speed is deﬁned by

(1)

i = 1 : q˜ = q<sub>i</sub> − µ<sub>i</sub> (q<sub>i</sub> − V<sub>inﬂow</sub>

)

i

where V<sub>inﬂow</sub> is the inﬂow-plane speed calculated directly from the inlet Mach number M<sub>1</sub>,

correcting for any radius and streamtube thickness changes between m where M is deﬁned

1

1

and the actual grid inﬂow plane location. This use of the inlet Mach number is the means by

which the extra incoming characteristic is set in supersonic-inlet ﬂows.

4\.8.6 Dissipation enhancement during convergence

The movement of a captured shock tends to be a slow process, with the movement distance per

Newton iteration being limited to the smeared shock’s thickness. An ideally-resolved shock can

therefore move at most one or two cells per iteration, which is deemed unacceptable. To speed

up this process, ISES takes the liberty of temporarily reducing M<sub>crit</sub> and/or increasing C<sub>µ</sub> to

smear the shock as much as possible so that it can move as fast as possible. The reduction of

M<sub>crit</sub> is especially eﬀective here, since this smears the subsonic side of the shock, allowing it

to move many cells per Newton iteration. The reduction is based on the maximum fractional

density change from the previous Newton step,

δ ρ¯ = | δρ/ρ |<sub>max</sub>

which typically is large at a moving shock.

The temporary value used for M<sub>crit</sub> is

ꢄ

ꢅ

(M<sub>crit</sub>)<sub>temp</sub> = 0.75 + (M<sub>crit</sub> − 0.75) exp −r<sup>2</sup>

δ ρ¯

where r =

or r =

ε

δ ρ¯<sup>3</sup>

ε (δ ρ¯<sup>2</sup> + (ε/4)<sup>2</sup>)

and ε = 0.15 is a reasonable threshold. When δ ρ¯ is large, the exponential is negligibly small,

giving the decreased value of (M<sub>crit</sub>)<sub>temp</sub>. As the shock reaches its destination, δ ρ¯ decreases, and

(M<sub>crit</sub>)<sub>temp</sub> reverts to its prescribed value, which then sharpens the shock. The two forms for

the argument r give nearly the same results, except very close to convergence where δ ρ¯ is small.

2

The simple form for r gives 1 − O(δ ρ¯ ) for the Gaussian, while the second more complicated

6

form gives 1 − O(δρ¯ ). Both forms give terminal quadratic convergence of the overall Newton

scheme, but the second form has a larger basin of attraction and is preferred for this reason.

The “old” fractional density change δ ρ¯ is saved in the idat.xxx ﬁle, so that the overall process

is automated, and is not outwardly visible to the user even if the Newton iteration is halted

and then restarted.

27

![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADABADASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAf/Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAGMDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAAAAAf/Z)

<a name="br29"></a> 

For future MISES releases, an automatic grid-sequencing procedure is planned as a replace-

ment for the current dissipation enhancement approach. This should also reduce overall com-

putation times for shocked ﬂow cases, by virtue of mostly eliminating the shock-propagation

bottleneck.

4\.9 Example ises.xxx input-ﬁle lines

Examples of selected lines of the ises.xxx input ﬁle are given below for a variety of situations.

A special feature of the input routine is that it looks for a “!” in the ﬁrst two lines containing

the global variable and constraint indices, causing the subsequent numbers to be ignored. For

example, the following two lines are equivalent:

1 2 5 ! 11 12

1 2 5

This is provided for convenience, since numerous indices must sometimes be repeatedly added

and deleted to reconﬁgure ISES for the various types of calculation cases.

4\.9.1 Lines 1–4. Variables, constraints, ﬂow conditions.

Quantities not used in the calculation are shown as zero, although MINLin should always be

set close to the anticipated inlet Mach, since it is sometimes used to initialize the ﬂowﬁeld for

iteration.

a) Speciﬁed subsonic inlet Mach and inlet slope (ﬁxed mass ﬂow), blunt leading edge.

1 2 5

|Global variables

1 4 3

|Global constraints

0\.80 0.

1\.50 -0.5

0\. 1.1

|Minl p1/po1 Sinl Xinl

|Mout p2/po1 Sout Xout

0\.

0\.

b) Speciﬁed outlet pressure and inlet slope (mass ﬂow unknown). Choking not permitted.

1 2 5 15

1 4 3 18

|Global variables

|Global constraints

0\.

0\.

0\.

1\.50 -0.5

1\.1

|Minl p1/po1 Sinl Xinl

|Mout p2/po1 Sout Xout

0\.50 0.

c) Same as b) above, choking permitted.

1 2 5 15 6

1 4 3 18 6

|Global variables

|Global constraints

0\.

0\.

0\.

1\.50 -0.5

1\.1

|Minl p1/po1 Sinl Xinl

|Mout p2/po1 Sout Xout

0\.50 0.

28



<a name="br30"></a> 

d) Supersonic/subsonic-axial inﬂow, speciﬁed inlet Mach, outlet pressure. The “unique-incidence”

condition will set the inlet slope.

1 2 5 15 6

15 4 3 18 6

1\.30 0.

|Global variables

|Global constraints

0\.

-0.5

1\.1

|Minl p1/po1 Sinl Xinl

|Mout p2/po1 Sout Xout

0\.

0\.50 0.

e) Same as d) above, but relative tangential inlet speed ratio v<sub>1</sub>/a<sub>o1</sub> is imposed rather than the

total inlet Mach number V /a .

1

1

1 2 15 6

9 4 18 6

1\.30 0.

|Global variables

|Global constraints

0\.

-0.5 0.88

1\.1

|Minl p1/po1 Sinl Xinl v1/ao1

|Mout p2/po1 Sout Xout

0\.

0\.50 0.

f) Same as d) above, but with sharp leading edge.

1 2 15 6

|Global variables

15 4 18 6

|Global constraints

1\.30 0.

0\.

-0.5

1\.1

|Minl p1/po1 Sinl Xinl

|Mout p2/po1 Sout Xout

0\. 0.50 0.

g) Mixed-Inverse design for case a) above.

1 2 5 11 12

1 4 3 11 12

|Global variables

|Global constraints

0\.80 0.

0\. 0.

1\.50 -0.5

0\.0 1.1

|Minl p1/po1 Sinl Xinl

|Mout p2/po1 Sout Xout

h) Modal camber design to attain speciﬁed outlet slope. Mode 9 is assumed to be a camber-

changing mode. Note that mode 9 is selected as a free variable, but the GMOD value is ignored,

since there is no (20) constraint speciﬁed. A negative mode index (e.g. -9) can be used to

selectively omit any mode from being imposed via constraint (20).

1 2 5 20

1 4 3 2

|Global variables

|Global constraints

0\.80 0.

1\.50 -0.5

0\.80 1.1

|Minl p1/po1 Sinl Xinl

|Mout p2/po1 Sout Xout

0\.

0\.

.

.

9

0\.0

|KMOD GMOD

i) Modal-Inverse design with 5 modes, driven by least-squares ﬁt to speciﬁed surface pressure

input via EDP .

29



<a name="br31"></a> 

1 2 5 20

1 4 3

|Global variables

|Global constraints

0\.80 0.

1\.50 -0.5

0\. 1.1

|Minl p1/po1 Sinl Xinl

|Mout p2/po1 Sout Xout

0\.

0\.

.

.

1

2

3

4

5

0\.0

0\.0

0\.0

0\.0

0\.0

|KMOD GMOD

|KMOD GMOD

|KMOD GMOD

|KMOD GMOD

|KMOD GMOD

j) Modal-optimization viscous design, 20 modes. Mode-sensitivity run.

1 2 5 20

1 4 3 20

0\.80 0.

|Global variables

|Global constraints

1\.50 -0.5

0\. 1.1

|Minl p1/po1 Sinl Xinl

|Mout p2/po1 Sout Xout

0\.

0\.

.

.

1

2

0\.0

0\.0

|KMOD GMOD

|KMOD GMOD

.

.

20 0.0

|KMOD GMOD

k) Parametric-Inverse, driven by least-squares ﬁt to speciﬁed surface pressure, with two user-

deﬁned constraints. Variable-index 40 takes parameter declarations from required ﬁle bplist.xxx

.

1 2 5 40

|Global variables

1 4 3 41 42

0\.80 0.

|Global constraints

1\.50 -0.5

0\. 1.1

|Minl p1/po1 Sinl Xinl

|Mout p2/po1 Sout Xout

0\.

0\.

l) Geometry-parameter sensitivity run.

1 2 5 40

1 4 3 40

|Global variables

|Global constraints

0\.80 0.

0\. 0.

1\.50 -0.5

0\. 1.1

|Minl p1/po1 Sinl Xinl

|Mout p2/po1 Sout Xout

4\.9.2 Lines 6–7. Viscous ﬂow parameters.

a) Inviscid analysis, or freeze current displacement thickness distributions.

30



<a name="br32"></a> 

0\.0E6 5.0

1\.0 1.0

|Re Ncrit

|Xtr1 Xtr2

b) Viscous analysis, free transition with Ncrit speciﬁed.

1\.3E6 5.0

1\.0 1.0

|Re Ncrit

|Xtr1 Xtr2

c) Viscous analysis, free transition with Turb speciﬁed in %.

1\.3E6 -1.0

1\.0 1.0

|Re -Turb

|Xtr1 Xtr2

d) Viscous analysis, BL trips at 50% chord on side 1, and 70% chord on side 2.

1\.3E6 5.0

0\.5 0.7

|Re Ncrit

|Xtr1 Xtr2

e) Two blades, viscous analysis, BL trips only on blade 1.

1\.3E6 5.0

|Re Ncrit

|Xtr1 Xtr2 Xtr3 Xtr4

0\.5 0.7 1.0 1.0

4\.9.3 Line 8. Isentropy and dissipation

a) Conserve S-momentum. Second-order dissipation.

1

0\.98 1.0

b) Conserve entropy (cannot be used with choked ﬂow).

0\.98 1.0 |Ismom Mcrit Mucon

c) Conserve entropy everywhere except at shocks (preferred).

0\.98 1.0 |Ismom Mcrit Mucon

d) Increased dissipation, shocks across sheared grid expected.

0\.90 1.5 |Ismom Mcrit Mucon

e) First-order dissipation, very strong shocks expected.

0\.95 -1.2 |Ismom Mcrit Mucon

|Ismom Mcrit Mucon

2

4

4

4

31



<a name="br33"></a> 

4\.9.4 Line 9. Streamtube thickness mode amplitudes

a) Streamtube thickness mode 1 amplitude speciﬁed directly.

1 2 5 7

|Global variables

1 4 3 7

|Global constraints

.

.

0\.1 0.

|B1 B2

b) Streamtube thickness mode 1 amplitude speciﬁed via outlet pressure, mode 2 amplitude

speciﬁed directly (note: mass ﬂow is ﬁxed here).

1 2 5 7

8

|Global variables

1 4 3 18 8

0\.30 0.

|Global constraints

1\.50 -0.5

0\.6 0. 1.1

|Minl p1/po1 Sinl Xinl

|Mout p2/po1 Sout Xout

0\.

.

.

0\. 0.1

|B1 B2

4\.10 Geometry perturbation mode speciﬁcation ﬁle modes.xxx

All the geometry deformation modes selected in ises.xxx must be described in the ﬁle modes.xxx.

This has the following format.

KMOD(1) IMODE(1) GWT(1) ST1(1) ST2(1) IEL(1)

KMOD(2) IMODE(2) GWT(2) ST1(2) ST2(2) IEL(2)

.

.

.

.

.

.

.

.

KMOD(N) IMODE(N) GWT(N) ST1(N) ST2(N) IEL(N)

where,

a) KMOD ties the mode shape to one of the modal geometry unknowns GMODn chosen via variable

(20). KMOD = 1 ties it to GMOD1, KMOD = 2 ties it to GMOD2, etc. A mode can be composed

of any number of disconnected pieces. Each piece of the mode is described on a separate line,

each beginning the same KMOD value.

b) IMODE speciﬁes the mode (or mode piece) shape. The shapes are implemented in FUNCTION

GFUN or variants thereof (in gmodes.f). The particular shapes currently implemented are:

32



<a name="br34"></a> 

IMODE n

n = 1 . . . 8

n = 20

mode shape

sin(nπs/s<sub>max</sub>

Local bump at leading edge, sized by local curvature

)

n = 21 . . . 40 T<sub>n−20</sub>(s/s<sub>max</sub>

)

n = 41 . . . ∞ sin<sup>a</sup>[π(s/s<sub>max</sub>) ]

b

sin<sup>a</sup>[π(1 − s/s<sub>max</sub>) ]

b

T<sub>n</sub> are Chebyshev polynomials, modiﬁed to be zero at the mode piece endpoints s/s<sub>max</sub> = 0, 1.

These are a good alternative to the simple sine waves (n = 1 . . . 19), since they give more

resolution at the endpoints. A good ﬁrst choice for mode shapes is to use a reasonably large

number of the Chebyshev modes (IMODE = 21, 22, . . . ) over the entire upper and/or lower

surface. This allows nearly arbitrary geometry variations. The mode shapes selected can be

plotted in EDP and in LINDOP .

c) GWT is the mode-piece scaling factor. If a mode consists of only one piece, GWT has no eﬀect,

as it merely rescales the mode amplitude variable. However, it is needed if the geometry mode

is composed of two or more pieces, and each piece must be scaled diﬀerently. An example is a

pure camber mode, where two identical shapes are placed on opposite sides of the airfoil, and

GWT is speciﬁed as +1.0 and -1.0 for the two pieces (a positive mode displacement is taken as

outward normal to the surface of the airfoil). GWT is also signiﬁcant in that it will alter the

convergence behavior (but not the ﬁnal answer) of the steepest-descent optimization process.

In this respect, it has been found advantageous to set the GWT factors for all the modes so that

the mode derivatives with respect to ST are all comparable in magnitude.

d) ST1, ST2 are the mode endpoint locations on the airfoil. These are the normalized arc lengths

s/s<sub>side</sub> from the airfoil nose (not the stagnation point like with Mixed-Inverse!), to the trailing

edge. For example, specifying ST1, ST2 = 0.0, 0.5 will result in the mode extending over

approximately the front half of the airfoil. There is also the option to specify ST1 and ST2 as

element x/c values. The necessary code is presently hibernating in SUBROUTINE GNSET (in

gnset.f) as comments, and only needs to be enabled.

e) IEL speciﬁes the target blade on which the mode piece acts.

4\.11 Design-parameter speciﬁcation ﬁle params.xxx

This ﬁle speciﬁes new design-parameter values, and is intended only for communicating with

the LINDOP optimization driver. The ﬁle contains values for...

SINLin

MINLin

P1PTin

V1ATin

33



<a name="br35"></a> 

SOUTin

MOUTin

P2PTin

V2ATin

REYNin

BVRNin n=1...NBVRN

GMODin n=1...NGMOD

GPARin n=1...NGPAR

MOVXin, MOVYin, SCALin, ROTAin

All these quantities overwrite the values in ises.xxx and bspec.xxx. A message is printed

when this occurs to warn the user.

Normally ﬁle params.xxx is written only by LINDOP . There is no reason to create it by

the user directly.

5 Program Descriptions

The descriptions for running the code on a UNIX system are given below. Similar, but diﬀerent

commands would be used for VAX/VMS systems. Starting from scratch, the usual program

execution sequence is

% iset xxx

% ises xxx

% iplot xxx

with the necessary input and output ﬁles for each step shown on the ISES Roadmap data ﬂow

diagram. The other programs are executed with the same command syntax, e.g.

% edp xxx

% iprint xxx

% bldset xxx

On UNIX, the execution of all these programs can be more conveniently and more rapidly

directed via the shell script run, For example, if one is computing case “xxx”, one would invoke

the shell script with

% run xxx

and any program can then be invoked for the xxx case with a quick menu selection.

34



<a name="br36"></a> 

5\.1 ISET

5\.1.1 Basic Initialization

ISET is the program which initializes the grid, the densities and a variety of other variables.

The required and optional input ﬁles to ISET are shown on the MISES Roadmap. The output

ﬁle is the main solution state ﬁle idat.xxx.

ISET is menu-driven to allow the user to iteratively generate a good initial grid by tweaking

a small number of gridding parameters. The top level ISET menu is

1

2

3

4

5

6

7

8

9

Specify new inlet slope and block off grid

Generate spacings and initialize grid

Elliptic grid smoother

Write idat.xxx

Plot grid

Plot Cp vs x/c

Modify grid parameters

Write grid parameters to gridpar.xxx

Change plot size

10 Read geometry from blade.xxx file

11 Read geometry from bparm.xxx file

Select grid generation option (0 to exit):

The generate the initial grid and write out the solution ﬁle, Options 1, 2, 3, 4 (in that order),

must be issued as a minimum. Normally, Option 1 is executed automatically when ISET is

started and can be skipped.

By default, ISET tries to read blade.xxx, but it can be forced to read and use the bparm.xxx

geometry parameter ﬁle from its menu. The default ﬁle is hardwired near the top of PROGRAM

ISET (in src/iset/iset.f), and can be easily changed if desired.

5\.1.2 Panel solution

Option 1 uses the speciﬁed inlet slope to generate an incompressible 2-D panel solution, which

is then used to trace a pair of stagnation streamlines just above and below each blade. Iso-

potentials emanating from all leading and trailing edges are also located. This divides up the

domain into blocks, which are then automatically displayed in a plot. These blocks form the

skeleton on which the grid is generated. The speciﬁed inlet slope here is of course somewhat

arbitrary, since it is only used for initial grid generation. It is a good idea, however, to specify

an inlet slope which avoids massive C<sub>p</sub> spikes on the leading edge. This minimizes the start-up

trauma with a subsequent viscous ISES solution. With small leading-edge radii, the range of

35



<a name="br37"></a> 

“tolerable” inlet slopes is quite small. Option 6 can be used to examine if the panel solution is

reasonable, and a new slope can be speciﬁed again with Option 1 if necessary.

To simplify ﬁnding a “reasonable” inlet slope for sharp or nearly-sharp leading edges, it is

also possible to set the inlet slope implicitly by choosing the leading-edge Kutta condition at

the Option 1 prompt.

Enter new inlet slope (or -999 to use LE Kutta):

Entering -999 will result in the inlet slope being set so that there is zero loading at the sharp

leading edge point, or at the nose tangency point (described in the EDP setion) for a blunt-

leading edge case. Enforcing the leading-edge Kutta condition will of course minimize any C<sub>p</sub>

spikes at the leading edge.

5\.1.3 Initial surface gridding

Option 2 distributes grid nodes along the streamlines on the blade surfaces. The local arc-length

increment ∆s between two surface grid nodes is determined from

1

∆s ∼

1 + a|κ|<sup>b</sup>

where κ is the local surface curvature. In regions of high curvature, the spacing is therefore

smaller, depending on the curvature exponent b and the coeﬃcient a. The exponent is speciﬁed

directly from the menu described below. A large exponent (b = 2, say), makes the spacing

small in high-curvature regions. A small exponent (b = 0.05, say), makes the spacing more

nearly uniform everywhere. The curvature coeﬃcient a is indirectly controlled by specifying ∆s

spacing at the leading edge (or stagnation point, to be more precise) of the blade. Note that

if κ is rather small at the stagnation point, the eﬀect of a is largely shut oﬀ in the expression

above. If necessary, a fudged additional curvature is added locally very near the stagnation

point to allow the spacing requirement to be met. A message is printed when this action is

taken. Fudged curvatures are also introduced near the trailing edge, and optionally at selected

local-reﬁnement zones on the upper and lower surfaces. The aim is also to control the local

spacing.

Since MISES v 2.1, the local/average spacing ratios ∆s/∆s<sub>avg</sub> are speciﬁed instead of the

actual spacings ∆s. The average spacing is deﬁned as ∆s<sub>avg</sub> = s<sub>side</sub>/N, where s<sub>side</sub> is 1/2 of

the airfoil perimeter and N is the number of airfoil-side points. Specifying ∆s/∆s<sub>avg</sub> is more

convenient than specifying ∆s itself, since the former automatically adjusts the spacing if the

number of points is changed.

The curvature exponent, and stagnation-point and trailing-edge spacing ratios and local-

reﬁnement parameters are altered from the following menu, which appears each time the spacing

is generated and plotted with Option 2:

36

![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADADkDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAf/9k=)

<a name="br38"></a> 

D sLE/dsAvg, dsTE/dsAvg spacing ratios

C urvature exponent

U pper side spacing refinement

L ower side spacing refinement

N umber of points

|

|

B lowup plot

R eset plot

Change what? (<Return> if spacing OK):

For each change request, the current values are displayed after the prompt. Selecting “U”, for

example, might produce

Enter new upper s/smax limits, local/avg density ratio 0.1000 0.2500 0.800

and just hitting <Return> will take the current values as the default input. One can change

only some of the required three inputs by using commas. Entering

0\.15

will only change the ﬁrst value from 0.1 to 0.15, while entering

,,0.5

will only change the third value from 0.8 to 0.5. The ﬁrst two “s/smax” values specify the

fractional arc length from stagnation point to the trailing edge where the local reﬁnement is to

be placed. Tick marks inside the airfoil element contour indicate this fractional arc length in

increments of 0.10. The local/avg density ratio speciﬁes the increase in local grid density over

the average density which would occur with all points spaced equally.

The point-number option “N” allows one to specify the average point number per side. With

oﬀset grids (described shortly), the upper and lower point numbers are diﬀerent.

After the spacing parameters are altered, the new distribution is generated and displayed.

The actual LE, TE, max, min, spacing ratios are also printed out. It must be mentioned that

only the stagnation point spacing ratio “dsLE/dsAvg” can be controlled precisely with the input

parameters. The other spacing ratios are approximate and may need to be iterated.

Once a good node distribution on each element is obtained, ISET proceeds to modify all the

spacings to resolve any conﬂicts between vacing blade surfaces. Basically, all the distributions

are automatically fudged to make spacings on element surfaces facing each other match as

closely as possible. This prevents massive grid shearing which would otherwise occur. No such

action is necessary for single-element airfoils.

Once all the node distributions are ﬁnalized, intermediate streamline nodes are generated in

the ﬂowﬁeld interior by simple linear interpolation from the stagnation and farﬁeld streamlines.

The resulting grid is not yet suitable for ISES calculations, but can be viewed with Option 5.

37



<a name="br39"></a> 

5\.1.4 Grid smoothing

Option 3 invokes an elliptic SLOR grid smoother to “clean up” the linearly interpolated grid.

This eliminates all kinks, overlaps, and also makes all the grid streamlines correspond to stream-

lines of incompressible inviscid ﬂow. This is then an excellent initial guess for the ISES solver.

5\.1.5 Initial solution ﬁle output

After the grid is smoothed, Option 4 can be issued to write out the initial solution ﬁle idat.xxx

which is then ready for the ISES solver. Before this is done, however, it is a good idea to view

the grid with Option 5. Options 2,3 or 1,2,3 can then be repeated if necessary to obtain an

acceptable grid before it is written out.

5\.1.6 Grid parameters

Option 7 puts up the menu

Current grid parameters:

G

I

O

S

X

A

t t inlet,outlet offset grid flags

30 number of inlet points

30 number of outlet points

15 number of streamlines

0\.500 x-spacing parameter

1\.500 aspect ratio of each cell at stagnation point

Change what (<return> if done)?:

in which each gridding parameter is displayed and can be immediately altered. The type of

grid topology is controlled by the inlet and outlet grid ﬂags. A ﬂag set to “t” speciﬁes an oﬀset

I-type grid, while “f” speciﬁes a periodic H-type grid. Diﬀerent grid types can be used over the

inlet and outlet. An oﬀset grid is very nearly orthogonal, but it increases the Newton matrix

bandwidth and thus signiﬁcantly increases the CPU requirements for a given grid size. On the

other hand, a nearly-orthognal grid requires much fewer streamlines, which mostly alleviates

this CPU penalty.

The great advantage of a nearly-orthognal grid is that it is much more tolerant of shock waves.

Above a certain amount grid shear, supersonic ﬂows are nearly incomputable. Conversely,

subsonic ﬂows cause little trouble even with strong grid shear, and are suitable in turbine

inlets, for example. The recommended grid topology for common cascade cases are as follows:

38



<a name="br40"></a> 

Flags

Case

f f low-speed, low-turning

f t high-speed turbine

t t transonic compressor tip

t f transonic compressor mid-section

The number of inlet and outlet points for periodic grids is set with options “I” and “O”. These

are ignored for oﬀset grids, since the inlet and/or outlet points cannot be set independently of

the surface points due to topological constraints. The number of streamlines is set with “S” for

all cases.

The “X-spacing” parameter controls the repelling-force between the quasi-normal grid lines

during the SLOR smoothing phase. There is little reason to adjust it from its default value.

The leading-edge cell aspect ratio controls the width of the streamtubes adjacent to surfaces.

Again, the default value suﬃces for most cases.

5\.1.7 Grid parameter saving, recall

Once a good set of gridding parameters is obtained, including the spacing parameters generated

with Option 2, they can all be saved to gridpar.xxx by specifying Option 8 at any time. If

this ﬁle already exists, it is overwritten. gridpar.xxx will then be automatically read when

ISET is executed again for that same xxx case, which causes all the gridding parameters to

take on their saved values. This allows rapid generation of grids for cases which diﬀer only

slightly (e.g. camber, inlet slope), since the same gridding parameters can then be used.

If the trailing edge of the blade is not closed, a constant-thickness wake gap is left extending

from the blade base. In inviscid calculations, this gap remains constant in width, but is free

to move up and down so that it sustains no pressure jump and hence no lift. For viscous

calculations, the gap will collapse down to the local wake displacement thickness, and is still free

to move up or down. A special treatment is used to correct for the dead air region immediately

behind the blunt base whose length and shape is set in ISET to match experimental correlations.

This special treatment results in an increase in momentum thickness downstream and accurately

accounts for “base drag”, which is also reﬂected in an increase in the mixed-out loss.

5\.1.8 Smoothing and writing the grid

ISET uses an elliptic grid generator to initialize the grid, which is invoked with Option 3 at

the top level menu. The grid can be plotted before or after smoothing with option 4 if desired.

If the overall grid is unsatisfactory, option 1 can be repeated as often as necessary.

Once a satisfactory smoothed grid is obtained, executing option 4 will ﬁrst initialize the

ﬂowﬁeld using hard-wired defaults, and then write out the unformatted state ﬁle idat.xxx.

The description of this ﬁle can be found from the comments in STATE.INC, which declares all

the variables in COMMON blocks.

39

![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCABnARgDASIAAhEBAxEB/8QAGQABAQEBAQEAAAAAAAAAAAAAAAkEBQgK/8QALxABAAAEAQwCAgEFAAAAAAAAAAECAwgEBQYHCRESMjlyeLTBEyIUIRUjUWFxgf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwD7+AAAAHLxcJfmjCpt+0JKlKFKX+tNPR27Yb23i+8Nz9fr7OoyYiWX7Txllmmk2TSxjDbsmlhHdj/zbHZ/sE4NBuk/SBlzWbXw6KcsZ35dyho/zB0AWd5y5qZj4zGTVc381cuZ8YvTxJnRljIlCMN2GMznlzdyPLlapLLJGMMj4LejPGEN2lSTFvNGlJrftYrTlpywlwtsNhX4+yH7pQxGNuVjVll/tLN8Un6/wrOAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAy4jgqdPqDUy4jgqdPqAJTW984HWO9sFgfm3MqxpOW984HWO9sFgfm3MqxgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAMuI4KnT6g1MuI4KnT6gCU1vfOB1jvbBYH5tzKsaTlvfOB1jvbBYH5tzKsYAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADLiOCp0+oNTLiOCp0+oAlNb3zgdY72wWB+bcyrGk5b3zgdY72wWB+bcyrGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAy4jgqdPqDUy4jgqdPqAJTW984HWO9sFgfm3MqxpOW984HWO9sFgfm3MqxgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAMuI4KnT6g1MuI4KnT6gCU1vfOB1jvbBYH5tzKsaTlvfOB1jvbBYH5tzKsYAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADiZVrTyy/FS+09SpSo/HGO5CrNWhNuU/l+3xQm3Zt6puT7uyH1jtAHiTQ9b3pDzQv8u6uTyxDIs2jbTZoXtczEzKrYfKs1fOCOW9EOI0x1c66WWMkRwlOTBYOnLn1kX+MxcuMrxx002MhGjh/x4fJ72AAAAAAAAAAAAAAAAAAAAAAAAAAAAAH/9k=)

<a name="br41"></a> 

5\.2 ISES

ISES is the main program that solves the Euler equations. It always reads the two required

input ﬁles, idat.xxx and ises.xxx, and then writes the output ﬁle back to idat.xxx. Thus

the input ﬁle idat.xxx can either be a restart ﬁle from an old calculation, or a new ﬁle created

by ISET . The optional input ﬁles will be accessed as needed, depending on the type of case

being run.

New for v 2.63:

Upon termination, ISES will also generate the small text ﬁle sens.xxx which contains the sen-

sitivities of various global variables, with respect to the speciﬁed parameters. These sensitivities

can be used for optimization procedures, solution interpolation and extrapolation, and a number

of other uses. The sensitivity calculation requires negligible CPU time, and is described in the

supplemental document MISES Constrained Least-Squares-Inverse Formulation and Sensitivity

Calculation Procedures.

ISES is run by typing run xxx and selecting the ises option, at which point the user is asked

for the number of Newton iterations “n” to be executed. The program response is as follows.

n = 0 writes solution to output ﬁle and terminates

n > 0 performs n iterations and repeats the question

n < 0 performs |n| iterations then writes solution to output ﬁle and terminates

After each iteration, the r.m.s. and maximum changes of the density, node position, and

viscous variables are displayed. Also displayed are the changes of various global variables.

Convergence to plotting accuracy occurs when the changes drop to about 0.1 × 10−<sup>3</sup> or so.

Convergence to machine accuracy is achieved when the changes refuse to go down further with

5

each iteration (about 0.1×10− ). ISES will terminate execution early if convergence is reached.

The convergence tolerances are speciﬁed in the include-ﬁle EPS.INC. If the Mach number is low

(below 0.1, say), or signiﬁcant ﬂow separation is present, the changes will not go down as far as

they would otherwise. This is due to the Newton matrix being less well-conditioned for these

cases. For incompressible viscous cases, it is recommended that a Mach number of at least 0.05

– 0.1 be used. This is eﬀectively incompressible.

ISES will automatically select the appropriate inﬂow/outﬂow boundary conditions based on

the local Mach number. Note that these BCs are distinct from the inlet/outlet mixed-out ﬂow

condition speciﬁcation at m<sup>′</sup> and m<sup>′</sup> . Rather, they determine the local nature of the solution

1

2

at the actual inﬂow and outﬂow grid planes. The details are given below.

5\.2.1 Inﬂow boundary conditions

The inﬂow BCs given in the table below use suitable combinations of the local ﬂow variables,

depending on the local streamwise and axial components of the mass-averaged inﬂow Mach

number.

40



<a name="br42"></a> 

Inﬂow type

Subsonic

Mach number

boundary condition







V<sub>1</sub> S<sub>1</sub>



q

M < 1 , M<sub>axial</sub> < 1 r(v + Ωr) = r<sub>1</sub>

\+ Ωr<sub>1</sub>

1 + S<sup>2</sup>

1

ꢄ

ꢅ

Supersonic/

Subsonic-axial

ˆ

ˆ

<sup>M > 1 , M</sup><sub>axial</sub> < 1 ν(M) ± β = ν(M ) ± arctan S

1

1





V<sub>1</sub> S<sub>1</sub>





q

Supersonic

M > 1 , M<sub>axial</sub> > 1 r(v + Ωr) = r<sub>1</sub>

\+ Ωr<sub>1</sub>

1 + S<sup>2</sup>

1

The subsonic and supersonic inﬂow BCs simply require that the angular momentum not

change between the inﬂow plane and the deﬁning station at m<sup>′</sup><sub>1</sub>. If both the radius and stream-

tube thickness are constant over the inlet region, then this is equivalent to imposing the inlet

streamline slope:

v

r(v + Ωr) = r (u S + Ωr ) −→

= S<sub>1</sub>

1

1

1

1

u

The supersonic/subsonic-axial BC relates the local ﬂow angle β = arctan(v/u) and the local

Mach number M via the Prandtl-Meyer function

s

arctan  <sup>sγ−1p</sup>M<sup>2</sup> − 1<sup>!</sup> − arctan <sup>ꢄp</sup>M<sup>2</sup> − 1<sup>ꢅ</sup>

γ+1

γ+1

ν(M) ≡

γ−1

ˆ

and thus allows waves to pass out of the grid inﬂow plane without reﬂection. The quantities M<sub>1</sub>

ˆ

and S are eﬀectively the same as M and S deﬁned at m<sup>′</sup> , but are corrected for any diﬀerence

1

1

1

1

between r,b at the inlet grid plane, and r ,b . The corrected ﬂow state ρˆ,uˆ,vˆ,pˆ, is related to

1

1

the m<sup>′</sup> ﬂow state ρ ,u ,v ,p , through conservation of mass ﬂow, absolute angular momentum,

1

1

1

1

1

constant rotation-corrected total pressure, and constant rotation-corrected total enthalpy (i.e.

rothalpy):

ρˆuˆ r b = ρ u r b

1 1 1 1

r (vˆ + Ωr) = r (v + Ωr )

1

1



!

−γ/γ−1

uˆ<sup>2</sup> + vˆ<sup>2</sup>

pˆ 1 −

= p<sub>oa</sub>

2I

γ

pˆ uˆ<sup>2</sup> + vˆ<sup>2</sup>

1

\+

= I + (Ωr)<sup>2</sup>

γ−1 ρˆ

2

2

The appropriate set of global variables and constraints for each of the three types of inlet

ﬂows is listed below (repeating some of the earlier input-ﬁle examples).

Inﬂow type Variables, constraints

1 15 6

1 18 6

Subsonic

Supersonic/

Subsonic-axial

1 15 6

15 18 6

1 15 6

9 18 6

or

1 15 6

15 18 1

Supersonic

41

![](data:image/jpeg;base64,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)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAAsDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAf/Z)![ref15]![ref11]![ref15]![ref11]![ref16]![ref16]![ref9]![ref11]![ref4]![ref9]![ref4]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCACAAUYDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAf/9k=)

<a name="br43"></a> 

In the supersonic/subsonic-axial case, the inlet slope constraint (1) is discarded in favor of

the inlet Mach constraint (15) or the inlet tangential speed constraint (9). Not specifying the

inlet slope explicitly is physically correct in light of the “unique-incidence” condition. In the

supersonic case, the inlet stagnation pressure constraint (6) is additionally discarded in favor

of the inlet slope constraint (1). Now the stagnation pressure will be implicitly constrained by

the fully-speciﬁed inlet ﬂow, which implies a speciﬁed mass ﬂow.

Using an inappropriate variable/constraint set for any inlet ﬂow type will usually still produce

a stable solution, but in general will result in some mismatch between the speciﬁed and resulting

quantity, e.g. M = MINLin , S = SINLin , p = 1/γ , etc.

1

1

oa

5\.2.2 Outﬂow boundary conditions

For the outﬂow, the boundary conditions again use some combination of local and global vari-

ables.

Outﬂow type

Subsonic

Mach number

boundary condition

<sup>h</sup>r(v<sup>isen</sup> + Ωr)<sup>i</sup> = 0 , S = S

∂

M < 1 , M<sub>axial</sub> < 1

<sup>M > 1 , M</sup><sub>axial</sub> < 1

exit

∂n

∂

h ꢄ

ꢅ

2

∂ p

i

Supersonic/

Subsonic-axial

isen

ν M

± β <sup>= 0 , S = S</sup>exit

∂n

Supersonic-axial M > 1 , M<sub>axial</sub> > 1

= 0

∂s<sup>2</sup>

where M<sup>isen</sup> is an isentropic Mach number deﬁned from the local static pressure p and radius

r, and the known inlet total pressure p<sub>o1</sub> and radius r<sub>1</sub>.







!



!

γ−1

γ

γ

ꢄ

ꢅ

isen

o

p

γ−1

1

2

2

2

2

 p

I + Ωr



<sub>M</sub>isen

isen

o

\=

− 1

p

= p<sub>o1</sub>

1

γ−1

I + Ωr<sup>2</sup>

2

1

The isentropic tangential velocity v<sup>isen</sup> is determined in a similar manner from the local static

pressure and the local ﬂow angle β. As with the inﬂow BCs, these outﬂow BCs are constructed

to be essentially invisible to the ﬂowﬁeld. In particular, the supersonic/subsonic-axial BC allows

waves to pass out without spurious reﬂection, although here the “transparency” isn’t quite as

perfect as at the inﬂow boundary if strong shock losses are present.

¯

The exit ﬂow slope S<sub>exit</sub> (which may be diﬀerent than the mixed-out ﬂow slope S ) is imposed

2

on one streamline at the outﬂow boundary, and then the ﬂow slopes of the remaining streamlines

are constrained by the ∂/∂n outﬂow BCs. S<sub>exit</sub> is typically chosen as a global variable DSLEX

(2), and it is normally implicitly constrained by specifying the trailing edge Kutta condition (4).

Setting the ﬂow slope explicitly via SOUTin using constraint (2) is not advised, since violating

the Kutta condition will result in a large ﬂow disturbance at the trailing edge. Constraint (2)

can be used in conjunction with a camber-changing mode variable (20), in a design case where

it is desired to modify the blade camber to attain a speciﬁed amount of turning.

42

![](data:image/jpeg;base64,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)![ref11]![ref17]![ref3]![ref10]![ref18]![ref10]![ref3]

<a name="br44"></a> 

The simple pressure gradient constancy condition for supersonic-axial ﬂow is in eﬀect a direct

extrapolation from the interior to the boundary, and does not inﬂuence the blade ﬂow solution

in any case. This outﬂow boundary condition doesn’t depend on S<sub>exit</sub>, and hence does not

require the outlet slope variable DSLEX (2) / TE Kutta condition (4) pair to be speciﬁed as a

global variable and constraint. These may be left in, however, with no ill eﬀects.

Besides the TE Kutta condition redundancy, supersonic-axial outﬂows are treated diﬀerently

in a number of other ways. Most notable is that the speciﬁed outlet pressure or Mach number

will generally not be what comes out of the ﬁnal solution. The local exit Mach number of an

axially-supersonic ﬂow is set entirely by the initial conditions at the blade row, and cannot be

prescribed. The ISES recognizes this and generally ignores any speciﬁed mixed-out constraint.

5\.3 IPLOT

IPLOT is the program which displays the solution in idat.xxx at any time whether the solution

is converged or not. It is executed by the command run xxx and selecting the iplot option.

Note that if the solution in idat.xxx is not converged, the results are physically meaningless.

The top-level IPLOT menu is shown below.

1

2

3

4

5

6

7

8

Blade surface plots

Streamtube plots

Contour/grid plots

Wake profile plots

r,b,ln(Po) vs m’ stream surface definition plots

Wheel View

Dump flowfield to text file

Dump BL quantities to text file

Select IPLOT option (0=Quit):

5\.3.1 Blade surface plots

The “Blade surface plots” menu brought up by the top-level option 1 allows plotting of most

of the airfoil surface and wake boundary layer variables:

1

3

Mach vs x

Hk vs x

2

Cp vs x

7

8

9

Ue vs x

A/Ao vs x

Ct vs x

4 s1 D,T vs x

5 s2 D,T vs x

6

Cf vs x

10 Rtheta vs x

43



<a name="br45"></a> 

11 Forces

12 Options

13 Change blade

14 Hardcopy toggle

15 Change x-axis coordinate type on BL plots

16 Change x-axis limits on BL plots

17 Change y-axis limits on current BL plot

18 Cursor blowup of current BL plot

19 Reset x,y-axis limits for BL plots

20 Annotation menu

21 Plot-page options

Select surface plot option for blade 1:

The Mach and menu item 1 is the isentropic Mach number calculated from the local pressure.

The menu items (3 . . . 10) display boundary layer quantities for at most one element at a

time. Items 3, 6 . . . 10 show one variable on both sides of the element, while items 4,5 show

δ∗ and θ together for one side only. Items 4,5 also show the total (top + bottom) thicknesses

for the wake as dotted lines, and also the inviscid-grid wall-oﬀset distance ∆n as a dashed line.

Normally, ∆n = δ∗, and the two curves will overlay, but only if the case is fully converged. If

the dashed ∆n curve can be discerned, the case is not converged. A number of plot coordinate

types can be selected with item 15. Items 16,17 allow rescaling of the BL variable plot axes to

zoom in on details of interest.

One important feature of IPLOT which needs some elaboration is the normalizing conven-

tions for C<sub>p</sub>, forces, losses, etc., which are listed with the “Forces” menu item 11. In general, all

ﬂow quantities are normalized with static isentropic reference values deﬁned from the reference

Mach number M<sub>ref</sub>, at the radius r<sub>ref</sub>. When IPLOT is started, M<sub>ref</sub> is initialized to the inlet

Mach number M<sub>1</sub>, and r<sub>ref</sub> is initialized to r(m<sup>′</sup><sub>1</sub>), so that all quantities are referenced to the

usual “1” quantities used for most of the normalization in ISES . However, M<sub>ref</sub> and/or r<sub>ref</sub>

can be changed to anything using Option 12 in the “Blade surface plots” menu.

Three types of C<sub>p</sub> can be displayed from the surface plots menu:

p − p

1

2

C<sub>p</sub>

\=

1

2

ρ V

1

1

<sup>p − p</sup>1

¯

C

\=

\=

p

p<sub>o1</sub> − p<sub>1</sub>

p<sub>o1</sub> − p

C<sub>po</sub>

p<sub>o1</sub> − p<sub>1</sub>

The default type is C<sub>p</sub>, but any type can be chosen using Option 12.

44

![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADACsDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAf/Z)![ref10]![ref6]![ref6]

<a name="br46"></a> 

Two types of loss coeﬃcients are deﬁned from the hypothetical mixed-out state ρ¯ , p¯ . . . at

2

2

m<sup>′</sup><sub>2</sub>, described earlier.

isen

p

p

− p¯

o2

o2

ω¯

\=

\=

p<sub>o1</sub> − p<sub>1</sub>

isen

− p¯

o2

o2

ζ

isen

p

− p¯<sub>2</sub>

o2

In addition, inviscid shock-loss and viscous-loss components are deﬁned at the grid outﬂow

plane as derived in H. Youngren’s thesis:

Z

ꢄ

ꢅ

1

dm˙

isen

o

ω<sub>i</sub>

\=

\=

p

− p<sub>o</sub>

p<sub>o1</sub> − p<sub>1</sub>

m˙

ꢀ

ρV <sup>2</sup>Θ b<sup>ꢁ</sup>

1

p<sub>o</sub> ρV

ω<sub>v</sub>

<sup>p</sup>o1 <sup>− p</sup>1 p m˙

exit

The derivation of these two components assumes a small wake velocity defect, or H −1 ≪ 1, and

hence is less rigorous than that of the mixed-out loss ω¯. Nevertheless, in most cases ω +ω ≃ ω¯

i

v

to within a few percent, which is reassuring! The cases for which the two approaches diﬀer

signiﬁcantly are supersonic-exit turbines, and radial-outﬂow blading in general. Here, there is

no rigorous way to deﬁne the viscous and inviscid loss components, since the momentum defect

doesn’t asymptote to a nearly-uniform value over the outlet region.

5\.3.2 Suction

The ﬂuid withdrawn through the suction slot has several parameters which are signiﬁcant for

the overall stage or machine performance. Besides the speciﬁed suction coeﬃcient C<sub>Q</sub> deﬁned

earlier, the total enthalpy h<sub>os</sub> and pressure p<sub>os</sub> are also important. The complication here is

the diﬃculty in deﬁning the loss of the blade row. Figure 5 shows the lateral displacement of

the ﬂow due to suction being applied at the blade surface. In the computation, the mass is not

actually removed from the streamtubes, but rather the “removed” ﬁctitious portion of the ﬂow

(shown shaded in the Figure) is made to overlap the real ﬂow downstream of the slot. The same

treatment is used in viscous cases, except that and additional shift of δ∗ is superimposed on

top of the suction-induced displacement. The net boundary condition applied to the streamline

adjacent to the surface at streamwise location s<sup>′</sup> is

Z

1

s<sup>′</sup>

∆n<sup>′</sup> r = δ∗ −

−ρ v br ds<sup>′</sup>

w

w

ρ V b

e

e

where ∆n<sup>′</sup> is the normal distance (in the m<sup>′</sup>–θ plane) from the wall to the streamline. The

integral is simply the total mass ﬂow sucked away upstream of the s<sup>′</sup> location, and reaches the

ultimate value m˙ <sub>suct</sub> downstream of the blade row.

5\.3.3 Streamtube plots

This menu allows the plotting of various ﬂowﬁeld quantities versus s<sup>′</sup> (arc length in the m<sup>′</sup>–θ

plane). It is mainly a diagnostic tool.

45

![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAFADASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAf/Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAFADASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAf/Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADADwDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAf/9k=)![ref1]![ref19]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADABADASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAf/Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADABYDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAf/9k=)![ref20]

<a name="br47"></a> 

**.**

**.**

*m − m*

suct

**.**

*m*

suct

**.**

*m*

*/*

*V b*

suct *e e*

**.**

"removed" flow

*m*

suct

**.**

*m*

**.**

**.**

*m − m*

suct

**.**

*\**

*m*

suct

**.**

*m*

*/*

*V b*

suct *e e*

"removed" flow

**.**

*m*

suct

**.**

*m*

Figure 5: Inviscid and viscous ﬂows with suction.

46

![](data:image/jpeg;base64,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)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCABhAAgDASIAAhEBAxEB/8QAGgABAAIDAQAAAAAAAAAAAAAAAAMGBAcJCv/EAC0QAAAFAwEHAwQDAAAAAAAAAAECAwQGAAUHCAkREhQ5ebgVISITFhkxMkFh/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/APYzkPKQz3VBFNJDc8viRjYdbapRyDEJP6FfCKYuzhAbAGOrnbPTHXPR6VmufIS9H1BAblHVrpZeBLnucRVrO4pAptcogiY6m78cmTeFQD7lUzq6lMVlFQh93xUKBhAhtwiX/aUEr7q8Q/ty5I8lsVUo+6vEP7cuSPJbFVKA+6vEP7cuSPJbFVKwnizgNrfHnaaYvUWmzwn7Vf6BfkzWV1GYtcptT8AnKd0+QTO7bIrmbnMiQxkwOmAqAoKvieGEwltSsxQ6OSKWyBhqn043TU5khed35aVu7ZN8f5ThuLYlGoMu4SbGi2PbHDJVdLcwiqAOUEVCMVyuA5UEzKur7q8Q/ty5I8lsVUoD7q8Q/ty5I8lsVUo+6vEP7cuSPJbFVKCK4Lpl2u0RMI/w2cmRxOUA9yb9S2KfY36DiAB4hKAiIFAR/qlSy+2vrHtSMQS+8255b4nM9HuU8TxG+g2XWtl0yc1y5DsmqxcizcigNLr9gROSSYyr4GrNa3Wp00TdneLINllB0mpSlB//2Q==)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCABhAAgDASIAAhEBAxEB/8QAGgABAAIDAQAAAAAAAAAAAAAAAAMEBQcJCv/EAC8QAAAFAwEHAwIHAAAAAAAAAAECAwQGAAUHCAkREhQ5ebgTFRkWISIjJTEyQWH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A9jOQ8pDPdUEU0kNzy+JGNh1tqlHIMQk/sV8Ipi7OEBsAY6uds9sdc9HpWa58hL0fcEBuUdWull4Eue5xFWs7ikCm1yiCJjqbvjkybwqAfcqmdXUpisoqEPu/CoUDCBDbhEv+0oJX3V4h/blyR5LYqpR91eIf25ckeS2KqUB91eIf25ckeS2KqVVfKn+W2LOjGTIky2e2Qrc8VT/PIi7f6j8Wq2xm5Eu4Gjt4iHNpIqDxqJEN6YHS4lSqCvj1mbFe0flWJLCBHtszfponOpTIkpkX65P7nObDmOG4+h9gbS1X0HLSBQSGym5R+HxjllULPb02aCK/A3ADKyT7q8Q/ty5I8lsVUoD7q8Q/ty5I8lsVUo+6vEP7cuSPJbFVKCK4Lpl2u0RMI/w2cmRxOUA+5N+pbFP2N+wcQAPEJQERAoCP9UqWX219Y9qRiCX3m3PLfE5no9ynieI30Gy61sumTmuXIdk1WLkWbkUBpdfoCJySTGVfA1ZrW61Omibs7xZBssoOk1KUoP/Z)

<a name="br48"></a> 

5\.3.4 Contour/grid plots

This menu allows plotting of the ﬂowﬁeld. The streamline grid, ﬂow variable contours, and

Mach waves can be displayed. There are also miscellaneous options for locating cell indices,

shading isentropic cells, plotting BL proﬁles, etc.

5\.3.5 Wake proﬁle plots

The wake proﬁle plots allow the display of various inviscid ﬂow quantities versus tangential

distance, much like a boundary layer proﬁle. This is mainly a diagnostic tool, useful for checking

whether the shock defect wake is adequately resolved, for instance.

5\.3.6 r,b vs m<sup>′</sup> stream surface deﬁnition plots

This gives a plot of r(m<sup>′</sup>) and b(m<sup>′</sup>) as deﬁned in the stream.xxx ﬁle, showing the individual

X(I), R(I) spline node values as well as the resulting analytic spline function. The individual

streamtube thickness modiﬁcation modes B b (m<sup>′</sup>) and B b (m<sup>′</sup>) are also displayed.

1

1

2

2

5\.3.7 Wheel view

This displays a nifty picture of the entire rotor from the side and along the axis of rotation.

Only geometries with signiﬁcant radial changes will make a meaningful picture.

5\.4 EDP

EDP is an interactive menu-driven program used to modify data in idat.xxx, primarily for

inverse calculations. Three types of inverse methods are provided: Mixed-Inverse, Modal-

Inverse, and Parametric-Inverse. The latter two are essentially the same, and diﬀer only in the

manner in which the blade shape modiﬁcation is deﬁned. The table below lists the key features

of each method.

Method

Geometry Description

Solution Method

Mixed-Inverse

Pointwise (arbitrary shape)

require p(s) = p + . . .

spec

R

Modal-Inverse Surface-normal displacement modes minimize (p − p<sub>spec</sub>)<sup>2</sup> ds

R

Parametric Inverse User-deﬁned geometric parameters minimize (p − p<sub>spec</sub>)<sup>2</sup> ds

As the input-ﬁle examples indicate, Mixed-Inverse is triggered by global variables (11),(12),

while Modal- and Parametric-Inverse are triggered anytime the number of global variables

exceeds the number of global constraints. With the least-squares minimization, all the global

variables are determined so as to obtain the best-possible ﬁt of the computed surface pressure

p(s) to the speciﬁed surface pressure p<sub>spec</sub> input in EDP , while simultaneously satisfying the

global constraints that are declared. These declared constraints would be the usual ones used

47

![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCABYAjUDASIAAhEBAxEB/8QAGgABAQEBAQEBAAAAAAAAAAAAAAkBBggHCv/EAC8QAQAABAQEBgEDBQAAAAAAAAABAgMGBAUHCRExcbUINjl2eLEhEhNBIyYyUWH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A/fwAAAAAAAAAAAAAAyb/ABm5co8+XL+f+NZNyj0j9A+d3tisTgLNu3HYGpWo5hg7YzyrLXw00KUcNXwGVYupSr0Zofq/rUqsstTDw4/mMsYwjCPF4+2o9SL41f24PBjqfqTdudX7fd86B2Hcd03pceImxWfXPnGY5XJWxecZvWnjNNNjsZUjGpXhNNPGE0fzNHi9i6jUKM1h3rCNOWH9n3LPCaEOE0s9HJcb+3NLH+IywnmhCPDlNF4A2VeEdp7wA1YSyyz1/DJprXqxlhwhPVqZNTjPPGH+5o8wVAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAZNyj0j9NZNyj0j9A4TUTyLevs26uy4pPvZU9Jrb7+L+mPZpFBNRPIt6+zbq7Lik+9lT0mtvv4v6Y9mkBUAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABk3KPSP01k3KPSP0DhNRPIt6+zbq7Lik+9lT0mtvv4v6Y9mkUE1E8i3r7NursuKT72VPSa2+/i/pj2aQFQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTco9I/TWTco9I/QOE1E8i3r7NursuKT72VPSa2+/i/pj2aRQTUTyLevs26uy4pPvZU9Jrb7+L+mPZpAVAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAZNyj0j9NZNyj0j9A4TUTyLevs26uy4pPvZU9Jrb7+L+mPZpFBNRPIt6+zbq7Lik+9lT0mtvv4v6Y9mkBUAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABk3KPSP0AOJvWhiM0s+6cFl+Hq4rE4+17jwWDp0oSxjXxWKyrEUcNRk4zQh+qvVnlkp8eEIxj+YweMNqHTO/dGdtvwV6U6pWtmlk6iWB4fbAti8rSzqnSpZtb+f5ZlUlHH5Xj6dGrWpSYnC1YRkqy06tSWE0PxNEAUGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB//Z)

<a name="br49"></a> 

*s*

*B*

*A*

*B’*

*B*

**SINL**

Figure 6: Deﬁnition of surface arc length for Mixed-Inverse (from the stagnation point A), and

for Modal-Inverse (from the nose tangency point B). Tangency point moves during modal shape

changes, but the original point B remains ﬁxed.

for analysis calculations — speciﬁed inlet slope (1), TE Kutta condition (4), etc. Any global

variable can therefore play the role of a design parameter if it doesn’t have a declared constraint

associated with it, although this can easily give ill-conditioned problems if such a “free” variable

has only a weak eﬀect on the surface pressure. The intent of the generalized least-squares

procedure is to drive geometry-related variables for inverse design where the blade shape is

restricted to a particular parametric description. In contrast, the Mixed-Inverse method is

aimed at eﬃciently eliminating local aerodynamic defects (C<sub>p</sub> spikes, etc.) in design problems

which permit arbitrary blade shapes.

5\.4.1 Surface parameterization

The Mixed-Inverse method and Modal-Inverse method use a slightly diﬀerent way to deﬁne the

fractional surface arc length parameter σ = s/s<sub>side</sub>, illustrated in Figure 6. Here, s is understood

to be the arc length in the m<sup>′</sup> − θ plane.

Since the surface grid points move along with the stagnation point A, Mixed-Inverse essentially

deﬁnes the surface pressure and geometry in the streamwise grid node index: p(i), m<sup>′</sup>(i), θ(i)

. In contrast, Modal-Inverse deﬁnes these quantities in terms of the arc length: p(s), m<sup>′</sup>(s),

θ(s), with s measured from the ﬁxed nose tangency point B. This nose point is set once and for

all in ISET at the location where the surface-tangent is perpendicular to the initial grid slope

SINL. If the blade is deformed by a camber mode which extends all the way to B, then the nose

tangency point will move from point B to B’ as shown in Figure 6. However, the initial point

B is always retained for deﬁning s, and subsequent points B’ are ignored.

Of particular importance are the endpoints of the target-segment, which is the part of the

surface which is to be modiﬁed. In Mixed-Inverse, these endpoints are actually speciﬁed as grid

point indices i ,i , while for Modal-Inverse, they are speciﬁed as the normalized arc lengths

0

1

48

![](data:image/jpeg;base64,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)

<a name="br50"></a> 

σ = s /s , σ = s /s<sub>side</sub>. The endpoints for either Inverse method are speciﬁed along with

0

0

side

1

1

the surface pressure distributions in EDP .

5\.4.2 EDP execution

EDP is executed by run xxx and selecting the edp option. EDP ’s top level menu is:

1

2

3

4

EDIT pressures

SET redesign flags

WRITE idat.xxx

10 Read flow(spec) data file

11 Write flow(spec) data file

12 Write x,y,z blade file

READ idat.xxx

13 Write m’,theta blade file

5

6

7

8

9

Change idat.xxx params 14 Annotate plot

Change r,b vs x

15 Change plot size

Change DPo vs x

16 Hardcopy current plot

Print flow parameters

Plot mode shapes

Select EDP option (0=QUIT):

Most of these options are self-explanatory. Options 6 and 7 permit changing the r(m<sup>′</sup>), b(m<sup>′</sup>),

and/or P(m<sup>′</sup>) distributions currently stored in idat.xxx . Likewise, option 5 permits changing

of several variables and ﬂags which are otherwise inaccessible. The option 5 sub-menu can be

easily customized to allow changing of any quantity present in the STATE.INC global include

ﬁle.

The primary purpose of EDP is of course the interactive input of speciﬁed surface pressures

for inverse calculations. MISES is really a blade redesign system rather than a pure design

code, and productive use of EDP necessitates that idat.xxx contain a previously converged

case. The surface pressure modiﬁcation is then done via option 1, which puts the user into the

pressure-editing sub-menu:

\------------------------------------------------------

I nitialize Mach (spec)=Mach (wall) on target element

M odify Mach (spec)

D emark inverse segment

S lope-matching at segment endpoints (toggle) ->

F

F low data select (Mach, Cp, P/P0a, P/P01)

T oggle plot type (flow vs s, or vector plot)

Z oom

R eset zoom to original size

49

![ref21]

<a name="br51"></a> 

L imits, set plot limits

A nnotation menu

H ardcopy current plot

Select edit option:

The speciﬁed pressure array CPspec read in from idat.xxx is originally zeroed out in ISET ,

so if this is a ﬁrst inverse editing session, the user must select option I to initialize CPspec to

the current wall pressure coeﬃcient array CPwall. Option M will then allow the user to “edit”

CPspec with the screen cursor. This can be done repeatedly if needed. If necessary, the inverse

target segment endpoints can be cursor-speciﬁed with option D. This is done either on the C<sub>p</sub>(s)

plot or on the geometry/vector plot, depending on what’s currently on the plot screen (toggle

with option T). The initial default target segment endpoints for Mixed-Inverse are the ﬁrst grid

point after the front stagnation point, and at the rear trailing edge point. For Modal-Inverse

they are the nose and trailing edge locations.

i<sub>0</sub> = i<sub>LE</sub> + 1

i<sub>1</sub> = i<sub>TE</sub>

Mixed-Inverse default

Modal-Inverse default

σ<sub>0</sub> = 0.0

σ<sub>1</sub> = 1.0

After the pressure-editing menu is exited, the following redesign-ﬂag menu comes up (it can

also be brought up with option 2). The intent is to verify the endpoint locations and several

other control ﬂags, and change them from the keyboard if necessary. If a Parametric-Inverse

case is to be run, this menu can be ignored.

Current redesign flags...

B lade

:

1

2

2

M oved-side:

P spec-side:

I endpoints:

20

78 (for Mixed-Inverse)

S endpoints: 0.0000 1.0000 (for Modal-Inverse)

Select flag to change (<return> if OK):

The target-blade NDES, moved-side KSMOVE, and Pspec-side KSPRES ﬂags determine what blade

side(s) get redesigned and where the surface pressures are imposed. Figure 7 shows the four

possible ways to change the target blade, corresponding to the possible moved-side KSMOVE ﬂag

settings.

-1 Both sides move opposite (camber preserved)

0

1

2

Both sides move together (thickness preserved)

Upper side moves

Lower side moves

With Mixed-Inverse, KSMOVE = -1 or 0 forces corresponding grid nodes on opposite sides of the

blade to move opposite or together. For oﬀset grids, it is not possible to identify “corresponding”

50



<a name="br52"></a> 

**KSMOVE =**

**-1**

**0**

**1**

**2**

<i>s<sub>1</sub></i>

*s*

*0*

**SINL**

Figure 7: Two-side and single-side geometry changes.

point pairs, and hence Mixed-Inverse with KSMOVE = -1, 0 can only be used on non-oﬀset

(periodic) grids. In any case, this option is not recommended on periodic grids either. Camber

and thickness changes are best performed using Modal-Inverse.

The geometry changes can be driven by specifying the pressure jump (i.e. the loading) across

the blade, or by specifying the pressures on an individual surface. This is controlled by the

Pspec-side KSPRES ﬂag, which can take on the following values.

0

1

2

delta(p) across blade is specified

Upper side pressures specified

Lower side pressures specified

The following combinations produce ill-posed inverse problems and must be avoided.

KSMOVE KSPRES

1

2

2

1

0

-1

The last combination is ill-posed since thickness changes have little eﬀect on loading to ﬁrst

order.

If a screen cursor is not available, then top-level option 11 can be used to dump the surface

pressures (or equivalent Mach numbers) into a formatted scratch ﬁle, which can then be edited

manually, or preferably with the user’s own software. This ﬁle can then be read back into

EDP with top-level option 12 . The scratch ﬁle also contains the indices of the requested

inverse segment endpoints and the index of the blade side which contains the inverse segment.

However the CPspec distribution is generated in EDP , it must be written out with the

idat.xxx ﬁle so that ISES can access the information. This is done with top-level option 3.

ISES is conﬁgured for a Mixed-Inverse case by selection of the global variables (11) and

(12) and the global constraints (11) and (12). If desired, variables (13) and/or (14) can also

51

![](data:image/jpeg;base64,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)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAH4DASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAAAAAAAf/9k=)

<a name="br53"></a> 

be speciﬁed, with corresponding constraints (13),(14). These produce a smoother pressure

distribution at the inverse segment endpoint(s), although they are rarely necessary in practice.

5\.4.3 Modal-Inverse

The Modal-Inverse method restricts changes in the blade airfoil to be a sum of some speciﬁed

set of geometry modes, which deform the surface normal to its current shape. If the modes are

smooth, then the modiﬁed airfoil is guaranteed to be smooth as well. Modal-Inverse is quite

robust, is particularly useful in ﬂows with shocks and/or separation, where the Mixed-Inverse

method might produce and irregular airfoil shape or simply fail due to ill-posedness. On the

other hand, Modal-Inverse will match CPwall to CPspec in only a least-squares sense, and is

intended for changing the overall airfoil shape, rather than removing small geometric defects.

If a Modal-Inverse or modal optimization case is to be run, then the global variables (20)

must be chosen. The corresponding geometry modes and their endpoints must be speciﬁed in

ﬁle modes.xxx (described in the Optimization section below). The geometry mode shapes are

set up in FUNCTION GMODES (in src/ises/gmodes.f), and can be altered as desired. Each

mode f(σˆ) is deﬁned over the interval 0..σˆ..1 , with σˆ being the fractional arc length over the

target segment.

σ − σ<sub>0</sub>

σˆ =

σ − σ

1

0

The current modes implemented in SUBROUTINE GMODES consist of modiﬁed Tchebyshev

polynomials, plus one “tail-wagging” mode for changing the overall blade camber. It is essential

that the target segment extend to the trailing edge (i.e. σ<sub>1</sub> = 1) if the “tail-wagging” mode

is used. The modes implemented in SUBROUTINE GMODES can be plotted in EDP with

option 10.

It must be mentioned that the least-squares pressure-matching condition cannot be fully

linearized, with the result that the convergence rate will be somewhat more sluggish than

usual. Nevertheless, only a few more iterations relative to a usual analysis or Mixed-Inverse

case will be required. On the other hand, a viscous Modal-Inverse case can typically handle

limited separated regions within the target segment. It is not necessary to temporarily “freeze”

the boundary layers as with Mixed-Inverse.

In practice, it is rarely necessary to converge a Mixed-Inverse or a Modal-Inverse case down to

the usual analysis convergence tolerance. If the iterations are halted before full convergence is

reached, then the new blade shape is still quite usable — it might have been changed 95% of the

way towards the “correct” geometry rather than the 99.999% obtainable with full convergence.

This slight diﬀerence is irrelevant in actual iterative design applications, since the “correct”

geometry rarely turns out to be exactly what is required.

52

![ref16]

<a name="br54"></a> 

5\.4.4 Parametric-Inverse

The Parametric-Inverse method is conceptually the same as Modal-Inverse, and diﬀers only

in the manner in which the geometry changes are handled. With each update of the user-

deﬁned geometry parameters, the entire blade shape is nonlinearly recomputed in ISES by

calling SUBROUTINE BLDGEN. This allows MISES to be used as an inverse-design engine

with any geometry-deﬁnition system. This ﬂexibility of the inverse-design system can be greatly

enhanced by user-supplied geometric constraints coded in SUBROUTINE BPCON, as described

earlier.

5\.4.5 Blade Translation, Scaling, Rotation

The blade rigid-body movement and scaling variables (31-34) are intended to complement the

Modal-Inverse formulation, which can only change the shape of the blade. Scaling variable SCAL

(33), for example, can be used to grow or shrink the blade. Since the pitch remains ﬁxed, this

in eﬀect modiﬁes the cascade’s solidity. For single-blade cases, the blade m<sup>′</sup>-translation variable

MOVX (31) is meaningful only if the radius and streamtube thickness distributions r(m<sup>′</sup>), b(m<sup>′</sup>)

are nonuniform, and the θ-translation variable MOVY (32) has no eﬀect. For multi-blade cases,

these variables can be used to move the blades relative to one another.

The rigid-body variables (31-34) can also be used to augment the Parametric-Inverse for-

mulation if the user’s blade parameterization does not permit overall blade motion or scaling.

Conversely, if the user’s parameterization does contain any such rigid-body mode, then the

corresponding global variable (31-34) must not be used, since the system will have redundant

variables and be very ill-conditioned.

The rigid-body variables (31-34) can be driven either directly, via their speciﬁed values in the

ises.xxx ﬁle, or indirectly, via some other constraint. For example, the blade rotation variable

ROTA (34) can be driven by specifying an outlet slope (2) or outlet tangential velocity (10).

The rigid-body variables can also be used for Modal-Inverse or Parametric-Inverse calculations

in conjunction with the wall pressures speciﬁed in EDP . Adding the rotation variable (34) to

the usual mode variables (20) can often produce a better ﬁt to these speciﬁed wall pressures,

since changing the blade shape alone may not be adequate.

5\.4.6 Modiﬁed-Blade Output

Option 5 writes out the current m<sup>′</sup>, θ coordinates in a new blade.yyy ﬁle, with the user being

prompted for the new ﬁlename. This can then be used as the input to a new ISET , ISES run.

Option 15 writes out the current user-deﬁned geometric parameters in a new bparm.yyy ﬁle.

This can likewise be used to start a new ISET , ISES run.

Option 13 creates a cartesian coordinate ﬁle corresponding to the current m<sup>′</sup> − θ coordinates.

This also uses the current radius function r(m<sup>′</sup>) which was speciﬁed in the stream.xxx ﬁle.

53



<a name="br55"></a> 

The cartesian coordinates x, y, z are calculated as follows.

x = r cos(θ + ∆θ)

y = r sin(θ + ∆θ)

Z

s

Z

ꢀ

ꢁ

2

dz

dr

z =

dm<sup>′</sup> =

±

r<sup>2</sup> −

dm<sup>′</sup>

dm<sup>′</sup>

dm<sup>′</sup>

As before, the wheel axis cooincides with the z axis, and the wheel lies in the x − y plane.

The transformation assumes that z(m<sup>′</sup>) is monotonic (i.e. the ﬂow doesn’t “curl back” along

the z axis). Otherwise, the sign of the root in the integrand for z would need to be switched

after the turn point. The information needed for this decision is not present in the idat.xxx

ﬁles, and the positive root is always assumed. The result may be a mirror image of the desired

blade. The integration for z is performed around the blade contour via the trapezoidal rule.

The constant circumferential angle oﬀset ∆θ is requested from the user. This allows all the

blades to be generated in turn by incrementing ∆θ by the pitch angle.

5\.4.7 ISES Parameter Changes

Options 6, 7, and 14 are simply convenient means to change some of the contents of idat.xxx

which are not normally accessible via the usual ises.xxx input ﬁle. In particular, the quasi-3D

stream surface geometry information can be altered using option 6, which prompts for a new

ﬁle stream.yyy. This will replace the previous stream.xxx contents which were read during

ISET execution. Note that the altered idat.xxx state ﬁle must be written out with option 3

to store any new information.

5\.4.8 Inverse Design Session

Below is a sample inverse design calculation sequence, starting from the seed case xxx. Program

executions as well as option selections within EDP are shown.

% iset xxx

% ises xxx

% edp xxx

(usual analysis run)

1 Edit Cp distributions

I nitialize CPspec

M odify Cp

B lowup

(optional)

M odify Cp

(repeat as needed)

.

.

D emark inverse segment

return

2 Set design ﬂags (if necessary)

54

![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADABsDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAf/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAFcDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAAAf/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADABsDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAf/9k=)

<a name="br56"></a> 

3 Write idat.xxx

0 Quit

% edit ises.xxx

% ises xxx

% iplot xxx

% edp xxx

(add 11,12 to input lines 1 and 2)

(inverse run)

(optional)

1 Edit Cp distributions

I nitialize CPspec

M odify Cp

.

(optional)

.

return

3 Write idat.xxx

0 Quit

% ises xxx

% edp xxx

(inverse run)

5 Write new airfoil coordinate ﬁle(if satisﬁed with design)

The sequence above can be repeated as often as needed. The C<sub>p</sub> plot in EDP displays both

the current surface pressures and the speciﬁed pressures. Any diﬀerence after convergence is due

to the “error” terms which were added to attain closure at the segment endpoints. IPLOT can

of course be used anytime to examine the design in more detail.

The example above shows a Mixed-Inverse calculation. A Modal-Inverse calculation would

be done by adding 21,22,. . . to input line 1 in ises.xxx, while a Parametric-Inverse calculation

would require the addition of 40 to line 1. The latter also does not require the demarkation of

the inverse segment.

5\.4.9 Parameter-Modiﬁcation Design Session

Below is a sample design calculation sequence where the geometry parameters G(k) are repeat-

edly modiﬁed. The idea here is to implement each parameter change by simply reconverging

an existing ﬂow solution maintained in idat.xxx with ISES only, rather than starting from

scratch each time with ISET . Just reconverging can easily be an order of magnitude faster if

only modest parameter changes are made. The sequence assumes that the initial ISET setup

is done with the starting parameter ﬁle bparm.xxx, and that bspec.xxx is initially the same as

bparm.xxx.

% iset xxx

% ises xxx

(usual analysis run, with 40,40 variable,constraint ﬂags)

% bldset xxx

MODI Modify contour shape

P arameter modiﬁcation

55



<a name="br57"></a> 

Enter k, Gk:

Enter k, Gk:

Enter k, Gk: (repeat as needed)

.

.

Enter k, Gk: 0 (done)

return

PSAV Write blade parameter ﬁle: bspec.xxx

QUIT

% ises xxx

(reconverge with new parameters in bspec.xxx)

(optional)

% iplot xxx

% bldset xxx

MODI Modify contour shape

P arameter modiﬁcation

etc.

The sequence above can be repeated as often as needed. Again, IPLOT can be used anytime.

When ﬁnished, bspec.xxx contains the new blade geometry description.

5\.5 POLAR

POLAR is a driver program for ISES which conveniently sweeps through a speciﬁed parameter

range, thus generating a loss curve, turning curve, etc. Because it takes full advantage of

the quadratic convergence of the Newton method, using POLAR is more eﬃcient (and much

easier!) than running a sequence of independent cases from scratch with ISES .

In addition to the usual ISES input ﬁles, POLAR also requires a spec.xxx ﬁle, which con-

tains the sequence of operating parameters, with an optional idat.xxx nn ﬁle save ﬂag for each

parameter:

KSPEC

SPEC(1) KSAVE(1)

SPEC(2) KSAVE(2)

.

.

.

SPEC(NA) KSAVE(NA)

The ﬁrst line contains the integer KSPEC, which indicates what all the SPEC(i) values repre-

sent. The following allowed values are currently implemented, matching the global-constraint

indices GCON(.) described earlier.

56

![ref7]

<a name="br58"></a> 

KSPEC SPEC(i)

1

7

SINLin

BVR1in

BVR2in

V1ATin

MINLin

P1PTin

MOUTin

P2PTin

8

9

15

16

17

18

For example, running POLAR with the following spec.xxx ﬁle

15

0\.70

0\.72

0\.74

1

2

is equivalent to running three separate ISES cases with MINLin = 0.70, 0.72, 0.74, with the

particular MINLin in ises.xxx being ignored. The 1,2 in the second column above instructs

POLAR to write out the ﬁles idat.xxx 01 and idat.xxx 02 for the two particular points, so

that they can perhaps be examined later.

Regardless of whether any KSAVE ﬂags are speciﬁed, POLAR always overwrites the basic

idat.xxx ﬁle each time it converges a SPEC(i) parameter value.

Everytime POLAR converges on a point to the tolerances in EPS.INC, it appends the in-

tegrated parameters to polar.xxx, appends the surface pressure and boundary layer variable

distributions to polarx.xxx, and overwrites idat.xxx with the converged solution. If KSAVE

for that point is nonzero, then the additional state ﬁle idat.xxx nn is also written out, as

described above.

Note: The polarx.xxx ﬁle is not currently written for MISES 2.5, since the necessary plotting

program is not yet available.

If POLAR fails to converge on any one point, it will restart from the previously-converged

point, and subdivide the oﬀending SPEC(i) increment. If that point fails, the SPEC(i) increment

will be subdivided further. If no convergence is achieved after ﬁve subdivisions, POLAR will

terminate with a “Severe convergence problem detected” message. Failure to converge on any

one point may be due to massive separation, or a very dramatic jump in a transition location.

Occasionally, a limit cycle occurs, with a transition location oscillating back and forth. This is

often caused by inadequate resolution near the transition point. Simply restarting POLAR may

ﬁx the problem, since such a restart will begin after the oﬀending point if the preceding interval

was subdivided previously.

If after POLAR execution it is found that the parameter sweep is not complete, points can

be added to the spec.xxx ﬁle, and POLAR restarted. The new points will be automatically

appended to the save and dump ﬁles.

57

![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAC2AIgDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAf/9k=)![ref7]![ref7]![ref7]

<a name="br59"></a> 

The results in polar.xxx can be perused directly, or plotted in an organized manner using

program POLPL . This is entirely menu-driven, and is simply executed with no argument.

% polpl

The menu options allow more than one polar ﬁle to be read into POLPL , and also allow any

pair of variables to be cross-plotted, e.g. M vs p /p , ω vs M<sub>2</sub>, etc. Experimental data can

1

2

o1

also be overlaid.

Plotting software for the binary polarx.xxx dump ﬁle is currently not provided. It is expected

to be available in future MISES versions.

5\.6 BLDSET

BLDSET is a geometry manipulation program for modifying geometry ﬁles blade.xxx and/or

bparm.xxx. It is fully menu-driven and fairly self-explanatory. It can be run with the case-

extension argument,

% bldset xxx

in which case it will ﬁrst try to read the geometry in one of four ways, in the following order of

precedence:

1\. geometry parameters from an idat.xxx ﬁle

2\. x,y coordinates from an idat.xxx ﬁle

3\. geometry parameters from a bparm.xxx ﬁle

4\. x,y coordinates from a blade.xxx ﬁle

Any one of these can also read in via a menu selection at any time. BLDSET can also write

out a blade.xxx or bparm.xxx-type ﬁle, the latter being possible only if the parameters are

available.

A typical use of BLDSET in a design session is to prepare a blade.xxx or bparm.xxx ﬁle

for an analysis run — shifting the blade, changing pitch, overlaying another blade, etc. A 2-D

inviscid panel method is provided to allow a quick sanity check on the geometry.

Another use of BLDSET is to read the geometry parameters from idat.xxx, modify them

somewhat, and then write them out in a new bspec.xxx ﬁle. These modiﬁed parameters can

then be easily and rapidly “implemented” in the ﬂow solution simply by reconverging it with

ISES . This typically requires only a few Newton iterations, and hence is much faster and more

reliable than using the new parameters with ISET to start a new case.

58



<a name="br60"></a> 

6 Optimization

Optimization capabilities are no longer part of MISES 2.5 itself, but are now implemented in

the external interactive design/optimization driver LINDOP , adapted from the isolated-airfoil

MISES version MSES .

To make use of LINDOP requires the specifying some number of geometry deformation

modes as additional global unknowns (20), and specifying their corresponding ﬁxing constraints

(20). ISES sets all the geometry mode amplitudes to zero during a calculation, but it will

still calculate the sensitivities of various quantities such as S , p /p , ω, etc. to the mode

2

2

o1

displacements. These sensitivities are written out to the unformatted ﬁle sensx.xxx, which is

then read in by LINDOP .

7 Graphics

The plot library used by all the MISES programs is Xplot11, (libPlt.a), which is aimed at

driving X-terminals. A PostScript ﬁle can be generated at any time from the plot visible on

the screen. The ﬁle plotlib/Doc contains much more information on this graphics package.

8 General Hints

8\.1 Viscous solutions

When a viscous case is executed with an initial idat.xxx ﬁle from ISET , two Newton iterations

will be performed in the inviscid mode before the boundary layers (BL) are initialized and the

viscous mode is turned on. This is to allow the inviscid ﬂow to settle down from the initial

guessed ﬂowﬁeld (set up in SUBROUTINE RQINIT), so that the BLs start with a better

initial guess. In some cases two iterations may not be suﬃcient, in which case the initial BL

solution might be quite bad, greatly increasing the number of Newton iterations required for

viscous convergence. Typical examples are shocked and/or choked ﬂows, and cases with strong

streamtube contraction and rotation eﬀects on the streamline pattern. For such cases, it may be

better to perform more iterations in inviscid mode (with REYNin temporarily set to zero) before

proceeding with the viscous calculation. Alternatively, one can alter the viscous-ﬂag deﬁnition

in PROGRAM ISES (in src/ises/ises.f)

LVISC = ICOUNT.GT.2 .AND. REYNIN.GT.0.0

to increase the number of initial inviscid iterations from the current 2 iterations.

59

![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAEwDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAf/Z)

<a name="br61"></a> 

8\.2 Inverse solutions

Care must be used when running the Mixed-Inverse mode with a viscous case. It is essential that

there is no separation or near-separation anywhere within the target segment. Since the surface

pressure imposed in Mixed-Inverse is also imposed on the BL, a physically ill-posed problem

results if the BL is separated. In practice, wild changes in the geometry will result under the

separated region. Usually, the blade shape will fold up and the calculation will collapse. A

fairly simple ﬁx to this problem is to temporarily “freeze” the BLs by specifying (REYNIN = 0)

when the inverse case is converged:

1 2 5 15 11 12

1 3 4 17 11 12

.

.

0\.0e6 ...

| REYNin

The BLs case then be “unfrozen” and the case reconverged in the usual analysis mode:

1 2 5 15 ! 11 12

1 3 4 17 ! 11 12

.

.

1\.0e6 ...

| REYNin

The resulting viscous p(s) will change slightly from the speciﬁed p<sub>spec</sub>(s) in the inverse calcula-

tion, but this is usually minor, and can be iterated if desired. Note that the use of “!” in the

ﬁrst two lines is convenient for quickly adding and removing the global variable and constraint

ﬂags.

The Modal-Inverse option generally does not have the problem with ill-posedness in separated

ﬂow, and can be usually be safely used in all situations.

The well-posedness of the Parametric-Inverse method depends a great deal on the action of

the user-deﬁned parameters. It is important that all the declared parameters be reasonably

independent in their inﬂuence on the blade shape, and that they all have an aerodynamically

signiﬁcant inﬂuence on the surface pressure distribution. One cannot least-squares ﬁt a param-

eter if that parameter has no inﬂuence on the ﬂow! If such a parameter is declared without

any constraints, the result will be a nearly-singular least-squares matrix, producing enormous

Newton changes and a certain solution failure.

8\.3 Grid resolution

Compared to most Euler solvers, ISES is usually quite insensitive to grid density. It behaves

more like a potential solver in this regard, especially if ISMOM=3 or ISMOM=4 are used. In

60



<a name="br62"></a> 

any case, it is a good idea to check that the leading edges are reasonably well resolved, And

that no spurious losses (or gains!) are being generated there. The “Streamtube plots” menu

of IPLOT can be used to plot the stagnation pressure variations along any streamtube. This

should be piecewise-constant, with sudden monotonic drops through the shocks.

Any separation bubble present in the ﬂow must be well-resolved. The default grid is usually

adequate for most cases, but maybe not if the bubble is close to the leading edge and very small

in streamwise extent. Moderate Reynolds numbers (1-3 million, say) require the ﬁnest grid,

since the bubbles are then still important, but very small. Fortunately, streamwise grid spacing

is ”cheap”, increasing the solution time only linearly, so it may be simplest to increase the grid

point number parameter N in ISET to 100 or more. Inadequate bubble resolution often results

in a ”ragged” or ”scalloped” loss vs. incidence curve, so this is easy to spot.

8\.4 Execution times

Using ISES or its variants requires substantial computational resources. A single point solution

for a typical grid size (160 × 24) will require anywhere between 3 Newton iterations for a

subcritical inviscid case, to ∼ 12 iterations for shocked transonic or low-Re viscous cases with

separated regions. Each iteration represents 2–10 CPU seconds on a RISC workstation. For

a given grid size, the oﬀset-periodic inlet grid topology roughly doubles the required time for

a solution. But since 2/3 as many streamlines are typically required on an oﬀset grid for the

same accuracy, the net CPU time required in practice is similar.

When a sequence of solutions is needed, the CPU requirements per solution can drop substan-

tially due to the quadratic convergence property of the Newton method used by ISES . Once a

converged solution is obtained, convergence to a “nearby” solution after a small inlet angle or

inlet Mach (say) change is quite rapid, requiring only 2 or 3 additional iterations. When inverse

design calculations are performed, a minor prescribed-pressure change will typically require only

about 2 Newton iterations to converge the new geometry, making interactive design even on a

low-end RISC workstation eﬀective.

While the Newton solution method is very eﬃcient for converging small parameter tweaks,

it intensely dislikes large changes. Trying to reconverge a solution after a drastic parameter

change, such as going from S = 0.5 to S = 1.0, is deﬁnitely not a good idea if a transonic

1

1

or supersonic viscous case is involved. Clearly, runs should be sequenced so that such large

changes are avoided.

61

![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADADoDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAf/9k=)


<a name="br63"></a> 

MISES Roadmap

blade.xxx

or

ꢃ

\-

\-

¨\*

BLDSET

bspec.xxx

bparm.xxx

gridpar.xxx

¡ꢁ

H

¨

¡

H

?

Hj

¡

¡

\-

¡

stream.xxx

loss.xxx

ISET

¡

¡

¡

¨\*

¨

¡

¡

spec.xxx

bparm.xxx

blade.xxx

¡

¡

¡

¡

@I

¡

¡

@

¡

¡

¡ꢂ ©

?

@

polarx.xxx

?

idat.xxx

6

¡

¨

ꢃ

ꢃ

\-

ꢄ

ꢃ

\-

EDP

POLAR

H

Hj

ꢄ

polar.xxx

?

6

ꢄ

ꢄ

ꢄ

ꢄ

ꢄ

ꢄ+

?

ISES

6

\-

ref. data

IPLOT

POLPL

?

ꢀ

¡ꢁ

¡

Xplot11

Graphics

¡

¡

bspec.xxx

suct.xxx

¡

ises.xxx

Program

¡

¡

Data ﬁle

Data ﬁle (optional)

![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAMMDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAAAAAAAAAAAAAAAAAf/Z)![ref22]![ref22]![ref22]![ref22]![](data:image/jpeg;base64,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)![](data:image/jpeg;base64,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)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAxALoDASIAAhEBAxEB/8QAGgABAQACAwAAAAAAAAAAAAAAAAgHCQUGCv/EADAQAAEDBAEBBQUJAAAAAAAAAAABAgUDBAYIBzgJESF3eDFBc7a3EhMXIjNEgbLC/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/APfwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAw/sHNSeOcD805BCSF1ETUFxPyLMxEtY1Pub2Lk4vEJe+spC0qp4surO5oUrm3cnilamzuVPakX6Y57mue6e6oZ1lk5XynKc01r4KyzJcnmpdVmcjyDI+L8WmJmelldRqKsnLyN5cyF+qveq3VxVVXuXxWu9oGtdrdz+jkRzfwT5WcrV8WuVmCTz2/aT3ojmovcvvRCIOz8g4W60L0kua0VHvrXGomtleq9bWl3vq1eGsMqVHL+X2ue5V/kDaoAAAAAAAAAAAAAAAAAAAAAAAAAAAAAwXs/03bAeSPLHyFPkZ9nr0CaO+j/Wj6L4UWZs/wBN2wHkjyx8hT5GfZ69Amjvo/1o+i+FAbQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGC9n+m7YDyR5Y+Qp8jPs9egTR30f60fRfCi2NjYuVm+AOboWCj68tNzPEfJEPERVpTfVvJOTlMOmLCwj7KlTa59a8vLu4o29tRa1z61aoykxquehIukeIZZgml+omD5jByeL5dhusHAWKZTjU3YV7Cax7I8e4pxOInIKXsa6Mr2UpEydndWF/aVmMq213b1qNRrXsciBsiAAAAAAAAAAAAAAAAAAAAAAAAAAAAAcdIftvj/4cdcrfq1fiP8A7KAB/9k=)![](data:image/jpeg;base64,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)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCABEAF8DASIAAhEBAxEB/8QAGwABAQACAwEAAAAAAAAAAAAAAAECBQYICQr/xAArEAACAQEFCAMAAwEAAAAAAAAAAQIDBAUGBwgJERkhUViRmTFB1hIicWH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A+/gAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAE/lHn/aPJtPmuTXyn/wBX2RSi0mpRafNNNNNdU0+YGQJvXVeUN66rygKCb11XlDeuq8oCgm9dV5Q3rqvKAoJvXVeUN66rygKDFyivmUVzS5tLm3uS+flvkl9sKcH8Si/8kn97vp9U1/q3AeZ2ZGzMymzZx3ijMK+s/tc+G70xtettv+8bly91oZ8YDwddlqrOM5WXDWE8OYssVy4duuLm1Tu267NQslOKjGNNKKOBXVshMka93WKrLU3tG4yqWenJxp6/tScIJtfEIRxqoxivpJbkABsOEBkf3ObR72A6lP2o4QGR/c5tHvYDqU/agAOEBkf3ObR72A6lP2o4QGR/c5tHvYDqU/agAOEBkf3ObR72A6lP2o4QGR/c5tHvYDqU/agAOEBkf3ObR72A6lP2o4QGR/c5tHvYDqU/agAaS/NkpkxddKxzs2pnaKzdpvGw2OorTr51IWmKpWm1UaM5QjVxo1CtGNSUqVWO6dKoo1INSimdqtMWinL/AEv31ia+8IZqanMwLRiq67Hdlssee+o3NHOu6ruo2C217TTr4fuzH1+3vZLittadWUbZbLup0a9qpKNKtOUIpAAf/9k=)![ref23]![ref24]![ref25]![ref26]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref28]![ref28]![ref28]![ref28]![ref23]![ref24]![ref25]![ref26]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref28]![ref28]![ref28]![ref28]![ref29]![ref24]![ref30]![ref26]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref28]![ref28]![ref28]![ref28]![ref23]![ref24]![ref25]![ref26]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref28]![ref28]![ref28]![ref28]![ref29]![ref24]![ref30]![ref26]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref28]![ref28]![ref28]![ref28]![ref23]![ref24]![ref25]![ref26]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref28]![ref28]![ref28]![ref28]![ref31]![ref31]![ref31]![ref31]![ref23]![ref24]![ref25]![ref26]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref28]![ref28]![ref28]![ref28]![ref23]![ref24]![ref25]![ref26]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref28]![ref28]![ref28]![ref28]![ref23]![ref24]![ref25]![ref26]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref27]![ref28]![ref28]![ref28]![ref28]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAfADYDASIAAhEBAxEB/8QAGgABAAIDAQAAAAAAAAAAAAAAAAQJBQgKA//EADMQAAECBAQDBAoDAQAAAAAAAAECAwQFBwgABhESCRQhEzJBURUYIiMxWZeh0ddCU3KB/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AO/jDDDAMMMMAxhI1s82l9UU+w2ygEJbIT2jytEob111caVu3KYATudAcKwE6YzeIERCsvRCHXE7lstqcaOp926Dt7VA/i4UFTRWNFFpa0a7VEYCkSyC+Wod3d9BVEpj8kU6grZboJREUrZzI9PJHE54oVftNKBwVRlByWSxKJvOss5adKm+xWZUxNHZSiJjm2BGOsaNcFdxZvNjUk6hVEOJIs/6TxeqgoH/AAJ6AeAwwHVFzqAhKyhWqkKWUJUhaklOmqDsUoFXXrtJA8SBoT7NPodQFpB2knaeh1HgfjqNfJQCh4gY1FuQodWWtAygqlF3NUrWGZCmdRM5cphkejmbRnUTL0dyaZmmq+Qc7Oy4SEQcTyqZQJcIv0o/zwijDwpY1obsUvMQkIhOLhdZCMJ9lEKzRezZ5tjTp2aHoygr0U6E/wBj7rjivipaj1wFqu8eR+35w3jyP2/OKrfUWvV+bzdj9E7Lf0Bh6i16vzebsfonZb+gMBalvHkft+cRlvIDykHcCIdSu6TqkuJRqNoPwUdNve8QCOuKtvUWvV+bzdj9E7Lf0BiE9Yrekt4MHi6XTrUU8woqorZ0iKLbeiEhtbVBUQw96Uey4NCkEkE9cBVdwVNVXnTAgEJbojxIm9xHRRVxeKgOA9O5oDoUL2rPeCSnrhi0CyPhnxNnFd5VUKV1gnFT5BLqC1fp1mWOzzK5LLqgZjqTWi56JuRn+boxnJsnkWSISUCImkzlrUFKZbCOodcYV2OxLisMB//Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAbADEDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAAAAAAAf/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAIAAUDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAeEAACAgICAwAAAAAAAAAAAAAFBgQHAAMBAhETJP/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwC4Wv7dmOsx9iT6usmr9SW8mU+ASsoauChNmRRvp2cWBWm5ZZ2uUVRDHbfz2GEWWGqHt3j6l2HzjGMD/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAIAAUDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAf/Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAIAAUDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAfEAABBAMAAwEAAAAAAAAAAAAGAQQFBwADCAITFiH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8ArEqoo6QI+weyXxUHW8I0PFRPP4nSOsy8Bz5MlJR7RbKXAZ1q3iyeVeKPT+1+A++RnGUBKSSs2a7Ydujb9Yxgf//Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAIAAUDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAgEAAABQQDAQAAAAAAAAAAAAACAwQFBgABBwgSExUU/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AK+WM7YRt3G2iUyF0kzfgVZjrWYrCpKxWnVxi8qRkZhtmMUcbgqe9uWdqiAe+ceUWFeLy7lXH8pnFSlB/9k=)![ref7]![ref7]![ref7]![ref7]![ref7]![ref7]![ref7]![ref7]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAHAAMDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAf/Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAHAAMDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAf/Z)![ref32]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEACMDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAUGCAr/xAAnEAAABQIFBAMBAAAAAAAAAAABAgMFBgAEBwgJERIZIViYFUHWcf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDrsnuQ5DEaZyKanziZ9IQeWOV69KRfDrMk4xWGMCivEwt8XYEo5dJs7UTlsnZEXWKUpShzHaqW1aatqq3Wag579SsgnRKbgnm0dSJl337EKEU7FD6Dcf7SlBIdNG088dS722dvylOmjaeeOpd7bO35SlKDc2H+GpMP4VGoUSdYkzEkaabZpLKZ9LFpLM30LYol+Rkj+pa2x3Z1X33ub01ukZYwAIkClKUH/9k=)![ref32]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEACMDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAUGCAr/xAAnEAAABQIFBAMBAAAAAAAAAAABAgMFBgAEBwgJERIZIViYFUHWcf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDrsnuQ5DEaZyKanziZ9IQeWOV69KRfDrMk4xWGMCivEwt8XYEo5dJs7UTlsnZEXWKUpShzHaqW1aatqq3Wag579SsgnRKbgnm0dSJl337EKEU7FD6Dcf7SlBIdNG088dS722dvylOmjaeeOpd7bO35SlKDc2H+GpMP4VGoUSdYkzEkaabZpLKZ9LFpLM30LYol+Rkj+pa2x3Z1X33ub01ukZYwAIkClKUH/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCABQAAQDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAQGBwr/xAAjEAAABQMFAQEBAAAAAAAAAAACAwQGCAABBQcJM3m1ERJB/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AOuiCmoT31GKmCY9HFnXeYyp3yP08bI8qtEcFts5r5FulYJo4kQgi+YXBgVHgQADYAQhPM+AD/VVHbYSFDTTuBb9lFkblMsk5RZQrgLAUVlWrYAQht9+fPt6UEnbV4J59l8tvWalKbavBPPsvlt6zUpQNtXgnn2Xy29ZqUptq8E8+y+W3rNSlA21eCefZfLb1mpSm2rwTz7L5bes1KUDbV4J59l8tvWalK06Fmjj40UTyotqInx+KO1Zmlr9rczi0GXQ5aypg6hLsCobK5WNIMQUS5SXj1N1GOP/ACpSXCGxwbXHalB//9k=)![ref32]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCABQAAQDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAQGBwr/xAAjEAAABQMFAQEBAAAAAAAAAAACAwQGCAABBQcJM3m1ERJB/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AOuiCmoT31GKmCY9HFnXeYyp3yP08bI8qtEcFts5r5FulYJo4kQgi+YXBgVHgQADYAQhPM+AD/VVHbYSFDTTuBb9lFkblMsk5RZQrgLAUVlWrYAQht9+fPt6UEnbV4J59l8tvWalKbavBPPsvlt6zUpQNtXgnn2Xy29ZqUptq8E8+y+W3rNSlA21eCefZfLb1mpSm2rwTz7L5bes1KUDbV4J59l8tvWalK06Fmjj40UTyotqInx+KO1Zmlr9rczi0GXQ5aypg6hLsCobK5WNIMQUS5SXj1N1GOP/ACpSXCGxwbXHalB//9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAB4DASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAQGBwr/xAAkEAAABQQCAwADAAAAAAAAAAABAgMEBQYHCBEAEgkTIRQxQf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDrsr7A5C41Y1JWp8xM9KIPVki9mlKWt1klI0rRkCor0MMfTEAlTjpOHiSd9JsiLrAQpSh3HXKXE+NVqrGslRzv8lZBO3TMJE8tJUiZdh+iFClNFKH8D7rjjgRp7xwoxkedyjnZ5JVje1BMU3WWEquiYqipCD2TGlS7Eu+xR38MAD91rm82HxAb2MqWXqNLJXLy7RpaFcQgwl8L4vbh07HdpNk+NLxUW4hI8rOcEzIGwyIKnN+G6eoevTgTFccD/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAeAAQDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAMGBwr/xAAkEAABAwMEAgMBAAAAAAAAAAACAQMFBAYIAAcREgkVEyJBFP/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDrnwU3CvfcZrMFy9LinbvcsrPDI/by2Tla1XhtuzrXkLdbgrSiSIT4hYQKp8KABQBEXnOGx/Wql42KRo6bO4E7tNseSnLKnabaJQbBpqVtVAERTnjjldNBN43HQpmc8BeLornkqyzeBOFLs05K2qoF9EJE5RF+q8En6iaa2/EnZi4djRyZauGWhZRd2cvN8d74T0QVSpH29uHWwj8XFyy1tHRL7yjCOdGR/mSqo+xt/BW1CdlFoP/Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCABQAAQDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAQGBwr/xAAjEAAABQMFAQEBAAAAAAAAAAACAwQGCAABBQcJM3m1ERJB/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AOuiCmoT31GKmCY9HFnXeYyp3yP08bI8qtEcFts5r5FulYJo4kQgi+YXBgVHgQADYAQhPM+AD/VVHbYSFDTTuBb9lFkblMsk5RZQrgLAUVlWrYAQht9+fPt6UEnbV4J59l8tvWalKbavBPPsvlt6zUpQNtXgnn2Xy29ZqUptq8E8+y+W3rNSlA21eCefZfLb1mpSm2rwTz7L5bes1KUDbV4J59l8tvWalK06Fmjj40UTyotqInx+KO1Zmlr9rczi0GXQ5aypg6hLsCobK5WNIMQUS5SXj1N1GOP/ACpSXCGxwbXHalB//9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEACMDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAUGCAr/xAAnEAAABQIFBAMBAAAAAAAAAAABAgMFBgAEBwgJERIZIViYFUHWcf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDrsnuQ5DEaZyKanziZ9IQeWOV69KRfDrMk4xWGMCivEwt8XYEo5dJs7UTlsnZEXWKUpShzHaqW1aatqq3Wag579SsgnRKbgnm0dSJl337EKEU7FD6Dcf7SlBIdNG088dS722dvylOmjaeeOpd7bO35SlKDc2H+GpMP4VGoUSdYkzEkaabZpLKZ9LFpLM30LYol+Rkj+pa2x3Z1X33ub01ukZYwAIkClKUH/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAGcDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAUGBwr/xAAnEAAABAQGAwADAAAAAAAAAAAAAgQFAQMGCQcIGViY1hESIRhBcf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDrsr3IdIxGrOoq1PnEz6UQerHJa9TKXw6zJONK0YwTJvqaLfS7BKpxVLZ2ont4loiT5xSlKWHvHwKW1W1Us1uRzI577lZInklN6S82jqSWXz5+ELClPhYfqHmP9AAEhpopN+Ny7ls7dUDTRSb8bl3LZ26oAAGmik343LuWzt1QNNFJvxuXctnbqgAAaaKTfjcu5bO3VA00Um/G5dy2duqAABpopN+Ny7ls7dUDTRSb8bl3LZ26oAAGmik343LuWzt1QNNFJvyuXctnbqgAA2rADK9+Pri+PcrMVmqxnmVMiI1TG7MDjMrxNaGWCBadQRwp1CpZmyDW6KIQglVrCzJsZ6SEJMSFhD2AAAf/2Q==)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAEcDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAUGBwr/xAAnEAAABQIGAwADAQAAAAAAAAABAgMEBQAGBwgJGViYERLWGCExQf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDrrv7IchiNeVxXqfOHnzsg92SL2aUtbDrMlI2rZkCoqBTDH2vAJW46Th4ont4TZEXWAhSlD3HxVMidNVqrGslRz36lZBO3TMJE82kqRMvkP4QoWp4KUP8AA/filKCQ20WnPHUu7bS3ylNtFpzx1Lu20t8pSlA20WnPHUu7bS3ylNtFpzx1Lu20t8pSlA20WnPHUu7bS3ylNtFpzy1Lu20t8pSlBtWAGV78fZGcm0sxWarGdS5mRIpSOzA4zO8TYiFBg9O4JIW6xcw0YEXKOAAGrt4VRUV2gAiJCgHtSlKD/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEACwDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAUGCAr/xAAnEAAABQMDBAIDAAAAAAAAAAABAgMEBQAGBwgREgkZWJgh1hZBcf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDrsv3Qchka87ivU+sTXpZB7skns0pa+OtSUjatmQCivEwx9rwCVuOk4eKJy2TZEXWKUpShzHaqXFdNVqrHM1B139SsgnRKbgnq0lSJl33+CFC1PgofoNx/tKUEh20Wnnj1LvbaW+qU7aLTzx6l3ttLfVKUoHbRaeePUu9tpb6pW0cH4gJg+w29iJ5Ly9lkreTkpP8AMM4XyvkO/FhkVCKCwcXI4Yx6isYx4cI1oLcAaJmOQDnA24KUH//Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEACwDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAUGCAr/xAAnEAAABQMDBAIDAAAAAAAAAAABAgMEBQAGBwgREgkZWJgh1hZBcf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDrsv3Qchka87ivU+sTXpZB7skns0pa+OtSUjatmQCivEwx9rwCVuOk4eKJy2TZEXWKUpShzHaqXFdNVqrHM1B139SsgnRKbgnq0lSJl33+CFC1PgofoNx/tKUEh20Wnnj1LvbaW+qU7aLTzx6l3ttLfVKUoHbRaeePUu9tpb6pW0cH4gJg+w29iJ5Ly9lkreTkpP8AMM4XyvkO/FhkVCKCwcXI4Yx6isYx4cI1oLcAaJmOQDnA24KUH//Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAVAAQDASIAAhEBAxEB/8QAFwABAAMAAAAAAAAAAAAAAAAAAAYICv/EACgQAAAFAwMBCQAAAAAAAAAAAAECAwQFBgcIABESCRMUFRYXIiM0Qf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDXRgpcKt7jJZgqVpUU9V6lFZ35H28pk0q9FYlN0fS8hTqcFSMQYwH2hIMjpdOPKAJlKVZTZMv61A+m8mAtM5iiqkyFPqQZWJ93TE4op8ZOlfrmApOaA7/GfiXkG/tDTQWrsNj/AB+PvrQ3jqkd1QS8WQFzb+OfEoxq08uvLluY1y4pWPBJdwDiMhhjSkZu1OwVXBY4nbJcQAWmmg//2Q==)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAQAA0DASIAAhEBAxEB/8QAFgABAQEAAAAAAAAAAAAAAAAACAcK/8QAHRAAAQUBAQEBAAAAAAAAAAAABQIDBAYHCAEUE//EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDUtsHQFph9yU/CgvenPWdMWkBC9hcnWvDptt122kZwE9JGnReiI0UHFijJ0qI2bjjF1Z734gk2F9fqpiZDF242vu42ij22HqmiZn0A7UL5Y6cK27JR8iiibkRrhScIu4MrmUiXZF0MrnFwhFKE9FTb7Ik16DWaS8O8leDmCRpc7uOjdj9SasMrV+13jqNm3PVbH4/UhEL3RSc8xXdQ80PTubTswwzHVcqMZjVAeez6e1XoFkatEa1LuQ+VSIgY+3OHRkXznKoqSN6uE+ulLY+5F7QLwjPQf6P2Qk8t62EYBc9GkDHlL9dq3qSTq0V1Y9DqGXPFNID/2Q==)

[ref1]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADABcDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAf/9k=
[ref2]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAEkDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAf/Z
[ref3]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADABUDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAf/9k=
[ref4]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAAoDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAf/Z
[ref5]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADABgDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAf/9k=
[ref6]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADADYDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAf/9k=
[ref7]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAAcDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAf/Z
[ref8]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADACQDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAf/Z
[ref9]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADADEDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAf/9k=
[ref10]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAAgDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAf/Z
[ref11]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAB8DASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAf/9k=
[ref12]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADABEDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAf/9k=
[ref13]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADABIDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAf/9k=
[ref14]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADABMDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAf/9k=
[ref15]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADACMDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAf/Z
[ref16]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADADIDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAf/9k=
[ref17]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADACADASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAf/9k=
[ref18]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAD0DASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAf/9k=
[ref19]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADADwDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAf/9k=
[ref20]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADACUDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAf/Z
[ref21]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADADcDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAf/9k=
[ref22]: data:image/jpeg;base64,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
[ref23]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAJAAYDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAf/Z
[ref24]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAKAAYDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAf/Z
[ref25]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAJAAYDASIAAhEBAxEB/8QAFgABAQEAAAAAAAAAAAAAAAAAAAkK/8QAIhAAAQQCAQQDAAAAAAAAAAAAAwQFBgcBAgAICRYXNIOy/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/ANYNzPF+Iu4FU7XEfamaVU9Ht4r5NqzBkhK6zZyS6qATQwjwpAHMd0m+8SVzbLMBUp0eyMvkJEYN0QnAg3KJunyC/R+CccD/2Q==
[ref26]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAKAAYDASIAAhEBAxEB/8QAFgABAQEAAAAAAAAAAAAAAAAAAAgK/8QAHhAAAgIDAAMBAAAAAAAAAAAABAYDBQECBwAIFBP/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A1sdl9jemKHvVy/gFRGvSc8bvVHs/YLnUyr2kuNnNC6/wxIXZh7TUvGIazah6Ex/QBkPf9yciEYKj+bMUrymGBJTD+2r7ccorBrWBy5rWwGcugqiWEJdsWtMsrCgEuphN7IalPsamrPMqoSdASja2vKng3nDGkjeB/9k=
[ref27]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAADAAkDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAf/Z
[ref28]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAJAAMDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAf/Z
[ref29]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAJAAYDASIAAhEBAxEB/8QAFgABAQEAAAAAAAAAAAAAAAAAAAcK/8QAIRAAAQQBAwUAAAAAAAAAAAAAAwIEBQcGAAEIEhQ2hLT/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A2i4nzH4uZ/btg0DidvYjlNwVWSVXYlfxi35JzFSREoxhZQso3OxC1H20pJsmJ9wPHG6TuRpT1p3UtLVBk/JpD2/oBpoP/9k=
[ref30]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAJAAYDASIAAhEBAxEB/8QAFgABAQEAAAAAAAAAAAAAAAAAAAkK/8QAIBAAAAYDAAMBAAAAAAAAAAAAAQMEBQYHAAIICTZytP/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDSDKewC7B8oyKhqZupVJ49V/KV8mXZXUdVOQMUQttmvPn9jiRkjSmoi0mssSRxzmiFGaSoPEUKp1DQQ022xld2r2Z5+T/0E4wP/9k=
[ref31]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAxAF8DASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAf/9k=
[ref32]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAjAAQDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAQGCAr/xAAnEAABAgUDAwUBAAAAAAAAAAACAwYAAQQFCAcJERIzeRQhMUG1Uf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDrnwV1De+oqWYKj1cV+d6jKzvyP08bJXWsNYW0zmvcW6lYmjaSIC4stjCqXCgAZAIiupwA/aKltsUqRU2dwD1JJoblOWSCSaRdIAkldWpIBlLifxKf9hASdtXsZ5+S/Lb9ZqQhtq9jPPyX5bfrNSEBvhh6dMbTwXeLIa1nbAvZ9uTUN2StFKNLK/vd0K06jhc9y6Zz9Rd7wpS0511VPgliRCc5e0IQgP/Z
